#pragma once
#include <fixi/backend.hpp>
#include <fixi/buckets.hpp>
#include <fixi/defines.hpp>
#include <fixi/utils.hpp>
#include <iostream>
namespace Fixi
{

class Rattle
{
    const std::vector<FixedBondLengthPair> pairs{};
    const std::vector<std::vector<FixedBondLengthPair>> buckets{};

    std::vector<std::vector<double>> pairwise_hg_ij{};
    std::vector<std::vector<Vector3>> pairwise_constraint_forces{};

    /** The scaled virial. Needs to be multiplied by timestep^2 to get the virial */
    Matrix3 scaled_virial{};

    int iteration{ 0 };

    Matrix3 compute_scaled_virial(
        const Eigen::Ref<Vectorfield> unadjusted_positions, const Vector3 & cell_lengths,
        const std::array<bool, 3> & pbc )
    {
        Matrix3 stress_tensor = Matrix3::Zero();

        // In order not to create any race conditions we iterate over the buckets serially
        for( int idx_bucket = 0; idx_bucket < buckets.size(); idx_bucket++ )
        {
            const auto & bucket = buckets[idx_bucket];

            auto cb = [&]( int idx_pair ) {
                const auto pair = bucket[idx_pair];
                const int i     = pair.i;
                const int j     = pair.j;

                // `s` is the current approximation for the vector displacement between atoms i and j
                const Vector3 rij
                    = mic( unadjusted_positions.row( i ) - unadjusted_positions.row( j ), cell_lengths, pbc );
                const double hg   = pairwise_hg_ij[idx_bucket][idx_pair];
                const Vector3 cij = -2.0 * hg * rij;

                // Save the pairwise constraint forces
                pairwise_constraint_forces[idx_bucket][idx_pair] = cij;

                const Matrix3 stress_tensor_ij = rij.transpose() * cij;
                return stress_tensor_ij;
            };

            auto red = [&]( auto & lhs, const auto & rhs ) { lhs += rhs; };

            stress_tensor += Backend::transform_reduce<Matrix3>( bucket.size(), cb, red, Matrix3::Zero() );
        }
        return stress_tensor;
    }

public:
    int maxiter{};
    double tolerance{};

    const std::vector<FixedBondLengthPair> & get_pairs() const
    {
        return pairs;
    }

    const Vectorfield get_constraint_forces( int n_atoms ) const
    {
        Vectorfield forces( n_atoms, 3 );
        forces.setZero();

        for( int idx_bucket = 0; idx_bucket < buckets.size(); idx_bucket++ )
        {
            const auto & bucket = buckets[idx_bucket];

#pragma omp parallel for
            for( int idx_pair = 0; idx_pair < bucket.size(); idx_pair++ )
            {
                const auto & pair   = bucket[idx_pair];
                const int i         = pair.i;
                const int j         = pair.j;
                const Vector3 & cij = pairwise_constraint_forces[idx_bucket][idx_pair];

                forces.row( i ) += cij;
                forces.row( j ) -= cij;
            }
        }
        return forces;
    }

    Matrix3 get_scaled_virial() const
    {
        return scaled_virial;
    }

    int get_iteration() const
    {
        return iteration;
    }

    Rattle( int maxiter, double tolerance, const std::vector<FixedBondLengthPair> & pairs )
            : pairs( pairs ), buckets( bucket_sorter( pairs ) ), maxiter( maxiter ), tolerance( tolerance )
    {
        pairwise_hg_ij.resize( buckets.size() );
        pairwise_constraint_forces.resize( buckets.size() );

        for( int idx_bucket = 0; idx_bucket < buckets.size(); idx_bucket++ )
        {
            pairwise_hg_ij[idx_bucket].resize( buckets[idx_bucket].size() );
            pairwise_constraint_forces[idx_bucket].resize( buckets[idx_bucket].size() );
        }
    }

    double adjust_positions(
        const Eigen::Ref<Vectorfield> unadjusted_positions, Eigen::Ref<Vectorfield> adjusted_positions,
        const Eigen::Ref<Scalarfield> masses, const Vector3 & cell_lengths, const std::array<bool, 3> & pbc )
    {
        iteration = 0;
        double hij_max{ 0.0 };

        // reset the pairwise_hg_ij to zero
        for( int idx_bucket = 0; idx_bucket < buckets.size(); idx_bucket++ )
        {
            const auto & bucket      = buckets[idx_bucket];
            const int n_pairs_bucket = bucket.size();
#pragma omp parallel for
            for( int idx_pair = 0; idx_pair < n_pairs_bucket; idx_pair++ )
            {
                pairwise_hg_ij[idx_bucket][idx_pair] = 0;
            }
        }

        do
        {
            hij_max = 0.0;
            // In order not to create any race conditions we iterate over the buckets serially
            for( int idx_bucket = 0; idx_bucket < buckets.size(); idx_bucket++ )
            {
                const auto & bucket = buckets[idx_bucket];

                auto cb = [&]( int idx_pair ) {
                    const auto & pair = bucket[idx_pair];

                    const int i      = pair.i;
                    const int j      = pair.j;
                    const double mi  = masses[i];
                    const double mj  = masses[j];
                    const double dij = pair.dij; // desired constrained bond length

                    // `s` is the current approximation for the vector displacement between atoms i and j
                    const Vector3 rij
                        = mic( unadjusted_positions.row( i ) - unadjusted_positions.row( j ), cell_lengths, pbc );
                    const Vector3 s
                        = mic( adjusted_positions.row( i ) - adjusted_positions.row( j ), cell_lengths, pbc );
                    const double s2   = s.squaredNorm();
                    const double dij2 = dij * dij;

                    const double hij = abs( s2 - dij2 ); // holonomic constraint

                    // if the constraint is violated, change adjusted_positions
                    // this is achieved by solving
                    //   | (adjusted_positions[i] - h g r_ij / mi) - (adjusted_positions[j] + h g r_ij / mj) |^2 = dij^2
                    // for g.
                    // When neglecting g^2 terms in the derivation, we get (also h cancels so we don't need the step size)
                    const double hg = ( s2 - dij2 ) / ( 2.0 * ( s.dot( rij ) ) * ( 1.0 / mi + 1.0 / mj ) );

                    adjusted_positions.row( i ) -= hg * rij / mi;
                    adjusted_positions.row( j ) += hg * rij / mj;

                    // accumulate the total lagrange multiplier
                    pairwise_hg_ij[idx_bucket][idx_pair] += hg;

                    return hij;
                };

                auto red = [&]( auto & hij_max_bucket, const auto & hij_max_thread ) {
                    hij_max_bucket = std::max( hij_max_bucket, hij_max_thread );
                };

                const double hij_max_bucket = Backend::transform_reduce<double>( bucket.size(), cb, red );
                hij_max                     = std::max( hij_max_bucket, hij_max );
            }
            iteration++;
        } while( ( hij_max > tolerance ) && ( iteration < maxiter ) );

        // Before returning, we also compute the scaled virial (aka the stress tensor)
        scaled_virial = compute_scaled_virial( unadjusted_positions, cell_lengths, pbc );

        return hij_max;
    }

    double adjust_velocities(
        const Eigen::Ref<Vectorfield> positions, Eigen::Ref<Vectorfield> adjusted_velocities,
        const Eigen::Ref<Scalarfield> masses, const Vector3 & cell_lengths, const std::array<bool, 3> & pbc )
    {
        iteration = 0;
        double hij_max{ 0.0 };
        do
        {
            hij_max = 0.0;

            // In order not to create any race conditions we iterate over the buckets serially
            for( int idx_bucket = 0; idx_bucket < buckets.size(); idx_bucket++ )
            {
                const auto & bucket = buckets[idx_bucket];
                auto cb             = [&]( int idx_pair ) {
                    const auto & pair = bucket[idx_pair];

                    const int i       = pair.i;
                    const int j       = pair.j;
                    const double mi   = masses[i];
                    const double mj   = masses[j];
                    const double dij  = pair.dij; // desired constrained bond length
                    const double dij2 = dij * dij;

                    const Vector3 rij = mic( positions.row( i ) - positions.row( j ), cell_lengths, pbc );
                    const Vector3 vij = adjusted_velocities.row( i ) - adjusted_velocities.row( j );

                    const double hij = abs( rij.dot( vij ) );
                    const double k   = rij.dot( vij ) / ( dij2 * ( 1.0 / mi + 1.0 / mj ) );

                    adjusted_velocities.row( i ) -= k * rij / mi;
                    adjusted_velocities.row( j ) += k * rij / mj;
                    return hij;
                };

                auto red = [&]( auto & hij_max_bucket, const auto & hij_max_thread ) {
                    hij_max_bucket = std::max( hij_max_bucket, hij_max_thread );
                };

                const double hij_max_bucket = Backend::transform_reduce<double>( bucket.size(), cb, red );
                hij_max                     = std::max( hij_max_bucket, hij_max );
            }
            iteration++;
        } while( ( hij_max > tolerance ) && ( iteration < maxiter ) );

        return hij_max;
    }
};

} // namespace Fixi