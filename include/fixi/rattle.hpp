#pragma once
#include <Eigen/Core>
#include <fixi/backend.hpp>
#include <fixi/buckets.hpp>
#include <fixi/defines.hpp>
#include <fixi/utils.hpp>
namespace Fixi
{

class Rattle
{
    const std::vector<FixedBondLengthPair> pairs{};
    const std::vector<std::vector<FixedBondLengthPair>> buckets{};
    std::vector<std::vector<double>> pairwise_hg_ij{};

    int iteration{ 0 };

public:
    int maxiter{};
    double tolerance{};

    const std::vector<FixedBondLengthPair> & get_pairs() const
    {
        return pairs;
    }

    int get_iteration() const
    {
        return iteration;
    }

    Rattle( int maxiter, double tolerance, const std::vector<FixedBondLengthPair> & pairs )
            : pairs( pairs ), buckets( bucket_sorter( pairs ) ), maxiter( maxiter ), tolerance( tolerance )
    {
        pairwise_hg_ij.resize( buckets.size() );
        for( int idx_bucket = 0; idx_bucket < buckets.size(); idx_bucket++ )
        {
            pairwise_hg_ij[idx_bucket].resize( buckets[idx_bucket].size() );
        }
    }

    double adjust_positions(
        const Eigen::Ref<Vectorfield> unadjusted_positions, Eigen::Ref<Vectorfield> adjusted_positions,
        const Eigen::Ref<Scalarfield> masses, const Vector3 & cell_lengths, const std::array<bool, 3> & pbc )
    {
        iteration = 0;
        double hij_max{ 0.0 };
        do
        {
            hij_max = 0.0;

            int idx_bucket = 0;
            // In order not to create any race conditions we iterate over the buckets serially
            for( const auto & bucket : buckets )
            {
                const int n_pairs_bucket = bucket.size();
#pragma omp parallel for
                for( int idx_pair = 0; idx_pair < n_pairs_bucket; idx_pair++ )
                {
                    pairwise_hg_ij[idx_bucket][idx_pair] = 0;
                }
                idx_bucket++;

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

                    // if the constraint is violated change adjusted_positions
                    // this is achieved by solving
                    //   | (adjusted_positions[i] - h g r_ij / mi) - (adjusted_positions[j] + h g r_ij / mj) |^2 = dij^2
                    // for g.
                    // When neglecting g^2 terms in the derivation, we get (also h cancels so we don't need the step size)
                    const double hg = ( s2 - dij2 ) / ( 2.0 * ( s.dot( rij ) ) * ( 1.0 / mi + 1.0 / mj ) );

                    adjusted_positions.row( i ) += -hg * rij / mi;
                    adjusted_positions.row( j ) -= -hg * rij / mj;

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
            for( const auto & bucket : buckets )
            {
                auto cb = [&]( int idx_pair ) {
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