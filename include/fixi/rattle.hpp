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
    const int maxiter{};
    const double tolerance{};
    const std::vector<FixedBondLengthPair> pairs{};
    // const std::array<bool, 3> pbc{};
    const std::vector<std::vector<FixedBondLengthPair>> buckets{};

public:
    int iteration{ 0 };

    Rattle( int maxiter, double tolerance, const std::vector<FixedBondLengthPair> & pairs )
            : maxiter( maxiter ), tolerance( tolerance ), pairs( pairs ), buckets( bucket_sorter( pairs ) )
    {
    }

    void adjust_positions(
        const Eigen::Ref<Vectorfield> unadjusted_positions, Eigen::Ref<Vectorfield> adjusted_positions,
        Eigen::Ref<Scalarfield> masses, const Vector3 & cell_lengths, const std::array<bool, 3> & pbc )
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

                    // if the constraint is violated adjust adjusted_positions
                    // this is achieved by iteratively solving
                    //   |adjusted_positions[i] - adjusted_positions[j]|^2 = dij^2
                    // Neglect g^2 term in derivation (also h cancels so we don't need the step size)
                    const double g = ( s2 - dij2 ) / ( 2.0 * ( s.dot( rij ) ) * ( 1.0 / mi + 1.0 / mj ) );

                    adjusted_positions.row( i ) += -g * rij / mi;
                    adjusted_positions.row( j ) -= -g * rij / mi;

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
    }

    void adjust_velocities(
        const Eigen::Ref<Vectorfield> positions, Eigen::Ref<Vectorfield> adjusted_velocities,
        Eigen::Ref<Scalarfield> masses, const Vector3 & cell_lengths, const std::array<bool, 3> & pbc )
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
    }
};

} // namespace Fixi