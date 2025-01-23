#pragma once
#include "fixi/defines.hpp"
#include <fixi/utils.hpp>
#include <iostream>
#include <span>

namespace Fixi
{

struct FixedBondLengthPair
{
    int i{};
    int j{};
    double dij{};
    double mass_i{ 1.0 };
    double mass_j{ 1.0 };
};

class Rattle
{
    int maxiter{};
    double tolerance{};
    std::vector<FixedBondLengthPair> pairs{};
    std::array<bool, 3> pbc{};

public:
    Rattle(
        int maxiter, double tolerance, const std::vector<FixedBondLengthPair> & pairs, const std::array<bool, 3> & pbc )
            : maxiter( maxiter ), tolerance( tolerance ), pairs( pairs ), pbc( pbc ){};

    void adjust_positions(
        const std::span<Vector3> unadjusted_positions, std::span<Vector3> adjusted_positions,
        const Vector3 & cell_lengths )
    {
        const int n_pairs = pairs.size();

        for( int idx_pair = 0; idx_pair < n_pairs; idx_pair++ )
        {
            const int i      = pairs[idx_pair].i;
            const int j      = pairs[idx_pair].j;
            const double mi  = pairs[idx_pair].mass_i;
            const double mj  = pairs[idx_pair].mass_j;
            const double dij = pairs[idx_pair].dij; // desired constrained bond length

            // `s` is the current approximation for the vector displacement between atoms i and j
            const Vector3 rij = mic( unadjusted_positions[i] - unadjusted_positions[j], cell_lengths, pbc );
            Vector3 s         = mic( adjusted_positions[i] - adjusted_positions[j], cell_lengths, pbc );
            double s2         = s.squaredNorm();
            const double dij2 = dij * dij;

            double hij = abs( s2 - dij2 ); // holonomic constraint

            // if the constraint is violated adjust adjusted_positions
            // this is achieved by iteratively solving
            //   |adjusted_positions[i] - adjusted_positions[j]|^2 = dij^2
            // Neglect g^2 term in derivation (also h cancels so we don't need the step size)
            int iteration{ 0 };
            while( hij > tolerance && iteration < maxiter )
            {
                const double g = ( s2 - dij2 ) / ( 2.0 * ( s.dot( rij ) ) * ( 1.0 / mi + 1.0 / mj ) );

                adjusted_positions[i] += g * rij / mi;
                adjusted_positions[j] -= g * rij / mi;

                // Update s and s2
                s   = mic( adjusted_positions[i] - adjusted_positions[j], cell_lengths, pbc );
                s2  = s.squaredNorm();
                hij = abs( s2 - dij2 );
                iteration++;
            }
        }
    }

    void adjust_momenta(
        const std::span<Vector3> positions, std::span<Vector3> adjusted_velocities, const Vector3 & cell_lengths )
    {
        const int n_pairs = pairs.size();

        for( int idx_pair = 0; idx_pair < n_pairs; idx_pair++ )
        {
            const int i       = pairs[idx_pair].i;
            const int j       = pairs[idx_pair].j;
            const double mi   = pairs[idx_pair].mass_i;
            const double mj   = pairs[idx_pair].mass_j;
            const double dij  = pairs[idx_pair].dij; // desired constrained bond length
            const double dij2 = dij * dij;

            const Vector3 rij = mic( positions[i] - positions[j], cell_lengths, pbc );
            Vector3 vij       = adjusted_velocities[i] - adjusted_velocities[j];

            double hij = abs( rij.dot( vij ) );

            int iteration{ 0 };
            while( hij > tolerance && iteration < maxiter )
            {
                const double k = rij.dot( vij ) / ( dij2 * ( 1.0 / mi + 1.0 / mj ) );

                adjusted_velocities[i] -= k * rij / mi;
                adjusted_velocities[i] += k * rij / mj;

                // Update
                vij = adjusted_velocities[i] - adjusted_velocities[j];
                hij = abs( rij.dot( vij ) );
                iteration++;
            }
        }
    }
};

} // namespace Fixi