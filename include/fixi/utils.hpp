#pragma once
#include <array>
#include <fixi/defines.hpp>
#include <vector>

namespace Fixi
{
inline Vector3 mic( const Vector3 & dr, const Vector3 & cell_lengths, const std::array<bool, 3> & pbc )
{
    Vector3 wrapped_dr = dr;
    // Use minimum image convention
    for( int k = 0; k < 3; ++k )
    {
        if( pbc[k] )
        {
            wrapped_dr( k ) -= cell_lengths( k ) * std::round( dr( k ) / cell_lengths( k ) );
        }
    } // end of getting the relative distance

    return wrapped_dr;
}

inline std::array<double, 2> check_constraints(
    const std::vector<FixedBondLengthPair> & pairs, const Eigen::Ref<Vectorfield> positions,
    const Eigen::Ref<Vectorfield> velocities, const Vector3 & cell_lengths, const std::array<bool, 3> & pbc )
{
    double max_hij   = {};
    double max_hij_v = {};

    const int n_pairs = pairs.size();
    for( int idx_pair = 0; idx_pair < n_pairs; idx_pair++ )
    {
        const int i      = pairs[idx_pair].i;
        const int j      = pairs[idx_pair].j;
        const double dij = pairs[idx_pair].dij; // desired constrained bond length

        Vector3 s         = mic( positions.row( i ) - positions.row( j ), cell_lengths, pbc );
        double s2         = s.squaredNorm();
        const double dij2 = dij * dij;
        const double hij  = abs( s2 - dij2 ); // holonomic constraint
        max_hij           = std::max( max_hij, hij );

        const Vector3 vij  = velocities.row( i ) - velocities.row( j );
        const Vector3 rij  = positions.row( i ) - positions.row( j );
        const double hij_v = abs( rij.dot( vij ) );
        max_hij_v          = std::max( max_hij_v, hij_v );
    }

    return { max_hij, max_hij_v };
}
} // namespace Fixi