#pragma once
#include <fixi/defines.hpp>

inline Fixi::Vector3
mic( const Fixi::Vector3 & dr, const Fixi::Vector3 & cell_lengths, const std::array<bool, 3> & pbc )
{
    Fixi::Vector3 wrapped_dr = dr;
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