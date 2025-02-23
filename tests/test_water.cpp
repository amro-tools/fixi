#include "catch2/matchers/catch_matchers.hpp"
#include "fixi/utils.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <fixi/defines.hpp>
#include <fixi/rattle.hpp>

TEST_CASE( "Test that the positions and velocities can be constrained using RATTLE", "[TestWater]" )
{
    int n_atoms     = { 27 };
    int n_molecules = { 9 };

    double r_OH = 1;
    double r_HH = 1.5;

    double mO = 16.0;
    double mH = 1.0;

    const int maxiter      = 5000;
    const double tolerance = 1e-6;
    const std::array<bool, 3> pbc{ false, false, false };
    const Fixi::Vector3 cell_lengths = { 20.0, 20.0, 20.0 };

    Fixi::Vectorfield positions( n_atoms, 3 );
    Fixi::Vectorfield velocities( n_atoms, 3 );
    Fixi::Scalarfield masses( n_atoms );
    std::vector<Fixi::FixedBondLengthPair> pairs{};

    // Define positions and masses
    for( int i = 0; i < n_molecules; i++ )
    {
        const int idx_O  = 3 * i;
        const int idx_H1 = idx_O + 1;
        const int idx_H2 = idx_O + 2;

        positions.row( idx_O )  = 2.0 * Fixi::Vector3::Random();
        positions.row( idx_H1 ) = positions.row( idx_O ) + 0.4 * Fixi::Vector3::Random();
        positions.row( idx_H2 ) = positions.row( idx_O ) + 0.4 * Fixi::Vector3::Random();

        masses( idx_O )  = mO;
        masses( idx_H1 ) = mH;
        masses( idx_H2 ) = mH;
    }

    // Define velocities
    for( int i = 0; i < n_molecules; i++ )
    {
        const int idx_O  = 3 * i;
        const int idx_H1 = idx_O + 1;
        const int idx_H2 = idx_O + 2;

        velocities.row( idx_O )  = 0.1 * Fixi::Vector3::Random();
        velocities.row( idx_H1 ) = 0.1 * Fixi::Vector3::Random();
        velocities.row( idx_H2 ) = 0.1 * Fixi::Vector3::Random();
    }

    // Define pairs
    for( int i = 0; i < n_molecules; i++ )
    {
        const int idx_O  = 3 * i;
        const int idx_H1 = idx_O + 1;
        const int idx_H2 = idx_O + 2;

        pairs.push_back( { idx_O, idx_H1, r_OH } );
        pairs.push_back( { idx_O, idx_H2, r_OH } );
        pairs.push_back( { idx_H1, idx_H2, r_HH } );
    }

    auto rattle = Fixi::Rattle( maxiter, tolerance, pairs );

    Fixi::Vectorfield adjusted_positions = positions;

    auto [hij_before, hij_v_before]
        = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );

    INFO( std::format( "hij_before {}\n", hij_before ) );

    rattle.adjust_positions( positions, adjusted_positions, masses, cell_lengths, pbc );
    INFO( std::format( "iterations {}\n", rattle.get_iteration() ) );

    auto [hij, hij_v] = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );
    INFO( std::format( "hij_after {}\n", hij ) );

    REQUIRE_THAT( hij, Catch::Matchers::WithinAbs( 0.0, tolerance ) );

    auto [hij2_before, hij_v2_before]
        = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );
    INFO( std::format( "hij_v_before {}\n", hij_v2_before ) );

    rattle.adjust_velocities( adjusted_positions, velocities, masses, cell_lengths, pbc );
    INFO( std::format( "iterations {}\n", rattle.get_iteration() ) );

    auto [hij2_after, hij_v2_after]
        = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );
    INFO( std::format( "hij_v_after {}\n", hij_v2_after ) );

    REQUIRE_THAT( hij_v2_after, Catch::Matchers::WithinAbs( 0.0, tolerance ) );
}