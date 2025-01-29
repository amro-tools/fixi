#include "catch2/matchers/catch_matchers.hpp"
#include "fixi/utils.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <cstdio>
#include <fixi/defines.hpp>
#include <fixi/rattle.hpp>
#include <format>
#include <iostream>

TEST_CASE( "Test that the positions and velocities of two atoms can be constrained using RATTLE", "[SingleBond]" )
{
    int n_atoms = { 2 };
    double rij  = 1.0;

    const int maxiter      = 20;
    const double tolerance = 1e-6;
    const double m         = 1.0;
    const std::array<bool, 3> pbc{ true, true, false };
    const Fixi::Vector3 cell_lengths = { 20.0, 20.0, 20.0 };

    Fixi::Vectorfield positions( n_atoms, 3 );
    Fixi::Vectorfield velocities( n_atoms, 3 );
    std::vector<Fixi::FixedBondLengthPair> pairs{};

    Fixi::Scalarfield masses = Fixi::Scalarfield( n_atoms );
    masses( 0 )              = m;
    masses( 1 )              = m;

    positions.row( 0 ) = Fixi::Vector3::Zero();
    positions.row( 1 ) = positions.row( 0 ) + 2.0 * rij * Fixi::Vector3::Random().normalized();

    velocities.row( 0 ) = Fixi::Vector3::Random();
    velocities.row( 1 ) = Fixi::Vector3::Random();

    pairs.push_back( { 0, 1, rij } );

    auto rattle = Fixi::Rattle( maxiter, tolerance, pairs );

    auto adjusted_positions = positions;

    auto [hij_before, hij_v_before]
        = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );

    INFO( std::format( "hij_before {}\n", hij_before ) );

    rattle.adjust_positions( positions, adjusted_positions, masses, cell_lengths, pbc );
    INFO( std::format( "iterations {}\n", rattle.iteration ) );

    auto [hij, hij_v] = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );
    INFO( std::format( "hij_after {}\n", hij ) );

    REQUIRE_THAT( hij, Catch::Matchers::WithinAbs( 0.0, tolerance ) );

    auto [hij2_before, hij_v2_before]
        = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );
    INFO( std::format( "hij_v_before {}\n", hij_v2_before ) );

    rattle.adjust_velocities( adjusted_positions, velocities, masses, cell_lengths, pbc );
    INFO( std::format( "iterations {}\n", rattle.iteration ) );

    auto [hij2_after, hij_v2_after]
        = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );
    INFO( std::format( "hij_v_after {}\n", hij_v2_after ) );

    REQUIRE_THAT( hij_v2_after, Catch::Matchers::WithinAbs( 0.0, tolerance ) );
}