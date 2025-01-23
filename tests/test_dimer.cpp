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
    double m    = 1.0;

    const int maxiter      = 20;
    const double tolerance = 1e-6;
    const std::array<bool, 3> pbc{ true, true, false };
    const Fixi::Vector3 cell_lengths = { 20.0, 20.0, 20.0 };

    std::vector<Fixi::Vector3> positions( n_atoms );
    std::vector<Fixi::Vector3> velocities( n_atoms );
    std::vector<Fixi::FixedBondLengthPair> pairs{};

    positions[0] = Fixi::Vector3::Zero();
    positions[1] = positions[0] + 2.0 * rij * Fixi::Vector3::Random().normalized();

    velocities[0] = Fixi::Vector3::Random();
    velocities[1] = Fixi::Vector3::Random();

    pairs.push_back( { 0, 1, rij, m, m } );

    auto rattle = Fixi::Rattle( maxiter, tolerance, pairs, pbc );

    auto adjusted_positions = positions;

    auto [hij_before, hij_v_before]
        = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );

    INFO( std::format( "hij_before {}\n", hij_before ) );

    rattle.adjust_positions( positions, adjusted_positions, cell_lengths );

    auto [hij, hij_v] = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );
    INFO( std::format( "hij_after {}\n", hij ) );

    REQUIRE_THAT( hij, Catch::Matchers::WithinAbs( 0.0, tolerance ) );

    auto [hij2_before, hij_v2_before]
        = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );
    INFO( std::format( "hij_v_before {}\n", hij_v2_before ) );

    rattle.adjust_velocities( adjusted_positions, velocities, cell_lengths );
    auto [hij2_after, hij_v2_after]
        = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );
    INFO( std::format( "hij_v_after {}\n", hij_v2_after ) );

    REQUIRE_THAT( hij_v2_after, Catch::Matchers::WithinAbs( 0.0, tolerance ) );
}