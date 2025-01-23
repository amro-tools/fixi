#include "fixi/utils.hpp"
#include <cstdio>
#include <fixi/defines.hpp>
#include <fixi/rattle.hpp>
#include <iostream>

int main()
{
    int n_atoms     = { 9 };
    int n_molecules = { 3 };

    double r_OH = 1;
    double r_HH = 1.5;

    double mO = 16.0;
    double mH = 1.0;

    const int maxiter      = 20;
    const double tolerance = 1e-6;
    const std::array<bool, 3> pbc{ true, true, true };
    const Fixi::Vector3 cell_lengths = { 20.0, 20.0, 20.0 };

    std::vector<Fixi::Vector3> positions( n_atoms );
    std::vector<Fixi::Vector3> velocities( n_atoms );
    std::vector<Fixi::FixedBondLengthPair> pairs{};

    // Define positions
    for( int i = 0; i < n_molecules; i++ )
    {
        const int idx_O  = 3 * i;
        const int idx_H1 = idx_O + 1;
        const int idx_H2 = idx_O + 2;

        positions[idx_O]  = 2.0 * Fixi::Vector3::Random();
        positions[idx_H1] = positions[idx_O] + 0.4 * Fixi::Vector3::Random();
        positions[idx_H2] = positions[idx_O] + 0.4 * Fixi::Vector3::Random();
    }

    // Define velocities
    for( int i = 0; i < n_molecules; i++ )
    {
        const int idx_O  = 3 * i;
        const int idx_H1 = idx_O + 1;
        const int idx_H2 = idx_O + 2;

        velocities[idx_O]  = 0.1 * Fixi::Vector3::Random();
        velocities[idx_H1] = 0.1 * Fixi::Vector3::Random();
        velocities[idx_H2] = 0.1 * Fixi::Vector3::Random();
    }

    // Define pairs
    for( int i = 0; i < n_molecules; i++ )
    {
        const int idx_O  = 3 * i;
        const int idx_H1 = idx_O + 1;
        const int idx_H2 = idx_O + 2;

        pairs.push_back( { idx_O, idx_H1, r_OH, mO, mH } );
        pairs.push_back( { idx_O, idx_H2, r_OH, mO, mH } );
        pairs.push_back( { idx_H1, idx_H2, r_HH, mH, mH } );
    }

    auto rattle = Fixi::Rattle( maxiter, tolerance, pairs, pbc );

    auto adjusted_positions = positions;

    auto [hij_before, hij_v_before]
        = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );

    std::cout << "hij_before " << hij_before << "\n";

    rattle.adjust_positions( positions, adjusted_positions, cell_lengths );

    auto [hij, hij_v] = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );
    std::cout << "hij_after " << hij << "\n";
}