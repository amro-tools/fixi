#ifdef _OPENMP
#include <omp.h>
#endif

#include <chrono>
#include <fixi/defines.hpp>
#include <fixi/rattle.hpp>
#include <fixi/utils.hpp>
#include <format>
#include <iostream>

int get_n_threads()
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1.0;
#endif
}

int main()
{
    int n_molecules = { 5000 };
    int n_atoms     = { 3 * n_molecules };

    double r_OH = 1;
    // double r_HH = 1.5;

    double mO = 16.0;
    double mH = 1.0;

    const int maxiter      = 500;
    const double tolerance = 1e-10;
    const std::array<bool, 3> pbc{ false, false, false };
    const Fixi::Vector3 cell_lengths = { 20.0, 20.0, 20.0 };

    Fixi::Vectorfield positions( n_atoms, 3 );
    Fixi::Vectorfield velocities( n_atoms, 3 );
    Fixi::Scalarfield masses( n_atoms );

    std::vector<Fixi::FixedBondLengthPair> pairs{};

    // Define positions
    for( int i = 0; i < n_molecules; i++ )
    {
        const int idx_O  = 3 * i;
        const int idx_H1 = idx_O + 1;
        const int idx_H2 = idx_O + 2;

        masses( idx_O )  = mO;
        masses( idx_H1 ) = mH;
        masses( idx_H2 ) = mH;

        positions.row( idx_O )  = 2.0 * Fixi::Vector3::Random();
        positions.row( idx_H1 ) = positions.row( idx_O ) + r_OH * 1.4 * Fixi::Vector3::Random().normalized();
        positions.row( idx_H2 ) = positions.row( idx_O ) + r_OH * 1.15 * Fixi::Vector3::Random().normalized();
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
        // pairs.push_back( { idx_H1, idx_H2, r_HH, mH, mH } );
    }

    std::cout << std::format( "{} threads\n", get_n_threads() );

    auto rattle = Fixi::Rattle( maxiter, tolerance, pairs );

    Fixi::Vectorfield adjusted_positions = positions;

    auto [hij_before, hij_v_before]
        = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );

    std::cout << "hij_before " << hij_before << "\n";

    auto t1 = std::chrono::high_resolution_clock::now();
    rattle.adjust_positions( positions, adjusted_positions, masses, cell_lengths, pbc );
    auto t2 = std::chrono::high_resolution_clock::now();
    /* Getting number of milliseconds as an integer. */
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 );
    std::cout << std::format( " duration = {}\n", ms_int );
    std::cout << std::format( " iterations = {}\n", rattle.iteration );

    auto [hij, hij_v] = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );
    std::cout << "hij_after " << hij << "\n";

    auto [hij2_before, hij_v2_before]
        = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );
    std::cout << "hij_v_before " << hij_v2_before << "\n";
    rattle.adjust_velocities( adjusted_positions, velocities, masses, cell_lengths, pbc );
    auto [hij2_after, hij_v2_after]
        = Fixi::check_constraints( pairs, adjusted_positions, velocities, cell_lengths, pbc );
    std::cout << "hij_v_after " << hij_v2_after << "\n";
}