#include <fixi/defines.hpp>
#include <fixi/rattle.hpp>
#include <fixi/utils.hpp>

#define PYBIND11_DETAILED_ERROR_MESSAGES

// Bindings
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>
#include <format>

// Namespaces
using namespace std::string_literals; // For ""s
using namespace pybind11::literals;   // For ""_a
namespace py = pybind11;              // Convention

PYBIND11_MODULE( fixicpp, m )
{
    py::class_<Fixi::FixedBondLengthPair>( m, "FixedBondLengthPair" )
        .def( py::init<>() )
        .def_readwrite( "i", &Fixi::FixedBondLengthPair::i )
        .def_readwrite( "j", &Fixi::FixedBondLengthPair::j )
        .def_readwrite( "dij", &Fixi::FixedBondLengthPair::dij )
        .def(
            "__str__",
            []( const Fixi::FixedBondLengthPair & p ) { return std::format( "i={}, j={}, dij={}", p.i, p.j, p.dij ); } )
        .def( py::pickle(
            []( const Fixi::FixedBondLengthPair & p ) { // __getstate__
                /* Return a tuple that fully encodes the state of the object */

                return py::make_tuple( p.i, p.j, p.dij );
            },
            []( py::tuple t ) { // __setstate__
                if( t.size() != 3 )
                    throw std::runtime_error( "Invalid state for FixedBondLengthPair object!" );

                const int i      = t[0].cast<int>();
                const int j      = t[1].cast<int>();
                const double dij = t[2].cast<double>();

                auto p = Fixi::FixedBondLengthPair();
                p.i    = i;
                p.j    = j;
                p.dij  = dij;
                return p;

            } ) );

    py::class_<Fixi::Rattle>( m, "Rattle" )
        .def( py::init<int, double, std::vector<Fixi::FixedBondLengthPair>>() )
        .def_readwrite( "maxiter", &Fixi::Rattle::maxiter )
        .def_readwrite( "tolerance", &Fixi::Rattle::tolerance )
        .def( "get_iteration", &Fixi::Rattle::get_iteration )
        .def( "get_pairs", &Fixi::Rattle::get_pairs )
        .def( "adjust_positions", &Fixi::Rattle::adjust_positions )
        .def( "adjust_velocities", &Fixi::Rattle::adjust_velocities )
        .def( py::pickle(
            []( const Fixi::Rattle & p ) { // __getstate__
                /* Return a tuple that fully encodes the state of the object */
                return py::make_tuple( p.maxiter, p.tolerance, p.get_pairs() );
            },
            []( py::tuple t ) { // __setstate__
                if( t.size() != 3 )
                    throw std::runtime_error( "Invalid state for Rattle object!" );

                int maxiter      = t[0].cast<int>();
                double tolerance = t[1].cast<double>();
                auto pairs       = t[2].cast<std::vector<Fixi::FixedBondLengthPair>>();

                return Fixi::Rattle( maxiter, tolerance, pairs );

            } ) );

    m.def( "check_constraints", &Fixi::check_constraints, "Checks constraints." );
    m.def( "set_num_threads", &Fixi::Backend::set_num_threads );
    m.def( "get_num_threads", &Fixi::Backend::get_num_threads );
}