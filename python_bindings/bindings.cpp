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
        .def_readwrite( "mass_i", &Fixi::FixedBondLengthPair::mass_i )
        .def_readwrite( "mass_j", &Fixi::FixedBondLengthPair::mass_j )
        .def_readwrite( "dij", &Fixi::FixedBondLengthPair::dij );

    py::class_<Fixi::Rattle>( m, "Rattle" )
        .def( py::init<int, double, std::vector<Fixi::FixedBondLengthPair>, std::array<bool, 3>>() )
        .def( "adjust_positions", &Fixi::Rattle::adjust_positions )
        .def( "adjust_velocities", &Fixi::Rattle::adjust_velocities );

    m.def( "check_constraints", &Fixi::check_constraints, "Checks constraints." );
}