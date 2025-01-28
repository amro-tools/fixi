#pragma once
#include <Eigen/Core>

namespace Fixi
{

using Vector3 = Eigen::Matrix<double, 1, 3>;

// Row-major for compatibility with numpy arrays. See https://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html
using Vectorfield = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::StorageOptions::RowMajor>;

struct FixedBondLengthPair
{
    int i{};
    int j{};
    double dij{};
    double mass_i{ 1.0 };
    double mass_j{ 1.0 };
};
} // namespace Fixi