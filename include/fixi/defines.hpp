#pragma once
#include <Eigen/Core>

namespace Fixi
{

using Vector3 = Eigen::Matrix<double, 1, 3>;
using Matrix3 = Eigen::Matrix<double, 3, 3>;

// Row-major for compatibility with numpy arrays. See https://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html
using Vectorfield = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::StorageOptions::RowMajor>;
using Scalarfield = Eigen::Matrix<double, Eigen::Dynamic, 1>;

struct FixedBondLengthPair
{
    int i{};
    int j{};
    double dij{};
};
} // namespace Fixi