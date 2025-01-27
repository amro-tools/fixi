#pragma once
#include <Eigen/Core>

namespace Fixi
{

using Vector3     = Eigen::Matrix<double, 3, 1>;
using Vectorfield = Eigen::Matrix<double, Eigen::Dynamic, 3>;

struct FixedBondLengthPair
{
    int i{};
    int j{};
    double dij{};
    double mass_i{ 1.0 };
    double mass_j{ 1.0 };
};
} // namespace Fixi