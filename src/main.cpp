#include <fixi/rattle.hpp>

int main()
{

    std::vector<Fixi::FixedBondLengthPair> pairs = { { 0, 1, 0.1 }, { 2, 0, { 0.3 } } };

    auto rattle = Fixi::Rattle( 100, 1e-3, pairs, { true, true, true } );

    // rattle.ad
}