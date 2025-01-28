#include <catch2/catch_test_macros.hpp>
#include <fixi/buckets.hpp>
#include <fixi/defines.hpp>
#include <format>
#include <iostream>

bool check_bucket( const std::vector<Fixi::FixedBondLengthPair> & bucket )
{
    for( int idx1 = 0; idx1 < bucket.size(); idx1++ )
    {
        const auto & p1 = bucket[idx1];
        for( int idx2 = 0; idx2 < bucket.size(); idx2++ )
        {
            // skip if pair is checking against itself
            if( idx1 == idx2 )
            {
                continue;
            }

            // make sure indices are completely disjoint
            const auto & p2 = bucket[idx2];
            if( p1.i == p2.i || p1.j == p2.j || p1.i == p2.j || p1.j == p2.i )
            {
                return false;
            }
            return true;
        }
    }
}

TEST_CASE( "Test pairs are correctly sorted into buckets", "[PairsInBuckets]" )
{
    std::vector<Fixi::FixedBondLengthPair> pairs
        = { { 1, 2 }, { 1, 3 }, { 1, 4 },   { 4, 5 },   { 4, 6 },      { 7, 8 },
            { 7, 9 }, { 8, 9 }, { 10, 11 }, { 10, 12 }, { 10, 11, 12 } };

    const int n_pairs  = pairs.size();
    const auto buckets = bucket_sorter( pairs );

    int n_buckets = buckets.size();

    INFO( std ::format( "There are {} buckets\n", n_buckets ) );

    for( int ibucket = 0; ibucket < n_buckets; ibucket++ )
    {
        std::cout << std::format( "ibucket = {}\n", ibucket ) << "\n";
        const auto & bucket = buckets[ibucket];

        for( int ipair = 0; ipair < bucket.size(); ipair++ )
        {
            const auto pair = bucket[ipair];
            std::cout << std::format( "    i = {}, j = {}\n", pair.i, pair.j ) << "\n";
        }
    }

    REQUIRE( n_buckets == 3 );

    int n_pairs_buckets = 0;
    for( const auto & bucket : buckets )
    {
        n_pairs_buckets += bucket.size();
        REQUIRE( check_bucket( bucket ) );
    }

    REQUIRE( n_pairs == n_pairs_buckets );
}