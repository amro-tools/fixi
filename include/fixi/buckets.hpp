#pragma once
#include <fixi/defines.hpp>
#include <fixi/utils.hpp>
#include <set>

namespace Fixi
{

// Sort into buckets so that buckets can be handled in parallel
inline std::vector<std::vector<FixedBondLengthPair>> bucket_sorter( const std::vector<FixedBondLengthPair> & pairs )
{
    std::vector<std::vector<FixedBondLengthPair>> buckets{};
    std::vector<std::set<int>> bucket_sets{};

    auto exists_in_bucket = [&]( const auto & pair, const auto & bucket_set ) {
        // for i
        if( bucket_set.find( pair.i ) != bucket_set.end() )
        {
            return true;
        }
        // for j
        if( bucket_set.find( pair.j ) != bucket_set.end() )
        {
            return true;
        }

        return false;
    };

    auto add_pair_to_bucket = [&]( const auto & pair, auto & bucket_set, auto & bucket ) {
        bucket.push_back( pair );
        // Update the set with i and j in pair
        bucket_set.insert( pair.i );
        bucket_set.insert( pair.j );
    };

    for( const auto & pair : pairs )
    {
        bool pair_added = false;

        for( size_t b_idx = 0; b_idx < buckets.size(); b_idx++ )
        {
            const bool pair_found = exists_in_bucket( pair, bucket_sets[b_idx] );
            // If the pair has not been found in the current bucket, add it
            if( !pair_found )
            {
                add_pair_to_bucket( pair, bucket_sets[b_idx], buckets[b_idx] );
                pair_added = true;
                break;
            }
        }

        // If the pair has not been added to any bucket, create a new one and add the pair
        if( !pair_added )
        {
            buckets.push_back( std::vector<FixedBondLengthPair>{ pair } );
            bucket_sets.push_back( std::set<int>{ pair.i, pair.j } );
        }
    }

    return buckets;
}

} // namespace Fixi