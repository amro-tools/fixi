#pragma once
#include <cstddef>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Fixi::Backend
{

inline void set_num_threads( int n_threads )
{
#ifdef _OPENMP
    omp_set_num_threads( n_threads );
#endif
}

inline int get_num_threads()
{
#ifdef _OPENMP
    return omp_get_max_threads();
#endif
    return 1;
}

// We define the reduction ourselves here (instead of declaring a custom omp reduction), because
// we don't want to have to declare a new reduction for each template instantiation and this
// approach works seamlessly with templates
template<typename result_t, typename CallbackT, typename ReductionT>
result_t transform_reduce( size_t n_max, const CallbackT & cb, const ReductionT & red, result_t zero = result_t{ 0 } )
{
    result_t result = zero;

#pragma omp parallel
    {
        // Each thread has its own private sum
        result_t private_result = zero;

// Distribute the loop iterations among threads
#pragma omp for nowait
        for( size_t m = 0; m < n_max; m++ )
        {
            red( private_result, cb( m ) );
        }

// Combine the private sums into the global sum
#pragma omp critical
        {
            red( result, private_result );
        }
    }

    return result;
}
} // namespace Fixi::Backend