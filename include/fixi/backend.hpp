#include <cstddef>
namespace Fixi::Backend
{
// We define the reduction ourselves here (instead of declaring a custom omp reduction), because
// we don't want to have to declare a new reduction for each template instantiation and this
// approach works seamlessly with templates
template<typename result_t, typename CallbackT, typename ReductionT>
result_t transform_reduce( size_t n_max, const CallbackT & cb, const ReductionT & red )
{
    result_t result{ 0.0 };

#pragma omp parallel
    {
        // Each thread has its own private sum
        result_t private_result{ 0.0 };

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