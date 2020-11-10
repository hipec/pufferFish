#include "hclib.hpp"
#include <iostream>

namespace hclib {

#ifndef NUMA_WS
template <typename T>
inline void _async_hinted(void* memory, size_t start_index, size_t end_index, T &&lambda) {
    MARK_OVH(current_ws()->id);
    typedef typename std::remove_reference<T>::type U;
    hclib_async(hclib::lambda_wrapper<U>, new U(lambda), nullptr, nullptr, nullptr, 0);
}
#define async_hinted _async_hinted
#endif

#ifdef NUMA_WS_ELASTIC
template <typename T>
inline void _async_elastic(void* memory, size_t start_index, size_t end_index, T &&lambda) {
    MARK_OVH(current_ws()->id);
    place_t* pl = place_with_max_pages(memory, start_index, end_index);
    if(serial_elision_enabled(pl)) {
        lambda();
    } else {
        typedef typename std::remove_reference<T>::type U;
        hclib_async(lambda_wrapper<U>, new U(lambda), nullptr, nullptr, pl, 0);
    }
}
#define async_hinted _async_elastic
#endif

template<typename T1, typename T2, typename T3>
void __divide2_1D(int low, int high, int threshold, void* array, T1&& start, T2&& end, T3 && lambda);

template<typename T1, typename T2, typename T3>
void __divide4_1D(int low, int high, int threshold, void* array, T1&& start, T2&& end, T3 && lambda){
    if(high-low > threshold && high-low>=4) {
	int chunk = (high - low) / 4;
        HCLIB_FINISH {
            async_hinted(array, start(low), end(low+chunk), [=]() {
                __divide2_1D(low, low+chunk, threshold, array, start, end, lambda);
	    });
            async_hinted(array, start(low+chunk), end(low+2*chunk), [=]() {
                __divide2_1D(low+chunk, low + 2 * chunk, threshold, array, start, end, lambda);
	    });
            async_hinted(array, start(low+2*chunk), end(low+3*chunk), [=]() {
                __divide4_1D(low + 2 * chunk, low + 3 * chunk, threshold, array, start, end, lambda);
	    });
            async_hinted(array, start(low+3*chunk), end(high), [=]() {
                __divide2_1D(low + 3 * chunk, high, threshold, array, start, end, lambda);
	    });
	}
    } else {
        for(int i=low; i<high; i++) {
            lambda(i);
	}
    }
}

template<typename T1, typename T2, typename T3>
void __divide2_1D(int low, int high, int threshold, void* array, T1&& start, T2&& end, T3 && lambda){
    if(high-low > threshold) {
        HCLIB_FINISH {
            async_hinted(array, start(low), end((high+low)/2), [=]() {
	        __divide4_1D(low, (low+high)/2, threshold, array, start, end, lambda);		    
	    });
            async_hinted(array, start((low+high)/2), end(high), [=]() {
	        __divide2_1D((low+high)/2, high, threshold, array, start, end, lambda);		    
	    });
	}
    } else {
        for(int i=low; i<high; i++) {
            lambda(i);
	}
    }
}

template<typename T1, typename T2, typename T3>
void irregular_recursion1D(loop_domain_t* loop, void* array, T1 && start, T2 && end, T3 && lambda, int mode) {
    assert(loop->stride == 1);
    __divide2_1D(loop->low, loop->high, loop->tile, array, start, end, lambda);
}

template<typename T1, typename T2, typename T3>
void __divide_1D(int low, int high, int threshold, void* array, T1&& start, T2&& end, T3 && lambda){
    if(high-low > threshold) {
        HCLIB_FINISH {
            async_hinted(array, start(low), end((high+low)/2), [=]() {
	        __divide_1D(low, (low+high)/2, threshold, array, start, end, lambda);		    
	    });
            async_hinted(array, start((low+high)/2), end(high), [=]() {
	        __divide_1D((low+high)/2, high, threshold, array, start, end, lambda);		    
	    });
	}
    } else {
        for(int i=low; i<high; i++) {
            lambda(i);
	}
    }
}

template<typename T1, typename T2, typename T3>
void recursion1D(loop_domain_t* loop, void* array, T1 && start, T2 && end, T3 && lambda, int mode) {
    assert(loop->stride == 1);
    __divide_1D(loop->low, loop->high, loop->tile, array, start, end, lambda);
}

template<typename T>
void firsttouch_random_initialization(int low, int high, T* array1d, bool decimal) {
    int numWorkers= hclib::num_workers();
    assert(high%numWorkers == 0);
    int batchSize = high / numWorkers;
    int chunk=0;
    HCLIB_FINISH {
    for(int wid=0; wid<numWorkers; wid++) { 
        int start = wid * batchSize;
        int end = start + batchSize;
        hclib::async_hinted(array1d, start, end-1, [=]() {
	    unsigned int seed = wid+1;
            for(int j=start; j<end; j++) {
	        int num = rand_r(&seed);
	        array1d[j] = decimal? T(num/RAND_MAX) : T(num);
            }
        });
    }
    }
}


}
