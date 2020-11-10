namespace hclib {

template<typename T>
void async(T &&lambda) {
    lambda();
}

template<typename T>
void async_hinted(void* mem, unsigned int start, unsigned int end, T &&lambda) {
    lambda();
}

#define HCLIB_FINISH
#define FORASYNC_MODE_RECURSIVE 1

typedef struct {
    int low;
    int high;
    int stride;
    int tile;
} loop_domain_t;

template<typename T1>
void __divide2_1D(int low, int high, int threshold, T1 && lambda);

template<typename T1>
void __divide4_1D(int low, int high, int threshold, T1 && lambda){
    if(high-low > threshold && high-low>=4) {
        int chunk = (high - low) / 4;
        __divide2_1D(low, low+chunk, threshold, lambda);
        __divide2_1D(low+chunk, low + 2 * chunk, threshold, lambda);
        __divide4_1D(low + 2 * chunk, low + 3 * chunk, threshold, lambda);
        __divide2_1D(low + 3 * chunk, high, threshold, lambda);
    } else {
        for(int i=low; i<high; i++) {
            lambda(i);
        }
    }
}

template<typename T1>
void __divide2_1D(int low, int high, int threshold, T1 && lambda){
    if(high-low > threshold) {
        __divide4_1D(low, (low+high)/2, threshold, lambda);              
        __divide2_1D((low+high)/2, high, threshold, lambda);             
    } else {
        for(int i=low; i<high; i++) {
            lambda(i);
        }
    }
}

template<typename T1, typename T2, typename T3>
void irregular_recursion1D(loop_domain_t* loop, void* array, T1 && start, T2 && end, T3 && lambda, int mode) {
    assert(loop->stride == 1);
    __divide2_1D(loop->low, loop->high, loop->tile, lambda);
}

template<typename T1>
void __divide_1D(int low, int high, int threshold, T1 && lambda){
    if(high-low > threshold) {
        __divide_1D(low, (low+high)/2, threshold, lambda);              
        __divide_1D((low+high)/2, high, threshold, lambda);             
    } else {
        for(int i=low; i<high; i++) {
            lambda(i);
        }
    }
}

template<typename T1, typename T2, typename T3>
void recursion1D(loop_domain_t* loop, void* array, T1 && start, T2 && end, T3 && lambda, int mode) {
    assert(loop->stride == 1);
    __divide_1D(loop->low, loop->high, loop->tile, lambda);
}

int num_workers() { return 1; }

int current_worker() { return 0; }

template<typename T>
void finish(T &&lambda) {
    lambda();
}

template<typename T>
void firsttouch_random_initialization(int low, int high, T* array1d, bool decimal) {
    unsigned int seed = 1;
    for(int j=low; j<high; j++) {
        int num = rand_r(&seed);
        array1d[j] = decimal? T(num/RAND_MAX) : T(num);
    }
}

}
