namespace hclib {

template<typename T>
T* _numa_malloc(size_t count) {
    return (T*) malloc(sizeof(T) * count);
}
template<typename T>
T** _numa_malloc(size_t row, size_t col) {
    T** var = (T**) malloc(sizeof(T*) * row);
    for(int j=0; j<row; j++) {
        var[j] = (T*) malloc(sizeof(T) * col);
    }
    return var;
}
template<typename T>
void _numa_free(T* mem, size_t row) {
    for(int j=0; j<row; j++) {
        free(mem[j]);
    }
    free(mem);
}
template<typename T>
void _numa_free(T* mem) {
        free(mem);
}
#define numa_malloc _numa_malloc
#define numa_interleave_malloc _numa_malloc
#define numa_free _numa_free
}
