#if defined(SEQUENTIAL)
#include "timer.h"
#include "sequential.h"
#elif defined(USE_CILK)
#include "timer.h"
#include "usecilk.h"
#else
#include "usehclib.h"
#endif

#ifndef NUMA_WS
#include "my_malloc.h"
#endif
