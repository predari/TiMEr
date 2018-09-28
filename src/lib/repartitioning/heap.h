//include dibapheader//
#ifndef HEAP_H
#define HEAP_H

#include <cstring>
#include <cstdlib>

#define MemoryAllocate(T, n) (T *)malloc((n) * sizeof(T))
#define MemoryRelocate(p, T, n) (T *)realloc((void *)p, ((n) * sizeof(T)))
#define MemoryFree(p, T, n) free((void *)p)
#define MemoryCopy(d, s, T, n) memcpy((void *)d, (void *)s, (n) * sizeof(T))
#define MemoryClear(p, T, n) memset((void *)p, 0, (n) * sizeof(T))
#define MemoryFill(p, T, n) memset((void *)p, -1, (n) * sizeof(T))
#define MemoryClearAllocate(T, n) (T *)calloc((n), sizeof(T))
#define SafeMemoryAllocate(p, T, n) { assert((p = MemoryAllocate(T, n)) != NULL); }
#define SafeMemoryFree(p, T, n) {if(p) {free((void *) p); (p)=NULL;}}
#define SafeMemoryClearAllocate(p, T, n) { assert((p = MemoryClearAllocate(T, n)) != NULL); }
#define SafeMemoryDelete(p) { if(p) { delete(p); (p)=NULL; }}

struct Heap
{
  int size;
  int * heap;
  int * location;
  float * priority;
  int maxsize;
};

namespace heap {

struct Heap * init(const int n);

void dispose(struct Heap *const h);

void insert(struct Heap *const h, const int x, const float p);

void remove(struct Heap *const h, const int x);

int size(const struct Heap *const h);

int top(struct Heap *const h);

float priority(struct Heap *const h, const int x);

}

#endif
