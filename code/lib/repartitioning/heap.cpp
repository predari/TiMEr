//include dibapheader//
#include "heap.h"
#include <cfloat>

namespace heap {

struct Heap * init(const int n)
{
  struct Heap * h = MemoryAllocate(struct Heap, 1);
  h->size = 0;
  h->heap = MemoryAllocate(int, n + 1);
  h->location = MemoryAllocate(int, n);
  h->priority = MemoryAllocate(float, n + 1);  
  h->priority++;
  h->heap[0] = -1;
  h->priority[-1] = FLT_MAX;
  MemoryFill(h->location, int, n);
  h->maxsize = n;
  return h;
}

void dispose(struct Heap *const h)
{
  h->priority--;
  MemoryFree(h->priority, float, h->maxsize + 1);
  MemoryFree(h->location, int, h->maxsize);
  MemoryFree(h->heap, int, h->maxsize + 1);
  MemoryFree(h, struct Heap, 1);
}

void insert(struct Heap *const h, const int x, const float p)
{
  h->size++;
  
  int i = h->size;
  int j = i >> 1;
  
  while (h->priority[h->heap[j]] < p)
    {
      h->heap[i] = h->heap[j];
      h->location[h->heap[i]] = i;
      i = j;
      j = i >> 1;
    }
  h->heap[i] = x;
  h->location[x] = i;
  h->priority[x] = p;
}

void remove(struct Heap *const h, const int x)
{
  int i = h->location[x];

  if (i < 0)
    return;

  const float p = h->priority[h->heap[h->size]];
  
  while (i <= (h->size >> 1))
    {
      int j = i + i;
      
      if (j < h->size && h->priority[h->heap[j]] < h->priority[h->heap[j + 1]])
	j++;
      if (p >= h->priority[h->heap[j]])
	break;
      h->heap[i] = h->heap[j];
      h->location[h->heap[i]] = i;
      i = j;
    }
  h->heap[i] = h->heap[h->size];
  h->location[h->heap[i]] = i;
  
  h->location[x] = -1;
  h->size--;
}

int size(const struct Heap *const h)
{
  return h->size;
}

int top(struct Heap *const h)
{
  return h->heap[1];
}

float priority(struct Heap *const h, const int x)
{
  return h->priority[x];
}

}
