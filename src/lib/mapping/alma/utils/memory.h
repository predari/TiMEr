#ifndef MEMORY_H_
#define MEMORY_H_

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cassert>

#define MemoryAllocate(T, n) (T *)malloc((n) * sizeof(T))
#define MemoryRelocate(p, T, n) (T *)realloc((void *)p, ((n) * sizeof(T)))
#define MemoryFree(p, T, n) free((void *)p)
#define MemoryClear(p, T, n) memset(p, 0, sizeof(T) * n);

#define MemoryCopy(d, s, T, n) memcpy((void *)d, (void *)s, (n) * sizeof(T))
#define MemoryFill(p, T, n) memset((void *)p, -1, (n) * sizeof(T))


#define xmalloc(p, T, n) { assert((p = MemoryAllocate(T, n)) != NULL); }
#define xfree(p, T, n) {if(p) {free((void *) p); (p)=NULL;}}
#define xcalloc(p, T, n) { assert(p = MemoryClearAllocate(T, n)) != NULL; }

#define xnew(p, T) { p = new T; assert(p != NULL); /* printf("new %p\n", p); */}
#define xnewarr(p, T, n) { p = new T[n]; assert(p != NULL); /* printf("new[] %p\n", p); */ }
#define xdel(p) { if (p != NULL) { /*printf("del %p\n", p);*/ delete p; (p)=NULL;} }
#define xdelarr(p) { if (p != NULL) { /*printf("del[] %p\n", p);*/ delete[] p; (p)=NULL;} }
#define xcleararr(p, n) { memset(p, 0, sizeof(*p) * n); }

#endif /* MEMORY_H_ */
