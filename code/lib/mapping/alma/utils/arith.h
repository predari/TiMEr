#ifndef ARITH_H_
#define ARITH_H_

#include<stdint.h>
#include<assert.h>

using namespace std;

int64_t floorLog2(int64_t x);
int64_t ceilLog2(int64_t x);
int64_t flipPosition(int64_t pos, int64_t l);
int64_t hammingDistance(int64_t v, int64_t w);
int64_t randomFunction(int64_t l);

#endif /* ARITH_H_ */
