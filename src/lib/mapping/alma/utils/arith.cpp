#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include <bitset>

#include "./arith.h"

int64_t floorLog2(int64_t x) {
    assert(x >= 1);
    int64_t y = x;
    int64_t exp = 0;
    while(y > 1) {
        y = y >> 1;
        exp++;
    }
    return(exp);
}

int64_t ceilLog2(int64_t x) {
    int64_t f = floorLog2(x);
    if(x == (1LL << f)) {
        return(f);
    } else {
        return(f + 1);
    }
}

int64_t flipPosition(int64_t pos, int64_t l) {
  int64_t ll = l;
  int64_t entryAtPosition = (ll >> pos) % 2;
  if(entryAtPosition == 1) {
    ll-= (1LL << pos);
  } else {
    ll+= (1LL << pos);
  }
  return ll;
}

int64_t hammingDistance(int64_t l1, int64_t l2) {
  std::bitset<64> b1(l1);
  std::bitset<64> b2(l2);
  return((b1 ^ b2).count());
}

int64_t randomFunction(int64_t l) {
    mt19937 generator((unsigned int)time(0));
    return (generator() % l);
}


