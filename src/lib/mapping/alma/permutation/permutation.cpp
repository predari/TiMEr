/*
 * permutation.cpp
 *
 *  Created on: Jan 02, 2017
 *      Author: Roland
*/

#include "permutation.h"
#include "../../alma_graph_access.h"
#include "../utils/arith.h"

#include <iostream>     // std::cout
#include <iterator>     // std::ostream_iterator

using namespace std;

/*****************************************/
/*****************************************/
void Permutation::ini(void) {
    for(int64_t pos = 0; pos < length; pos++) {
        p[pos] = pos;
    }
}

/*****************************************/
/*****************************************/
void Permutation::invert(Permutation pi) {
    for(int64_t pos = 0; pos < length; pos++) {
        p[pi.entry(pos)] = pos;
    }
}

/*****************************************/
/*****************************************/
void Permutation::shuffle(void) {
    for(int64_t i = length - 1; i >= 1; i--) {
        int64_t j = randomFunction(i+1);
        int64_t p_i = p[i];
        p[i] = p[j];
        p[j] = p_i;
    }
}
    
/*****************************************/
/*****************************************/
int64_t Permutation::permute(int64_t l) {
    int64_t lengthDec = length - 1;
    int64_t ret = 0;
    for(int64_t i = length - 1; i >= 0; i--) {
        ret+= (l % 2) << (lengthDec- p[i]);
        l = l >> 1;
    }
    return(ret);
}
/*****************************************/
/*****************************************/
void Permutation::display(void) {
    copy(p.begin(), p.end(), ostream_iterator<int>(cout, " "));
    cout << endl;
}

/*****************************************/
/*****************************************/
Permutation::Permutation(int64_t l) {
    length = l;
    p.resize(l);
}

/*****************************************/
/*****************************************/
Permutation::~Permutation(void) {
}
