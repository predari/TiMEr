/*
 * permutation.h
 *
 *  Created on: Dec 21, 2016
 *      Author: Roland
 */

#ifndef PERMUTATION_H_
#define PERMUTATION_H_

#include <stdlib.h>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm> 
#include <iostream>

using namespace std;

class Permutation {
    private:
    vector<int64_t> p;
    int64_t length;

    public:

    Permutation(int64_t length);

    ~Permutation();
    
    /**
     * returns p[pos]
     * @param[in] position pos.
     * @return p[pos].
     */
    inline int64_t entry(int64_t pos) {
        return(p[pos]);
    }
    
    /**
     * Initialize permutation to identity reversed,
     * Least significant bit mapped to 0 and
     * most significant bit mapped to length - 1.
     */
    void ini(void);
    
    /**
     * Sets permutation to the inverse of pi.
     * @param[in] permutation pi.
     * @return w with its digits permuted as determined by rho.
     */
    void invert(Permutation pi);
    
    /**
     * Creates random permutation through shuffling.
     * @param[in] vector p.
     * @param[out] shuffled p.
     */
    void shuffle(void);
    
    /**
     * permutes digits of int64_t w as determined by permutation.
     * @param[in] int64_t w.
     * @return w with its digits permuted.
     */
    int64_t permute(int64_t v);

    /**
     * Displays a permutation on standard out.
     * @param[in] permutation pi.
     */
    void display(void);

};

#endif /* PERMUTATION_H_ */
