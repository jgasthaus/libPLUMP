/*
 * Copyright 2008-2016 Jan Gasthaus
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef RANDOM_H_
#define RANDOM_H_

#include <gsl/gsl_rng.h>
#include <vector>
#include <cassert>
#include <algorithm> // for lower_bound
// only for debugging
#include <iostream>
#include "libplump/utils.h"

namespace gatsby { namespace libplump {

/*******************************************************************************
 * RANDOM NUMBER GENERATION AND SAMPLING BASED ON GSL
 ******************************************************************************/

/**
 * Initialize global RNG.
 *
 * This function has to be called before using any of the sampling functions.
 *
 * The RNG type and seed are determined by the environment variables
 * GSL_RNG_TYPE and GSL_RNG_SEED.
 */
void init_rng();

/**
 * Free the global RNG.
 */
void free_rng();

/**
 * Returns true with probability true_prob.
 */
bool coin(double true_prob);

/**
 * Returns a uniform integer between 0 and max-1.
 */
long int uniform_int(long int max);

/**
 * Sample from a discrete distribution on 0,...,MAX with the given PDF.
 *
 * The probability if returning x for x \in 0,...,MAX is given by 
 *   pdf[x] / (\sum_{i=0,...,end_pos} pdf[i])
 * i.e. pdf is normalized so that the sum of all elements up to and including
 * element end_pos is 1.
 *  
 * Algorithm:
 *   1) Compute CDF; normalizing constant Z = CDF[end_pos]
 *   2) Sample z ~ Uniform(0,Z)
 *   3) find the smallest element x of CDF that is larger than i using binary
 *      search
 *   4) return x
 *
 * Complexity: O(log MAX)
 */
int sample_unnormalized_pdf(std::vector<double> pdf, int end_pos = 0);




////////////////////////////////////////////////////////////////////////////////
//////////////////////    INLINE FUNCTION DEFINITIONS   ////////////////////////
////////////////////////////////////////////////////////////////////////////////

extern gsl_rng* global_rng;

/**
 * Returns true with probability true_prob.
 */
inline bool coin(double true_prob) {
    return (true_prob>gsl_rng_uniform(global_rng));
}

/**
 * Returns a uniform integer between 0 and max-1.
 */
inline long int uniform_int(long int max) {
    return gsl_rng_uniform_int(global_rng, max);
}

/**
 * Sample from a discrete distribution on 0,...,MAX with the given PDF.
 *
 * The probability if returning x for x \in 0,...,MAX is given by 
 *   pdf[x] / (\sum_{i=0,...,end_pos} pdf[i])
 * i.e. pdf is normalized so that the sum of all elements up to and including
 * element end_pos is 1.
 *  
 * Algorithm:
 *   1) Compute CDF; normalizing constant Z = CDF[end_pos]
 *   2) Sample z ~ Uniform(0,Z)
 *   3) find the smallest element x of CDF that is larger than i using binary
 *      search
 *   4) return x
 *
 * Complexity: O(log MAX)
 */
inline int sample_unnormalized_pdf(std::vector<double> pdf, int end_pos) {
    assert(pdf.size() > 0);
    assert(end_pos < pdf.size());
    assert(end_pos >= 0);

    // if end_pos == 0, use entire vector
    if (end_pos == 0) {
        end_pos = pdf.size()-1;
    }
    
    // compute CDF (inplace)
    for (int i = 0; i < end_pos; ++i) {
        assert(pdf[i] >= 0);
        pdf[i+1] += pdf[i];
    }

    assert(pdf[end_pos] > 0);

    // sample pos ~ Unigorm(0,Z)
    double z = gsl_rng_uniform_pos(global_rng)*pdf[end_pos];

    assert((z >= 0) && (z <= pdf[end_pos]));
    
    // Perform binary search for z using std::lower_bound.
    // lower_bound(begin, end, x) returns the first element within [begin,end)
    // that is equal or larger than x.
    int x = std::lower_bound(pdf.begin(), pdf.begin() + end_pos + 1, z)
         - pdf.begin();

    assert(x == 0 || pdf[x-1] != pdf[x]);
    assert(x < pdf.size());
    return x;

}

} } // namespace gatsby::libplump

#endif // RANDOM_H_
