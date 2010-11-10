/*
 * Copyright 2009, 2010 Jan Gasthaus (j.gasthaus@gatsby.ucl.ac.uk)
 * 
 * This file is part of libPLUMP.
 * 
 * libPLUMP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * libPLUMP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with libPLUMP.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef RANDOM_H_
#define RANDOM_H_

#include <gsl/gsl_rng.h>
#include <vector>
#include <cassert>
#include <algorithm> // for lower_bound

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
    double z = gsl_rng_uniform(global_rng)*pdf[end_pos];

    assert((z >= 0) && (z <= pdf[end_pos]));
    
    // Perform binary search for z using std::lower_bound.
    // lower_bound(begin, end, x) returns the first element within [begin,end)
    // that is equal or larger than x.
    return std::lower_bound(pdf.begin(), pdf.begin() + end_pos + 1, z)
         - pdf.begin();

}

} } // namespace gatsby::libplump

#endif // RANDOM_H_
