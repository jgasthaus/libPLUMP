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

/** 
 * Functions for computing and tabulating the generalized stirling numbers 
 * of type (-1,-d,0). 
 *
 */

#ifndef STIRLING_H_
#define STIRLING_H_
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <gsl/gsl_sf_gamma.h>


namespace gatsby { namespace libplump {

typedef std::vector<double> d_vec;
typedef std::vector<d_vec> d_vec_vec;

/**
 * computes log(exp(a) + exp(b)) while avoiding numerical
 * instabilities.
 */
inline double fast_logsumexp(double a, double b) {
  if (a>b) {
    return log(1+exp(b-a)) + a;
  } else {
    return log(1+exp(a-b)) + b;
  }
}

/**
 * computes log(exp(a) - exp(b)) while avoiding numerical
 * instabilities.
 */
inline double fast_logminusexp(double a, double b) {
  if (a>b) {
    return log(1-exp(b-a)) + a;
  } else {
    return log(1-exp(a-b)) + b;
  }
}

/**
 * Compute S_d(c,t) recursively by directly applying the recursion
 * S_d(c,t) = S_d(c-1,t-1) + (c-1 - t*d)S_d(c-1,t)
 *
 * This does a _lot_ of duplicate work and should only be used for debugging!
 */
double gen_stirling_recursive(double d, int c, int t);


/**
 * Compute S_d(c,t) recursively in log space by directly applying the recursion
 * S_d(c,t) = S_d(c-1,t-1) + (c-1 - t*d)S_d(c-1,t)
 *
 * This does a _lot_ of duplicate work and should only be used for debugging!
 */
double log_gen_stirling_recursive(double d, int c, int t);


/**
 * Compute S_d(c,t) directly in log space using the recursion
 * S_d(c,t) = S_d(c-1,t-1) + (c-1 - t*d)S_d(c-1,t)
 */
double log_gen_stirling_direct(double d, int c, int t);


double log_gen_stirling_ratio(double d, int c, int t);


d_vec_vec log_gen_stirling_table(double d, int c);


void log_gen_stirling_table_extend(double d, int c, d_vec_vec& table);


double log_get_stirling_from_table(d_vec_vec& table, int c, int t);


class stirling_generator_recompute_log {
  public:
    stirling_generator_recompute_log(double d, int c, int t);

    double ratio(int c, int t);

    static std::string statsToString();

  private:
    double d;
};


class stirling_generator_fast_log {

  public:
    stirling_generator_fast_log(double d, int c, int t);

    double get(int c, int t);

    double ratio(int c, int t);

    static std::string statsToString();

  private:
    void incC();
    void decC();

    double d;
    int current_c;
    d_vec row, col, prev_col;
};

/**
 * Class to encapsulate the generation of ratios of stirling numbers
 * for removing customers from restaurants.
 * 
 * One of these objects should be constructed for each restaurant, 
 * providing the current number of tables, and the total number of customers.
 */
class stirling_generator_full_log {
  public:
    stirling_generator_full_log(double d, int c, int t);

    double ratio(int c, int t);
    
    double getLog(int c, int t);

    static std::string statsToString();

  private:

    static int global_c_max, num_ratio_calls, num_extends, num_construct;

    d_vec_vec table;
    int c_max;
    double d;
};

inline double log_get_stirling_from_table(d_vec_vec& table, int c, int t) {
  // c and t must be non-negative
  assert(c >= 0 && t>= 0);

  if ((c == 1  && t == 1) || (c == 0 && t == 0))
    return 0;
  if (c == 0 || t == 0)
    return -INFINITY;
  if (t > c)
    return -INFINITY;
  if(c == t)
    return 0;
  if (t > (int)table.size() || (c - t) > (int)table[t-1].size()) {
    std::cout<< c << ", " << t << std::endl;
  }
  return table[t-1][c-t-1];
}


inline double log_stirling_asymptotic(double d, int c, int t) {
  return   gsl_sf_lngamma(c) 
         - gsl_sf_lngamma(1 - d) 
         - gsl_sf_lngamma(t)
         - (t-1) * log(d)
         - d * log(c);
}

// class log_gen_stirling_table {
//     private:
//         std::vector<std::vector<double> > table;
//         size_t c_cur, t_cur;
//        
//         /*
//          * Grow the arrays to the required size.
//          */
//         void _grow(size_t c_max, size_t t_max) {
//             for(int i=t_cur;i<t_max;i++) {
//                 std::vector<double> tmp(c_cur,0);
//                 table.push_back(std::vector<double>());
//             }
//         }
// 
//         /*
//          * Tabulate up to the specified limits
//          */
//         void _tabulate(size_t c_max, size_t t_max) {
//             for t in xrange(1,max_t):
//         for c in xrange(t,max_c):
//             a[c,t] = a[c-1,t-1] + (c-1-d*t)*a[c-1,t]
// 
//         }
// 
//     public:
//         double d;
//         gen_stirling_table(double d) : d(d), c_cur(0), t_cur(0)  {}
// 
// };






}} // namespace gatsby::libplump


#endif /* STIRLING_H_ */
