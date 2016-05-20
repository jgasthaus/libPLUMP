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
#include "libplump/stirling.h"
#include <cmath>
#include <sstream>
#include <gsl/gsl_sf_gamma.h>

namespace gatsby { namespace libplump {


/**
 * Compute S_d(c,t) recursively by directly applying the recursion
 * S_d(c,t) = S_d(c-1,t-1) + (c-1 - t*d)S_d(c-1,t)
 *
 * This does a _lot_ of duplicate work and should only be used for debugging!
 */
double gen_stirling_recursive(double d, int c, int t) {
    //tracer << "gen_stirling_recursive(" << d << ", " << c << ", " << t << ")" << std::endl;
    // base cases
    if ((c==1 && t==1) || (c==0 && t==0))
        return 1.;
    if (c==0 || t==0)
        return 0.;
    if (t>c)
        return 0.;
    // recursion
    return gen_stirling_recursive(d,c-1,t-1) + (c-1 - t*d)*gen_stirling_recursive(d,c-1,t);
}

/**
 * Compute S_d(c,t) recursively in log space by directly applying the recursion
 * S_d(c,t) = S_d(c-1,t-1) + (c-1 - t*d)S_d(c-1,t)
 *
 * This does a _lot_ of duplicate work and should only be used for debugging!
 */
double log_gen_stirling_recursive(double d, int c, int t) {
    //tracer << "log_gen_stirling_recursive(" << d << ", " << c << ", " << t << ")" << std::endl;
    // base cases
    if ((c==1 && t==1) || (c==0 && t==0))
        return 0;
    if (c==0 || t==0)
        return -INFINITY;
    if (t>c)
        return -INFINITY;
    // recursion
    return fast_logsumexp(log_gen_stirling_recursive(d,c-1,t-1),log((c-1 - t*d)) + log_gen_stirling_recursive(d,c-1,t));
}

/**
 * Compute S_d(c,t) directly in log space using the recursion
 * S_d(c,t) = S_d(c-1,t-1) + (c-1 - t*d)S_d(c-1,t)
 */
double log_gen_stirling_direct(double d, int c, int t) {
    //tracer << "log_gen_stirling_direct" << d << ", " << c << ", " << t << ")" << std::endl;
    // base cases
    if ((c==1 && t==1) || (c==0 && t==0))
        return 0;
    if (c==0 || t==0)
        return -INFINITY;
    if (t>c)
        return -INFINITY;

    // column to the left of the one that we are currently computing 
    // in the c x t matrix
    d_vec left_col(c,-INFINITY);
    d_vec cur_col(c,0);
    left_col[0] = 0;


    for (int cur_t = 1;cur_t <= t; ++cur_t) {
        for(int j = 1; j < (c - t) + cur_t; ++j) {
            // c = j + cur_t
            cur_col[j] = fast_logsumexp(left_col[j], log(j+cur_t-1 - cur_t*d) + cur_col[j-1]);
        }

        for (int i = 1; i < (c - t) + cur_t; ++i) {
            left_col[i] = cur_col[i];
        }
    }
    return cur_col[c-t];
}

double log_gen_stirling_ratio(double d, int c, int t) {
    //XXX:Make sure these are correct!
    // base cases
    if ((c==1 && t==1) || (c==0 && t==0))
        return 0;
    if (c==0 || t==0)
        return -INFINITY;
    if (t>c)
        return -INFINITY;

    typedef std::vector<double> d_vec;
    // column to the left of the one that we are currently computing 
    // in the c x t matrix
    d_vec left_col(c,-INFINITY);
    d_vec cur_col(c,0);
    left_col[0] = 0;


    for (int cur_t = 1; cur_t <= t; ++cur_t) {
        for(int j = 1; j < (c - t) + cur_t; ++j) {
            // c = j + cur_t
            cur_col[j] = fast_logsumexp(left_col[j], log(j+cur_t-1 - cur_t*d) + cur_col[j-1]);
        }
        if (cur_t != t) {
            for (int i = 1; i < (c - t) + cur_t; ++i) {
                left_col[i] = cur_col[i];
            }
        }
    }
    return left_col[c-t] - cur_col[c-t];
}

d_vec_vec log_gen_stirling_table(double d, int c) {
    d_vec_vec out;
    out.reserve(c - 1);
    // column to the left of the one that we are currently computing 
    // in the c x t matrix
    d_vec left_col(c, -INFINITY);
    d_vec cur_col(c, 0);
    left_col[0] = 0;
    const int t = c;
 
    for (int cur_t = 1; cur_t < t; ++cur_t) {
        for(int j = 1; j < (c - cur_t) + 1; ++j) {
            // c = j + cur_t
            cur_col[j] = fast_logsumexp(left_col[j],
                                        log(j + cur_t - 1 - cur_t*d) 
                                        + cur_col[j-1]);
        }
        if (cur_t != t) {
            for (int i = 1; i < c - cur_t + 1; ++i) {
                left_col[i] = cur_col[i];
            }
        }
        d_vec tmp(c - cur_t, 0);
        for (int i = 0; i < c - cur_t; ++i) {
            tmp[i] = cur_col[i + 1];
        }
        out.push_back(tmp);
    }
    return out;
}

void log_gen_stirling_table_extend(double d, int c, d_vec_vec& table) {
    // column to the left of the one that we are currently computing 
    // in the c x t matrix
    int prev_c = table.size() + 1;
    const int t = c;

    // treat the first column specially
    for (int j = prev_c - 1; j < c - 1; ++j) {
        table[0].push_back(log(j + 1 - d) + table[0][j - 1]);
    }

    for (int cur_t = 2; cur_t < t; ++cur_t) {
        if (cur_t>prev_c-1) {
            table.push_back(d_vec());
        }
        for(int j = std::max<int>((int)prev_c - (int)cur_t, 0);
            j < (c - cur_t); ++j) {
            // c = j + cur_t + 1
            if (j == 0) {
                table[cur_t-1].push_back(fast_logsumexp(table[cur_t - 2][j], 
                                                        log(cur_t - cur_t*d)));
            } else {
                table[cur_t-1].push_back(
                    fast_logsumexp(table[cur_t - 2][j],
                                   log(j+cur_t - cur_t*d)
                                   + table[cur_t - 1][j - 1]));
            }
        }
    }
}


////////////////////////// stirling_generator_recompute_log ///////////////////
stirling_generator_recompute_log::stirling_generator_recompute_log(double d, 
                                                                   int c,
                                                                   int t) 
    : d(d) {}

double stirling_generator_recompute_log::ratio(int c, int t) {
    return log_gen_stirling_ratio(d,c,t);
}


std::string stirling_generator_recompute_log::statsToString() {
    return std::string();
}

////////////////////////// stirling_generator_fast_log ////////////////////////
void stirling_generator_fast_log::incC() {
    d_vec new_col(current_c-1,0);
    new_col[0] = (fast_logsumexp(row.back(),log(current_c - 2*d) + col.front()));
    row.push_back(log(current_c - d) + row.back());
    for(int i = 1; i < current_c - 2; ++i) {
       new_col[i] = fast_logsumexp(col[i-1], log(current_c - (i+2)*d) + col[i]); 
    }
    new_col[current_c-2] = (fast_logsumexp(col.back(), log(current_c - current_c*d)));
    prev_col = col;
    col = new_col;
    current_c++;
}

void stirling_generator_fast_log::decC() {
    col = prev_col;
    row.pop_back();
    prev_col.pop_back();
    if (prev_col.size() == 0) {
        current_c--;
        return; // decrease from 4 to 3
    }
    prev_col[0] = fast_logminusexp(col[0],row[current_c-3]) - log(current_c-2 - 2*d); 
    for(int i = 1; i < current_c - 4; ++i) {
        prev_col[i] = fast_logminusexp(col[i],prev_col[i-1]) - log(current_c-2 - (i+2)*d);  
    }
    current_c--;
}

stirling_generator_fast_log::stirling_generator_fast_log(double d, int c, int t) : d(d) {
    row.push_back(0); // 1 customer, 1 table
    row.push_back(log(1-d)); // 2 customers, one table
    row.push_back(log((2-d)*(1-d))); // 3 customers, one table
    col.push_back(fast_logsumexp(log(1-d), log(2-2*d))); // 3 customers, 2 tables
    current_c = 3;
}

double stirling_generator_fast_log::get(int c, int t) {
    if (t==1) {
        return row[c-1];
    }
    if (c == current_c) {
        return col[t-2];
    }
    if (c == current_c - 1) {
        return prev_col[t-2];
    }
    return 0; // should never happen
}

double stirling_generator_fast_log::ratio(int c, int t) {
    if (c<2 || t<2)
        return NAN;
    if (t>c)
        return NAN;
    if(c==t)
        return 1;
    
    while(current_c < c) {
        //std::cout << "incrementing c";
        incC();
    }
    while(current_c > c) {
        decC();
    }

    //std::cout << "current_c: " << current_c << std::endl;
    //std::cout << iterableToString(prev_col) << std::endl <<iterableToString(col) << std::endl << iterableToString(row) << std::endl;
   
    return exp(get(c-1,t-1) - get(c,t));

}

std::string stirling_generator_fast_log::statsToString() {
    return std::string();
}


////////////////////////// stirling_generator_full_log ////////////////////////
stirling_generator_full_log::stirling_generator_full_log(double d, int c, int t) : d(d) {
    // generate the smallest possible table for a start
    // TODO: do something more clever
   table = log_gen_stirling_table(d,2);
   c_max = 2;
   num_construct++;
}

double stirling_generator_full_log::ratio(int c, int t) {
    num_ratio_calls++;
    if (t==1) {
        return 0;
    }
    if (c<2 || t<2)
        return NAN;
    if (t>c)
        return NAN;
    if(c==t)
        return 1;
    if (c>c_max) {
        num_extends++;
        log_gen_stirling_table_extend(d,c,table);
        c_max = c;
        if (c_max > global_c_max) {
            global_c_max = c_max;
        }
    }
    if (c==0 || t==0) {
        std::cout << "Error: calling ratio(0,0)" << std::endl;
    }
    //for (int i=0; i<table.size();i++) {
    //    std::cout << iterableToString(table[i]) << std::endl;
    //}
    return exp(log_get_stirling_from_table(table, c-1, t-1) - log_get_stirling_from_table(table, c, t));
}


double stirling_generator_full_log::getLog(int c, int t) {
    if (c>c_max) {
        num_extends++;
        log_gen_stirling_table_extend(d,c,table);
        c_max = c;
        if (c_max > global_c_max) {
            global_c_max = c_max;
        }
    }
    return log_get_stirling_from_table(table, c, t);
}

std::string stirling_generator_full_log::statsToString() {
    std::ostringstream out;
    out << "Constructed: " << num_construct << ", extended: " << num_extends << ", calls: " << num_ratio_calls << ", c_max " << global_c_max;
    return out.str();
}


int stirling_generator_full_log::global_c_max = 0;
int stirling_generator_full_log::num_ratio_calls = 0;
int stirling_generator_full_log::num_extends = 0;
int stirling_generator_full_log::num_construct = 0;

}} // namespace gatsby::libplump
