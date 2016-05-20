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

#ifndef PYP_SAMPLE_H_
#define PYP_SAMPLE_H_

#include <vector>

// only for debugging
#include <iostream>
#include "libplump/utils.h"

namespace gatsby { namespace libplump {

////////////////////////////////////////////////////////////////////////////////
//        PYP SAMPLING FUNCTIONS FOR GENERATING SEATING ARRANGEMENTS          //
////////////////////////////////////////////////////////////////////////////////

std::vector<int> sample_crp_z_fb(double d, int c, int t);

std::vector<int> sample_crp_z_bf(double d, int c, int t);

std::vector<int> sample_crp_given_z(double d, std::vector<int>& z);

inline std::vector<int> sample_crp_ct_fb(double d, int c, int t) {
  std::vector<int> z = sample_crp_z_fb(d, c, t);
  return sample_crp_given_z(d, z);
}

inline std::vector<int> sample_crp_ct_bf(double d, int c, int t) {
  std::vector<int> z = sample_crp_z_bf(d, c, t);
  return sample_crp_given_z(d, z);
}

/**
 * Sample a seating arrangement with c customers around t tables.
 * Uses backward filtering, forward sampling (stable, preferred).
 *
 * Returns the number of customers at each table in the resulting sample.
 * 
 * Runtime: O(c x t)
 */
inline std::vector<int> sample_crp_ct(double d, int c, int t) {
  return sample_crp_ct_bf(d, c, t);
}

/**
 * Sample a seating arrangement of c customers from a CRP with discount d and
 * concentration parameter a.
 *
 * @returns A vector that contains an element for each table in the arrangement
 *          whose value is the number of customers at that table.
 */
std::vector<int> sample_crp_c(double d, double a, int c);

}} // namespace gatsby::libplump

#endif
