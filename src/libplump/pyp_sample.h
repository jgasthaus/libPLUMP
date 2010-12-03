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
