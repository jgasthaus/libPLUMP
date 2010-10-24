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

#include "libplump/pyp_sample.h"

#include <cmath>
#include <cassert>
#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>

#include "libplump/random.h"
#include "libplump/utils.h"

namespace gatsby { namespace libplump {

////////////////////////////////////////////////////////////////////////////////
//        PYP SAMPLING FUNCTIONS FOR GENERATING SEATING ARRANGEMENTS          //
////////////////////////////////////////////////////////////////////////////////

/**
 * Sample a seating arrangement with c customers around t tables.
 *
 * Returns the number of customers at each table in the resulting sample.
 * 
 * Runtime: O(c x t)
 */
std::vector<int> sample_crp_ct(double d, int c, int t) {
  tracer << "sample_crp_ct(" << d << ", " << c << ", " << t << ")" << std::endl;
  assert(c!=0);
  assert(c>=t);

  if (t==1) {
    std::vector<int> out(1,c); // all c customers sit at the first table
    return out;
  }

  typedef boost::multi_array<double, 2> array_type;
  //typedef array_type::index index;
  typedef int index;
  array_type grid(boost::extents[t][c]);
  for(index i = 0; i < t; ++i) 
    for(index j = 0; j <  c; ++j)
      grid[i][j] = 0;

  // forward filtering
  grid[0][0] = 1;
  grid[t-1][c-1] = 1;
  for (index i = 1; i < c - 1; ++i) {
    for (index j = 0; j < t; ++j) {
      if (j>0) {
        grid[j][i] += grid[j-1][i-1];
      }
      grid[j][i] += grid[j][i-1] * (i - (j+1)*d);
    }
  }

  // backward sampling
  std::vector<int> z(c,0);
  int cur = t-1;
  z[c-1] = t - 1;
  d_vec probs(t,0);
  for(index i = c - 2; i >= 0; --i) {
    probs[cur] = (i+1-(cur+1)*d) * grid[cur][i];
    probs[cur-1] = grid[cur-1][i];
    cur = sample_unnormalized_pdf(probs);
    z[i] = cur;
    if (cur == 0) {
      break;
    }
    probs.assign(t,0);
  }

  // sample customers given z
  std::vector<int> arrangement(t,0);
  arrangement[0] = 1;
  for (int i=1; i<c; ++i) {
    probs.assign(t,0);
    if (z[i] == z[i-1] + 1) {
      probs[z[i]] = 1;
    } else {
      for(int j=0;j<z[i]+1;++j) {
        probs[j] = (arrangement[j] - d)/(i - (z[i]+1)*d);
      }
    }
    ++arrangement[sample_unnormalized_pdf(probs,z[i])];
  }

  assert(sum(arrangement) == c && (int)arrangement.size() == t);
  for(int i = 0; i < (int)arrangement.size(); ++i) {
    assert(arrangement[i] != 0);
  }
  return arrangement;
}


std::vector<int> sample_crp_c(double d, double a, int c) {
  d_vec probs(c,0);
  std::vector<int> arrangement;
  arrangement.push_back(1); // first customer at first table
  for(int i = 1; i < c; ++i) {
    for(int j = 0; j < (int)arrangement.size(); ++j) {
      probs[j] = arrangement[j] - d;
    }
    probs[arrangement.size()] = a + arrangement.size()*d;
    int sample = sample_unnormalized_pdf(probs,arrangement.size());
    if (sample == (int)arrangement.size()) {
      arrangement.push_back(1); // new table
    } else {
      ++arrangement[sample];
    }
  }
  return arrangement;
}

}} // namespace gatsby::libplump
