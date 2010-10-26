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

/**
 * The classes in this file represent and store the parameters for the
 * HPYP models (discount and concentration parameters).
 **/

#ifndef HPYP_PARAMETERS_H_
#define HPYP_PARAMETERS_H_

#include "libplump/config.h"
#include "libplump/context_tree.h" // for WrappedNodeList
#include "libplump/utils.h" // for WrappedNodeList

namespace gatsby { namespace libplump {

typedef std::vector<double> d_vec;

class IParameters {
  public:

    virtual ~IParameters() {}
    /**
     * Get discount parameters for each node in the node list.
     */
    virtual d_vec getDiscounts(const WrappedNodeList& path) = 0; 
    
    virtual void extendDiscounts(const WrappedNodeList& path, 
                                 d_vec& discount_path) = 0;

    virtual d_vec getConcentrations(const WrappedNodeList& path, 
                                    const d_vec& discounts) = 0;

    virtual void extendConcentrations(const WrappedNodeList& path, 
                                      const d_vec& discounts, 
                                      d_vec& concentration_path) = 0;

    virtual double getDiscount(l_type parent_length, l_type this_length) = 0;

    virtual double getConcentration(double discount,
                                    l_type parentLength,
                                    l_type thisLength) = 0;

    virtual double getDiscount(l_type level) = 0;
};

/**
 * A basic implementation of a class that stores per-level 
 * discount and concentration parameters.
 */
class SimpleParameters : public IParameters {
  public:
    
    SimpleParameters();

    SimpleParameters(d_vec s) : discounts(s), alpha(0) {}

    SimpleParameters(d_vec s, double alpha) : discounts(s), alpha(alpha) {}

    /**
     * Get discount parameters for each node in the node list.
     */
    d_vec getDiscounts(const WrappedNodeList& path);

    void extendDiscounts(const WrappedNodeList& path, d_vec& discount_path);

    d_vec getConcentrations(const WrappedNodeList& path, 
                            const d_vec& discounts);
    
    void extendConcentrations(const WrappedNodeList& path, 
                              const d_vec& discounts, 
                              d_vec& concentration_path);

    double getDiscount(l_type parent_length, l_type this_length);
    
    double getConcentration(double discount,
                            l_type parentLength,
                            l_type thisLength);

    double getDiscount(l_type level);
  
  private:
    d_vec discounts;
    double alpha;

    DISALLOW_COPY_AND_ASSIGN(SimpleParameters);
};

}} // namespace gatsby::libplump

#endif
