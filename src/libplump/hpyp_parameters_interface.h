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
 * The classes in this file represent and store the parameters for the
 * HPYP models (discount and concentration parameters).
 **/

#ifndef HPYP_PARAMETERS_INTERFACE_H_
#define HPYP_PARAMETERS_INTERFACE_H_

#include "libplump/context_tree.h" // for WrappedNodeList
#include "libplump/hpyp_restaurant_interface.h" // for IAddRemoveRestaurant

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

    virtual void accumulateParameterGradient(
      const IAddRemoveRestaurant& restaurant,
      const WrappedNodeList& path, 
      const d_vec& prob_path, 
      const d_vec& discount_path, 
      const d_vec& concentration_path, 
      e_type obs) = 0;

    virtual void stepParameterGradient(double stepSize) = 0;
};

}} // namespace gatsby::libplump

#endif
