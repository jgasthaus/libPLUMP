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

#include <limits>

#include "libplump/config.h"
#include "libplump/context_tree.h" // for WrappedNodeList
#include "libplump/utils.h" // for WrappedNodeList
#include "libplump/hpyp_parameters_interface.h"

namespace gatsby { namespace libplump {


/**
 * Compute the gradient of the PYP predictive probability
 * with respect to the discount parameter.
 */
double PYPPredictiveGradientDiscount(int cw, int tw, int c, int t, 
    double parentProbability, double discount, double concentration,
    double parentGradientDiscount) {

    return   (-tw + parentProbability*t + (discount*t + concentration)*parentGradientDiscount)
           / (c + concentration);

}


double PYPPredictiveGradientIndividualDiscount(int cw, int tw, int c, int t, 
      double parentProbability, double totalDiscount, double individualDiscount, double concentration,
      double DParentProbabilityWrtIndividualDiscount = 0.0,
      int individualDiscountMultiplicity = 1) 
{
    double DTotalDiscountWrtIndividualDiscount = individualDiscountMultiplicity * totalDiscount / individualDiscount;
    double DConcentrationWrtIndividualDisount = concentration / individualDiscount;
    double denom = (concentration + c);
    double denom2 = denom * denom;
    double term1 = DConcentrationWrtIndividualDisount      * (c*parentProbability + totalDiscount*(tw - t*parentProbability) - cw);
    double term2 = DTotalDiscountWrtIndividualDiscount     * (-tw + parentProbability*t);
    double term3 = DParentProbabilityWrtIndividualDiscount * (totalDiscount*t + concentration);
    return term1/denom2 + (term2 + term3)/denom; // the last factor is te same as above
}


/**
 * Compute the gradient of the PYP predictive probability
 * with respect to the (underlying) concentration parameter \alpha_0.
 * The given concentration parameter is assumed to be the actual concentration
 * parameter for the predictive distribution, so that
 *   concentration = \alpha_0 * discount
 */
double PYPPredictiveGradientConcentration(int cw, int tw, int c, int t, 
      double parentProbability, double discount, double concentration,
      double DParentProbabilityWrtConcentration) 
{
    double DConcentrationWrtAlpha0 = discount;
    double denom = (concentration + c);
    double denom2 = denom * denom;
    double term1 = (-cw + tw*discount + (c - t*discount)*parentProbability)*DConcentrationWrtAlpha0 / denom2;
    double term2 = (concentration + t*discount)*DParentProbabilityWrtConcentration / denom;
    return term1 + term2;
}


/**
 * A basic implementation of a class that stores per-level 
 * discount and concentration parameters.
 */
class SimpleParameters : public IParameters {
  public:
    
    SimpleParameters();

    SimpleParameters(d_vec s) : discounts(s), alpha(0) {}

    SimpleParameters(d_vec s, double alpha) : discounts(s), alpha(alpha)  {}

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
    
    void accumulateParameterGradient(
      const IAddRemoveRestaurant& restaurant,
      const WrappedNodeList& path, 
      const d_vec& prob_path, 
      const d_vec& discount_path, 
      const d_vec& concentration_path, 
      e_type obs);
    
    void stepParameterGradient(double stepSize);
  
    d_vec discounts;
    double alpha;

  private:
    DISALLOW_COPY_AND_ASSIGN(SimpleParameters);
};


/**
 * A basic implementation of a class that stores per-level 
 * discount and concentration parameters. Stores the parameters
 * in a form that is suitable for unconstraint optimization via
 * stochastic gradient descent.
 */
class GradientParameters : public IParameters {
  public:
    
    GradientParameters() : sigmoid_discounts(), 
        log_alpha(-std::numeric_limits<double>::infinity()),
        sigmoid_discount_gradient(), log_alpha_gradient(0) {
      sigmoid_discounts.push_back(logit(0.5));
      sigmoid_discount_gradient.push_back(0);
    }

    GradientParameters(d_vec s) : sigmoid_discounts(s.size()),
                                  log_alpha(-std::numeric_limits<double>::infinity()),
                                  sigmoid_discount_gradient(s.size(),0.0),
                                  log_alpha_gradient(0.0) {
      this->setDiscounts(s);
    }

    GradientParameters(d_vec s, double alpha) : 
        sigmoid_discounts(s.size()), 
        log_alpha(log(alpha)),
        sigmoid_discount_gradient(s.size(),0.0),
        log_alpha_gradient(0.0) {
      this->setDiscounts(s);
    }

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
    
    void accumulateParameterGradient(
      const IAddRemoveRestaurant& restaurant,
      const WrappedNodeList& path, 
      const d_vec& prob_path, 
      const d_vec& discount_path, 
      const d_vec& concentration_path, 
      e_type obs);
    
    void stepParameterGradient(double stepSize);

    std::pair<d_vec, double> approximateParameterGradient(
        const IAddRemoveRestaurant& restaurant,
        const WrappedNodeList& path, 
        const d_vec& prob_path, 
        const d_vec& discount_path, 
        const d_vec& concentration_path, 
        e_type obs);
  
  private:

    // set the internal discounts given "real discounts"
    void setDiscounts(d_vec& d) {
      for (size_t i = 0; i < d.size(); ++i) {
        sigmoid_discounts[i] = logit(d[i]);
      }
    }

    d_vec sigmoid_discounts;
    double log_alpha;
    d_vec sigmoid_discount_gradient;
    double log_alpha_gradient;
    static const int mini_batch_size = 100;
    static const double default_step_size = 1e-5;

    DISALLOW_COPY_AND_ASSIGN(GradientParameters);
};



}} // namespace gatsby::libplump

#endif
