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

#include <cmath> // for pow
#include <algorithm> // for fill
#include "libplump/hpyp_parameters.h"
#include "libplump/context_tree.h" // for WrappedNodeList
#include "libplump/utils.h" // for tracer


namespace gatsby { namespace libplump {


SimpleParameters::SimpleParameters() : discounts(), alpha(0) {
  discounts.push_back(0.5);
}

/**
 * Get discount parameters for each node in the node list.
 */
d_vec SimpleParameters::getDiscounts(const WrappedNodeList& path) {
  d_vec discount_path;
  int parent_length = -1;
  int max_context_len = discounts.size()-1;
  for(WrappedNodeList::const_iterator it = path.begin(); 
      it != path.end(); ++it) {
    double discount = 1;
    int this_length = it->end - it->start;
    
    // we need to loop over all depths that lie between parent_length+1 and
    // this_length (inclusive)
    for(int i = parent_length + 1; i <= std::min(this_length,max_context_len);
        ++i) {
      discount *= discounts[i];
    }
    
    if(this_length > max_context_len) { 
      // our context is longer than we have explicit discounts
      int additional_depth =   this_length 
                             - std::max(parent_length,max_context_len);
      discount *= std::pow(discounts.back(),additional_depth); 
    }
    
    discount_path.push_back(discount);
    parent_length = this_length;
  }
  return discount_path;
}


void SimpleParameters::extendDiscounts(const WrappedNodeList& path, 
                                       d_vec& discount_path) {
  int max_context_len = discounts.size()-1;
  WrappedNodeList::const_iterator it = path.begin();
  for (int i = 0; i < (int)discount_path.size() - 1; ++i) {
    it++;
  }
  int parent_length = (it->end - it->start); 
  it++;
  for(; it!= path.end(); it++) {
    double discount = 1;
    int this_length = it->end - it->start;
    
    // we need to loop over all depths that lie between parent_length+1 and 
    // this_length (inclusive)
    for(int i = parent_length+1; 
        i <= std::min(this_length,max_context_len);
        ++i) {
      discount *= discounts[i];
    }

    if(this_length > max_context_len) { 
      // our context is longer than we have explicit discounts
      int additional_depth =   this_length 
                             - std::max(parent_length,max_context_len);
      discount *= std::pow(discounts.back(),additional_depth); 
    }

    discount_path.push_back(discount);
    parent_length = this_length;
  }
}


d_vec SimpleParameters::getConcentrations(const WrappedNodeList& path, 
                                          const d_vec& discounts) {
  d_vec concentration_path;
  double current = alpha;
  for (d_vec::const_iterator it = discounts.begin();
      it != discounts.end(); ++it) {
    concentration_path.push_back(current);
    current *= *it;
  }
  return concentration_path;
}


void SimpleParameters::extendConcentrations(const WrappedNodeList& path,
                                            const d_vec& discounts,
                                            d_vec& concentration_path) {
  double current = concentration_path.back();
  for (int i = concentration_path.size(); i < (int)discounts.size(); i++) {
    current *= discounts[i];
    concentration_path.push_back(current);
  }
}


double SimpleParameters::getDiscount(int parent_length, 
                                     int this_length) {
  tracer << "SimpleParameters::getDiscount(" << parent_length << ", " 
         << this_length << ")" << std::endl;
  int max_context_len = discounts.size()-1;
  double discount = 1;

  // we need to loop over all depths that lie between parent_length+1 and 
  // this_length (inclusive)
  for(int i = parent_length+1;
      i <= std::min(this_length,max_context_len);
      i++) {
    discount *= discounts[i];
  }

  if(this_length > max_context_len) { 
    // our context is longer than we have explicit discounts
    int additional_depth =   this_length
                           - std::max(parent_length,max_context_len);
    discount *= std::pow(discounts.back(),additional_depth); 
  }

  return discount;
}
    

double SimpleParameters::getConcentration(double discount,
                                           l_type parentLength,
                                           l_type thisLength) {
  return this->alpha * discount;
}


double SimpleParameters::getDiscount(int level) {
  if (level < (int)discounts.size()) {
    return discounts[level];
  } else {
    return discounts.back();
  }
}


void SimpleParameters::accumulateParameterGradient(
    const IAddRemoveRestaurant& restaurant,
    const WrappedNodeList& path, 
    const d_vec& prob_path, 
    const d_vec& discount_path, 
    const d_vec& concentration_path, 
    e_type obs) {
  // DO NOTHING FOR NOW
}
    

void SimpleParameters::stepParameterGradient(double stepSize) {
  // DO NOTHING FOR NOW
}



/////// GradientParameters //////////////////////////

/**
 * Get discount parameters for each node in the node list.
 */
d_vec GradientParameters::getDiscounts(const WrappedNodeList& path) {
  d_vec discount_path;
  int parent_length = -1;
  int max_context_len = sigmoid_discounts.size()-1;
  for(WrappedNodeList::const_iterator it = path.begin(); 
      it != path.end(); ++it) {
    int this_length = it->end - it->start;
    
    // we need to loop over all depths that lie between parent_length+1 and
    // this_length (inclusive)
    double discount = 1; 
    for(int i = parent_length + 1; i <= std::min(this_length,max_context_len);
        ++i) {
      discount *= sigmoid(sigmoid_discounts[i]);
    }
    
    if(this_length > max_context_len) { 
      // our context is longer than we have explicit discounts
      int additional_depth =   this_length 
                             - std::max(parent_length,max_context_len);
      discount *= std::pow(sigmoid(sigmoid_discounts.back()),additional_depth); 
    }
    
    discount_path.push_back(discount);
    parent_length = this_length;
  }
  return discount_path;
}


void GradientParameters::extendDiscounts(const WrappedNodeList& path, 
                                       d_vec& discount_path) {
  int max_context_len = sigmoid_discounts.size()-1;
  WrappedNodeList::const_iterator it = path.begin();
  for (int i = 0; i < (int)discount_path.size() - 1; ++i) {
    it++;
  }
  int parent_length = (it->end - it->start); 
  it++;
  for(; it!= path.end(); it++) {
    int this_length = it->end - it->start;
    
    // we need to loop over all depths that lie between parent_length+1 and 
    // this_length (inclusive)
    double discount = 1;
    for(int i = parent_length+1; 
        i <= std::min(this_length,max_context_len);
        ++i) {
      discount *= sigmoid(sigmoid_discounts[i]);
    }

    if(this_length > max_context_len) { 
      // our context is longer than we have explicit discounts
      int additional_depth =   this_length 
                             - std::max(parent_length,max_context_len);
      discount *= std::pow(sigmoid(sigmoid_discounts.back()),additional_depth); 
    }

    discount_path.push_back(discount);
    parent_length = this_length;
  }
}


d_vec GradientParameters::getConcentrations(const WrappedNodeList& path, 
                                          const d_vec& discounts) {
  d_vec concentration_path;
  double current = exp(log_alpha);
  for (d_vec::const_iterator it = discounts.begin();
      it != discounts.end(); ++it) {
    concentration_path.push_back(current);
    current *= *it;
  }
  return concentration_path;
}


void GradientParameters::extendConcentrations(const WrappedNodeList& path,
                                            const d_vec& discounts,
                                            d_vec& concentration_path) {
  double current = concentration_path.back();
  for (int i = concentration_path.size(); i < (int)discounts.size(); i++) {
    current *= discounts[i];
    concentration_path.push_back(current);
  }
}


double GradientParameters::getDiscount(int parent_length, 
                                     int this_length) {
  tracer << "GradientParameters::getDiscount(" << parent_length << ", " 
         << this_length << ")" << std::endl;
  int max_context_len = sigmoid_discounts.size()-1;

  // we need to loop over all depths that lie between parent_length+1 and 
  // this_length (inclusive)
  double discount = 1;
  for(int i = parent_length+1;
      i <= std::min(this_length,max_context_len);
      i++) {
    discount *= sigmoid(sigmoid_discounts[i]);
  }

  if(this_length > max_context_len) { 
    // our context is longer than we have explicit discounts
    int additional_depth =   this_length
                           - std::max(parent_length,max_context_len);
    discount *= std::pow(sigmoid(sigmoid_discounts.back()),additional_depth); 
  }

  return discount;
}
    

double GradientParameters::getConcentration(double discount,
                                           l_type parentLength,
                                           l_type thisLength) {
  return exp(this->log_alpha) * discount;
}


double GradientParameters::getDiscount(int level) {
  if (level < (int)sigmoid_discounts.size()) {
    return sigmoid(sigmoid_discounts[level]);
  } else {
    return sigmoid(sigmoid_discounts.back());
  }
}


void GradientParameters::accumulateParameterGradient(
    const IAddRemoveRestaurant& restaurant,
    const WrappedNodeList& path, 
    const d_vec& prob_path, 
    const d_vec& discount_path, 
    const d_vec& concentration_path, 
    e_type obs) {
  // XXX: THERE IS A BUG IN THE GRADIENT COMPUTATION SOMEWHERE, 
  //      SO FALLING BACK TO FINITE DIFFERENCES FOR NOW
  // int parent_length = -1;
  // int max_context_len = sigmoid_discounts.size()-2;
  // double last_discount_gradient = 0;
  // int j = 0;
  // for(WrappedNodeList::const_iterator it = path.begin(); 
  //     it != path.end(); ++it) {
  //   int cw = restaurant.getC(it->payload, obs);
  //   int tw = restaurant.getT(it->payload, obs);
  //   int c  = restaurant.getC(it->payload);
  //   int t  = restaurant.getT(it->payload);
  //   int this_length = it->end - it->start;
  //   for(int i = parent_length + 1; i <= std::min(this_length, max_context_len); ++i) {
  //     double s = sigmoid(this->sigmoid_discounts[i]);
  //     // update explicit discount gradients
  //     this->sigmoid_discount_gradient[i] += PYPPredictiveGradientIndividualDiscount(cw, tw, c, t, 
  //     prob_path[j], discount_path[j], s, 
  //     concentration_path[j]) / prob_path[j+1] * s * (1 - s);
  //   }
  //   if(this_length > max_context_len) { 
  //      // our context is longer than we have explicit discounts
  //     int additional_depth =   this_length 
  //                              - std::max(parent_length, max_context_len);
  //     double s = sigmoid(this->sigmoid_discounts.back());
  //     // update explicit discount gradients
  //     last_discount_gradient = PYPPredictiveGradientIndividualDiscount(
  //               cw, tw, c, t, 
  //               prob_path[j], discount_path[j], s, concentration_path[j],
  //               last_discount_gradient, additional_depth
  //               ) / prob_path[j+1] * s * (1 - s);
  //     this->sigmoid_discount_gradient[this->sigmoid_discount_gradient.size() - 1] += last_discount_gradient;
  //   }
  //   //  
  //   parent_length = this_length;
  //   ++j;
  // }
    std::pair<d_vec, double> approx_grad = this->approximateParameterGradient(
      restaurant,
      path, 
      prob_path, 
      discount_path, 
      concentration_path, 
      obs);


  // std::cout << "Grad: " << iterableToString(this->sigmoid_discount_gradient) << std::endl; 
  // std::cout << "Approx-grad: " << iterableToString(approx_grad) << std::endl;

  fill(this->sigmoid_discount_gradient.begin(), this->sigmoid_discount_gradient.end(), 0.0); 
  add_vec(this->sigmoid_discount_gradient, approx_grad.first);
  this->log_alpha_gradient += approx_grad.second;
}

std::pair<d_vec, double> GradientParameters::approximateParameterGradient(
    const IAddRemoveRestaurant& restaurant,
    const WrappedNodeList& path, 
    const d_vec& prob_path, 
    const d_vec& discount_path, 
    const d_vec& concentration_path, 
    e_type obs) {
  const double eps = 10e-7;
  d_vec ret(this->sigmoid_discounts.size(), 0.0);

  for (int i=0; i < this->sigmoid_discounts.size(); ++i) {
    double orig_discount = sigmoid_discounts[i];
    this->sigmoid_discounts[i] += eps;
    d_vec d_path = this->getDiscounts(path);
    d_vec c_path = this->getConcentrations(path, d_path);
    int j = 0;
    double pp = 1./256.;
    double pm = pp;
    for(WrappedNodeList::const_iterator it = path.begin(); 
        it != path.end(); ++it) {
      //std::cout << pp << std::endl;
      pp = restaurant.computeProbability(it->payload, 
                                                 obs,
                                                 pp,
                                                 d_path[j],
                                                 c_path[j]);
    ++j;
    }
    this->sigmoid_discounts[i] = orig_discount - eps;
    d_path = this->getDiscounts(path);
    c_path = this->getConcentrations(path, d_path);
    j = 0;
    for(WrappedNodeList::const_iterator it = path.begin(); 
        it != path.end(); ++it) {
      pm = restaurant.computeProbability(it->payload, 
                                                 obs,
                                                 pm,
                                                 d_path[j],
                                                 c_path[j]);
    ++j;
    }
    ret[i] = (log(pp) - log(pm))/(2*eps);
    this->sigmoid_discounts[i] = orig_discount;
  }

  // compute approx gradient wrt. log-alpha
  double orig_log_alpha = this->log_alpha;
  this->log_alpha += eps;
  d_vec d_path = this->getDiscounts(path);
  d_vec c_path = this->getConcentrations(path, d_path);
  int j = 0;
  double pp = 1./256.;
  double pm = pp;
  for(WrappedNodeList::const_iterator it = path.begin(); 
      it != path.end(); ++it) {
    pp = restaurant.computeProbability(it->payload, 
                                               obs,
                                               pp,
                                               d_path[j],
                                               c_path[j]);
  ++j;
  }
  this->log_alpha = orig_log_alpha - eps;
  d_path = this->getDiscounts(path);
  c_path = this->getConcentrations(path, d_path);
  j = 0;
  for(WrappedNodeList::const_iterator it = path.begin(); 
      it != path.end(); ++it) {
    pm = restaurant.computeProbability(it->payload, 
                                               obs,
                                               pm,
                                               d_path[j],
                                               c_path[j]);
  ++j;
  }
  double log_alpha_gradient = (log(pp) - log(pm))/(2*eps);
  this->log_alpha = orig_log_alpha;

  return make_pair(ret, log_alpha_gradient);
}

void GradientParameters::stepParameterGradient(double stepSize) {
  std::cout << "Grad: " << iterableToString(this->sigmoid_discount_gradient) << std::endl; 
  mult_vec(this->sigmoid_discount_gradient, stepSize);
  add_vec(this->sigmoid_discounts, this->sigmoid_discount_gradient);
  fill(this->sigmoid_discount_gradient.begin(), this->sigmoid_discount_gradient.end(), 0.0); 
  std::cout << "Disc: " << iterableToString(this->sigmoid_discounts) << std::endl; 

  this->log_alpha += stepSize * this->log_alpha_gradient;
  this->log_alpha_gradient = 0.0;

}


}} // namespace gatsby::libplump
