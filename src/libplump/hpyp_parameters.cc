#include "libplump/hpyp_parameters.h"
#include <cmath> // for pow
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
  static int num_gradients = 0;
  if (num_gradients == mini_batch_size) {
    this->stepParameterGradient(default_step_size);
    num_gradients = 0;
  }
  
  num_gradients++;
}
    

void SimpleParameters::stepParameterGradient(double stepSize) {
  for (size_t i = 0; i < this->discounts.size(); ++i) {
    discounts[i] += stepSize * discount_gradient[i];
    discount_gradient[i] = 0.0;
  }
  alpha += stepSize * alpha_gradient;
  alpha_gradient = 0.0;
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
  // do nothing
}
    

void GradientParameters::stepParameterGradient(double stepSize) {
  // do nothing
}


}} // namespace gatsby::libplump
