/**
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
 **/

#ifndef UTILS_H_
#define UTILS_H_

#include <string>
#include <cmath>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <gsl/gsl_sf_gamma.h>

////////////////////////////////////////////////////////////////////////////////
/////////////////   SOME USEFUL MACROS   ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
#define tracer if (0) ; else std::cerr
#define DBG if (1)
#else
#define tracer if (1) ; else std::cerr
#define DBG if (1) ; else
#endif

  // some syntactic sugar for implementing interfaces using (multiple) inheritance
#define interface class
#define implements public

  // A macro to disallow the copy constructor and operator= functions
  // This should be used in the private: declarations for a class
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)


namespace gatsby { namespace libplump {

/**
 * Alias: std::vector<double>
 */
typedef std::vector<double> d_vec;

/**
 * Alias: std::vector<std::vector<double> >
 */
typedef std::vector<d_vec> d_vec_vec;

/**
 * Alias: std::vector<unsigned int>
 */
typedef std::vector<unsigned int> ui_vec;

/**
 * Alias: std::vector<std::vector<unsigned int> >
 */
typedef std::vector<ui_vec> ui_vec_vec;



/**
 * Compute log2(x) -- the log2 provided by GCC is somewhat weird.
 */
inline double log2(double x){
  static const double LOG2 = std::log(2.);
  return std::log(x)/LOG2;
}


/**
 * Read items of type E from a file and push them into a sequence of type
 * S (that has to support push_back(E item). 
 *
 * Items are read from the file sizeof(E) bytes at a time. No error checking
 * is performed.
 *
 * @tparam E Type of element to read from file, elements will be read sizeof(E)
 *           bytes at a time.
 * @tparam S Type of sequence to push items into; 
 *           must support push_back(E item)
 * @param fileName Name of the file to read from 
 * @param seq Sequence to push items into
 * @param limit Maximum number of items to push (default: 0, no limit)
 */
template<typename E, class S>
  inline void pushFileToVec(const std::string& fileName, S& seq, int limit = 0){
    std::ifstream in;
    in.open(fileName.c_str(),std::ios::binary | std::ios_base::in);
    char buffer[sizeof(E)];
    in.read(buffer,sizeof(E));
    int j = 0;
    if (limit==0)
      limit = -1;
    while (!in.eof() && j!=limit){
      seq.push_back(*((E*)buffer));
      in.read(buffer,sizeof(E));
      j++;
    }
    in.close();
  }

/**
 * Push each character of string s individually into container seq of type S, 
 * which has to support push_back(char c);
 *
 * @tparam S type of container to push into
 * @param s input string to push
 * @param seq output container to push into
 */
template<class S>
  inline void pushStringToVec(const std::string& s, S& seq){
    for (size_t i=0;i<s.length();i++) {
      seq.push_back(s[i]);
    }
  }

/**
 * Convert a map<T,T> where T is an integer type with 
 * range 0...MAX to a vector v of length MAX+1 so that 
 * v[x] = m[x] for x=1,...,MAX. v[x] = 0 if m.count(x) = 0.
 *
 * @tparam T integer type
 * @param m input map
 * @param v output vector
 * @returns pointer to histogram vector on the heap
 */
template<typename T>
  inline void mapHistogramToVectorHistogram(const std::map<T,T>& m, std::vector<T>& v) {
    T max = (*max_element(m.begin(), m.end())).first;
    v.clear();
    v.assign(max + 1,0);
    for (typename std::map<T,T>::iterator i = m.begin(); i != m.end(); ++i) {
      v[(*i).first] = (*i).second; 
    }
    return v;
  }

/** 
 * Count the number of times the elements k=0...max-1 occur in vec and 
 * return the result in ret[k].
 */
template<typename T>
  inline std::vector<T> vec2hist(std::vector<T>& vec, T max) {
    std::vector<T> ret(max,0);
    for (unsigned int i = 0; i < vec.size(); ++i) {
      ++ret[vec[i]];
    }
    return ret;
  }

/**
 * Compute the mean of the elements in a sequence.
 *
 * @tparam T iteratable sequence type; must have a const_iterator member
 *           as well as begin() and end() methods.
 * @param in sequence to compute the mean of
 * @returns the mean of the elements in the sequence
 */
template<typename T>
  inline double mean(const T& in) {
    double mean = 0.0;
    size_t j = 0;
    for (typename T::const_iterator i = in.begin(); i != in.end(); ++i) {
      mean += (((double)*i)-mean)/(++j);
    }
    return mean;
  }

/**
 * Sum the elements in a sequence.
 *
 * @tparam T iteratable sequence type; must have a const_iterator member
 *           as well as begin() and end() methods.
 * @param in sequence to compute the sum of
 * @returns the sum of the elements in the sequence
 */
template<typename T>
  inline typename T::value_type sum(const T& in) {
    typename T::value_type sum = 0;
    for (typename T::const_iterator i = in.begin(); i != in.end(); ++i) {
      sum += *i;
    }
    return sum;
  }

inline void log2_vec(d_vec& in) {
  for (size_t i=0;i<in.size();i++){
    in[i] = log2(in[i]);
  }
}

/**
 * Multiply all elements of a vector by a constant.
 */
template<typename elem_t>
  inline void mult_vec(std::vector<elem_t>& in, elem_t mult) {
    for (size_t i=0;i<in.size();i++){
      in[i] = in[i] * mult;
    }
  }

/**
 * Add a constant to all elements of a vector.
 */
template<typename elem_t>
  inline void add_vec(std::vector<elem_t>& in, elem_t add) {
    for (size_t i=0;i<in.size();i++){
      in[i] = in[i] + add;
    }
  }


/**
 * Add two vectors elementwise.
 */
template<typename elem_t>
  inline void add_vec(std::vector<elem_t>& inout, std::vector<elem_t>& add) {
    std::transform(inout.begin(), inout.end(), add.begin(), inout.begin(),
                  std::plus<elem_t>());
  }

/**
 * Elementwise multiplication. Results is returned in the first argument.
 */
template<typename elem_t>
  inline void mult_vec(std::vector<elem_t>& inout, std::vector<elem_t> in) {
    assert(inout.size() == in.size());
    for (size_t i=0;i<inout.size();i++){
      inout[i] *= in[i];
    }
  }

template<typename elem_t>
  inline double prob2loss(std::vector<elem_t> in) {
    double out = 0.0;
    for (size_t i=0;i<in.size();i++){
      out -= log2(in[i]);
    }
    return out/in.size();
  }

inline bool closeTo(double value, double target, double tol=10e-4) {
  return std::abs(value-target) < tol;
}

template<typename Iterable>
  inline std::string iterableToString(const Iterable input) {
    std::ostringstream output;
    for(typename Iterable::const_iterator it = input.begin(); it!=input.end();it++) {
      output << *it << ", ";
    }
    return output.str();
  }

template<typename Iterable>
  void iterableToCSVFile(const Iterable input, std::string fn) {
    std::ofstream output(fn.c_str());
    output.precision(10);
    for(typename Iterable::const_iterator it = input.begin(); it!=input.end();it++) {
      output << *it << ", ";
    }
  }

/**
 * computes log(exp(a) + exp(b)) while avoiding numerical
 * instabilities.
 */
inline double logsumexp(double a, double b) {
  // choose c to be the one that is largest in abs value
  double c = (a>b)?a:b;
  return (log(exp(a-c) + exp(b-c))) + c;
}


inline double sigmoid(double x) {
  return 1.0/(1.0 + exp(-x));
}

// logit = inverse sigmoid
inline double logit(double x) {
  return log(x) - log(1-x);
}


inline double logKramp(double base, double inc, double lim) {
  if (lim <= 0) {
    return 0;
  } else {
    return lim*log(inc) + gsl_sf_lnpoch(base/inc, lim);
  }
}


inline double kramp(double base, double inc, double lim) {
  return exp(logKramp(base, inc, lim));
}


static clock_t global_clock;

/*
 * Start timing
 */
inline void tic() {
  global_clock = clock();
}

/**
 * Return number of seconds since last tic()
 */
inline double toc() {
  return (clock() - global_clock)/(double)CLOCKS_PER_SEC; 

}


inline std::string makeProgressBarString(double percentDone, int total=10) {
  std::ostringstream out;

  int numDone = (int)floor(total * percentDone);

  out << "[";
  for (int i=0; i<numDone; i++) {
    out << "=";
  }
  out << ">";
  for (int i=numDone; i<total; i++) {
    out << " ";
  }
  out << "]";
  return out.str();
}


}} // namespace gatsby::libplump

#endif /* UTILS_H_ */
