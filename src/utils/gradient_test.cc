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

#include <iostream>
#include <fstream>
#include <libgen.h>
#include <cmath>
#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/scoped_ptr.hpp>

#include <libplump/libplump.h>
#include <libplump/switching_restaurant.h>

#define DISCOUNT sigmoid(d0)*sigmoid(d1)
#define CONCENTRATION sigmoid(d0)*sigmoid(d1)*exp(a)

using namespace std;
using namespace gatsby::libplump;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

static unsigned int num_types = 2;


IParameters* getParameters(po::variables_map& vm) {
  return new SimpleParameters(vm["disc"].as<d_vec>(), vm["alpha"].as<double>());
}



double testSingleGradients(po::variables_map& vm) {
  seq_type seq;

  boost::scoped_ptr<IParameters> parameters(new SimpleParameters(vm["disc"].as<d_vec>(), vm["alpha"].as<double>()));
  boost::scoped_ptr<IAddRemoveRestaurant> restaurant(new HistogramRestaurant());
  void* payload = restaurant->getFactory().make();

  double d0 = logit(0.5);
  double d1 = logit(0.5);
  double a = 0.;
  double eps = 10e-7;

  for (int i = 0; i < 100; ++i) {
    restaurant->addCustomer(payload, 0, 1/2., DISCOUNT, CONCENTRATION);
    restaurant->addCustomer(payload, 1, 1/2., DISCOUNT, CONCENTRATION);
    restaurant->addCustomer(payload, 1, 1/2., DISCOUNT, CONCENTRATION);
  }
  cout << "c0: " <<   restaurant->getC(payload, 0)
       << ", t0: " << restaurant->getT(payload, 0)
       << ", c: "  << restaurant->getC(payload)
       << ", t: "  << restaurant->getT(payload)
       << endl;

  for (int i = 0; i < 100; ++i) {

    double p0 =  restaurant->computeProbability(payload, 0, 1/2., DISCOUNT, CONCENTRATION);
    double p1 =  restaurant->computeProbability(payload, 1, 1/2., DISCOUNT, CONCENTRATION);

    cout << "disc: " << DISCOUNT<< ", alpha: " << CONCENTRATION<< ", prob of 0: " << p0 << ", prob of 1: " << p1 << endl;

    double gd0 = PYPPredictiveGradientIndividualDiscount(restaurant->getC(payload, 0),
                                          restaurant->getT(payload, 0),
                                          restaurant->getC(payload),
                                          restaurant->getT(payload),
                                          1/2.,
                                          DISCOUNT,
                                          sigmoid(d0),
                                          CONCENTRATION,
                                          0.0)*sigmoid(d0)*(1-sigmoid(d0))/p0;
    double gd1 = PYPPredictiveGradientIndividualDiscount(restaurant->getC(payload, 1),
                                          restaurant->getT(payload, 1),
                                          restaurant->getC(payload),
                                          restaurant->getT(payload),
                                          1/2.,
                                          DISCOUNT,
                                          sigmoid(d0),
                                          CONCENTRATION,
                                          0.0)*sigmoid(d0)*(1-sigmoid(d0))/p1;
    double ga0 = PYPPredictiveGradientConcentration(restaurant->getC(payload, 0),
                                          restaurant->getT(payload, 0),
                                          restaurant->getC(payload),
                                          restaurant->getT(payload),
                                          1/2.,
                                          DISCOUNT,
                                          CONCENTRATION,
                                          0.0)*exp(a)/p0;
    double ga1 = PYPPredictiveGradientConcentration(restaurant->getC(payload, 1),
                                          restaurant->getT(payload, 1),
                                          restaurant->getC(payload),
                                          restaurant->getT(payload),
                                          1/2.,
                                          DISCOUNT,
                                          CONCENTRATION,
                                          0.0)*exp(a)/p1;
    double p0p =  log(restaurant->computeProbability(payload, 0, 1/2., sigmoid(d0+eps)*sigmoid(d1),sigmoid(d0+eps)*sigmoid(d1)*exp(a) ));
    double p0m =  log(restaurant->computeProbability(payload, 0, 1/2., sigmoid(d0-eps)*sigmoid(d1),sigmoid(d0-eps)*sigmoid(d1)*exp(a) ));
    double gd0_est = (p0p - p0m)/(2*eps); 
    p0p = log(restaurant->computeProbability(payload, 0, 1/2., DISCOUNT, DISCOUNT*exp(a+eps)));
    p0m = log(restaurant->computeProbability(payload, 0, 1/2., DISCOUNT, DISCOUNT*exp(a-eps)));
    double ga0_est = (p0p - p0m)/(2*eps); 
    cout << "Discount: grad of 0: " << gd0
         << ", estimate: " << gd0_est
         << ", grad of 1: " << gd1 << endl 
         << "Alpha grad of 0: " << ga0
         << ", estimate: " << ga0_est
         << ", grad of 1: " << ga1 
         << endl;
    //restaurant->addCustomer(payload, 0, 1/2., 0.5, 1);
    //restaurant->addCustomer(payload, 0, 1/2., sigmoid(d), exp(a));
    d0 += gd0;
    a += ga0;
  }

  restaurant->getFactory().recycle(payload);

}


int main(int argc, char* argv[]) {
  const double sm_disc[] = {.62, .69, .74, .80, .95};
  d_vec default_discounts;
  default_discounts.assign(sm_disc,&sm_disc[5]);
  // Declare the supported options.
  po::options_description generic("Generic options");
  generic.add_options()
    ("alpha,a", po::value<double>()->default_value(5), "Concentration parameter") 
    ("disc,d", po::value<d_vec>()->default_value(default_discounts,"..."), "Discount parameter(s)") 
    ("num-types", po::value<int>()->default_value(256), "Number of types")
    ;


  po::options_description cmdline_options;
  cmdline_options.add(generic);

  po::options_description visible_options;
  visible_options.add(generic);


  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
      options(cmdline_options).run(), vm);
  po::notify(vm);

  string usage = "Usage: gradient_test [OPTIONS]...\n";


  num_types = vm["num-types"].as<int>();

  init_rng();

  testSingleGradients(vm);

  free_rng();
}
