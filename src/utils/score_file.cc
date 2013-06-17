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

using namespace std;
using namespace gatsby::libplump;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

static unsigned int num_types = 256;

/**
 * Predict probabilities on test file.
 */
d_vec predict(po::variables_map& vm, HPYPModel& m, int start_pos, seq_type seq) {
  d_vec predictive;
  for(int i = start_pos + 1; i < (int)seq.size(); ++i) {    
    if (vm.count("sum")) {
        d_vec dist = m.predictiveDistribution(start_pos, i);
        if (!closeTo(sum(dist), 1)) {
          std::cout << "probs don't sum to one: " << iterableToString(dist) << std::endl;
        }
    }
    switch (vm["fragment"].as<int>()) {
      case 1:
        predictive.push_back(m.predict(start_pos, i, seq[i]));
        break;
      case 2:
        predictive.push_back(m.predictWithFragmentation(start_pos, i, seq[i]));
        break;
      case 3:
        predictive.push_back(m.predictBelow(start_pos, i, seq[i]));
        break;
    }
  }
  return predictive;
}

IParameters* getParameters(po::variables_map& vm) {
  switch(vm["parameters"].as<int>()) {
    case 0:
      cerr << "getParameters(): Using SimpleParameters" << endl;
      return new SimpleParameters(vm["disc"].as<d_vec>(), vm["alpha"].as<double>());
    case 1:
      cerr << "getParameters(): Using GradientParameters" << endl;
      return new GradientParameters(vm["disc"].as<d_vec>(), vm["alpha"].as<double>());
  }
}

IAddRemoveRestaurant* getRestaurant(po::variables_map& vm) {
  switch(vm["restaurant"].as<int>()) {
    case 0:
      cerr << "getRestaurant(): Using KneserNeyRestaurant" << endl;
      return new KneserNeyRestaurant();
    case 1:
      cerr << "getRestaurant(): Using SimpleFullRestaurant" << endl;
      return new SimpleFullRestaurant();
    case 2:
      cerr << "getRestaurant(): Using HistogramRestaurant" << endl;
      return new HistogramRestaurant();
    case 3:
      cerr << "getRestaurant(): Using ReinstantiatingCompactRestaurant" << endl;
      return new ReinstantiatingCompactRestaurant();
    case 4:
      cerr << "getRestaurant(): Using StirlingCompactRestaurant" << endl;
      return new StirlingCompactRestaurant();
    case 5:
      cerr << "getRestaurant(): "
           << "Using SwitchingRestaurant(SimpleFullRestaurant, 10)"
           << endl;
      return new SwitchingRestaurant(new SimpleFullRestaurant(), 10);
    case 6:
      cerr << "getRestaurant(): Using PowerLawRestaurant" << endl;
      return new PowerLawRestaurant();
    case 7:
      cerr << "getRestaurant(): Using FractionalRestaurant" << endl;
      return new FractionalRestaurant();
    case 8:
      cerr << "getRestaurant(): Using ExpectedTablesCompactRestaurant" << endl;
      return new ExpectedTablesCompactRestaurant();
  }
  cout << "Unknown restaurant type (--restaurant)!";
  exit(1);
}

INodeManager* getNodeManager(po::variables_map& vm,
    const IPayloadFactory& payloadFactory) {
  return new SimpleNodeManager(payloadFactory);
}

void pushFileToSeq(po::variables_map& vm, std::string filename, seq_type& seq) {
  if (vm.count("read-int32")) {
    pushFileToVec<int>(filename, seq, vm["head"].as<int>());
  } else {
    pushFileToVec<unsigned char>(filename, seq, vm["head"].as<int>());
  }
}

void runSampler(po::variables_map& vm, HPYPModel& model, int train_length) {
  if (vm["sampler"].as<int>() == 1) {
    model.runGibbsSampler();
  } else {
    model.removeAddSweep(1, train_length);
  }
}


double score_file(po::variables_map& vm) {
  string filename = vm["input-file"].as<string>();
  fs::path input_path(filename);
  if (!fs::exists(input_path)) {
    cerr << "Input file not found!" << endl;
    exit(1);
  }

  seq_type seq;
  pushFileToSeq(vm, filename, seq);
  cout << "Sequence length: " << seq.size() << endl;
  cout << "Number of types:  " << num_types << endl;
  cout << seq[0] << endl;

  boost::scoped_ptr<IParameters> parameters(getParameters(vm));
  boost::scoped_ptr<IAddRemoveRestaurant> restaurant(getRestaurant(vm));
  boost::scoped_ptr<INodeManager> nodeManager(
      getNodeManager(vm, restaurant->getFactory()));

  HPYPModel model(seq, *nodeManager, *restaurant, *parameters, num_types);

  d_vec losses;
  if (vm.count("load-serialized-nodes")) {
    Serializer nodeSerializer(vm["load-serialized-nodes"].as<string>());
    nodeSerializer.loadNodesAndPayloads(*nodeManager, restaurant->getFactory());
  } else {
    int lag = vm["lag"].as<int>();
    if (lag == 0) {
      losses = model.computeLosses(0, seq.size());
    } else {
      losses = model.computeLossesWithDeletion(0, seq.size(), lag);
    }
  }
  
  if (vm.count("debug")) {
    model.checkConsistency();
  }

  cout << "Training loss: " << mean(losses) << endl;

  if (vm.count("print-tree")) {
    cout << model.toString() << endl;
  }


  if (vm.count("dump-losses")) {
    iterableToCSVFile(losses,"losses");
  }


  if (vm.count("test-file")) {
    d_vec current_sample_losses;
    int start_pos = seq.size();
    pushFileToSeq(vm, vm["test-file"].as<string>(), seq);
    cout << "Test sequence length: " << seq.size()-start_pos << endl;
    int mode = vm["mode"].as<int>();
    ostringstream dump_fn;
    dump_fn << "losses_" << fs::basename(input_path); // <<  "_"  << m.parameters.alpha;
    dump_fn << "_" << vm["restaurant"].as<int>()  << "_" << vm["mode"].as<int>() << "_" << vm["fragment"].as<int>() << "_" << vm["prefix"].as<string>() << ".csv";
    ofstream losses_f(dump_fn.str().c_str());
    current_sample_losses.push_back(prob2loss<double>(predict(vm, model, start_pos, seq)));
    losses_f << current_sample_losses.back() << ", ";
    cout << "loss: " << current_sample_losses.back() << endl;
    if (vm.count("burn-in")) {
      for (int i = 0; i < vm["burn-in"].as<int>(); ++i) {
        cout << "Burn-in iteration: " << i << endl;
        runSampler(vm, model, start_pos);
        if (vm.count("debug")) {
          if (model.checkConsistency()) {
            cout << "Model consistent." << endl;
          } else {
            cout << "Consitency check failed." << endl;
          }
        }
        current_sample_losses.push_back(prob2loss<double>(predict(vm, model, start_pos, seq)));
        losses_f << current_sample_losses.back() << ", ";
        cout << "loss: " << current_sample_losses.back() << endl;
        losses_f.flush();
        if (vm.count("debug") && vm.count("print-tree")) {
          cout << model.toString() << endl;
        }
      }
      d_vec_vec sample_predictions;
      for (int i = 0; i < vm["samples"].as<int>(); ++i) {
        cout << "sampling iteration iteration: " << i << endl;
        runSampler(vm, model, start_pos);
        sample_predictions.push_back(predict(vm, model, start_pos, seq));
        cout << "loss (this sample): " << prob2loss<double>(sample_predictions.back()) << endl;
        cout << "loss (avg): " << prob2loss<double>(average(sample_predictions)) << endl;
      }
    }
  }

  if (vm.count("save-serialized-nodes")) {
    Serializer nodeSerializer(vm["save-serialized-nodes"].as<string>());
    nodeSerializer.saveNodesAndPayloads(*nodeManager, restaurant->getFactory());
  }
  cout << "Discounts: ";
  for (int i=0; i < vm["disc"].as<d_vec>().size(); ++i) {
    cout << parameters->getDiscount(i) << ", ";
  }
  cout << endl;


  return mean(losses);
}


int main(int argc, char* argv[]) {
  const double sm_disc[] = {0.05, 0.7, 0.8, 0.82, 0.84, 0.88, 0.91, 0.92, 0.93, 0.94, 0.95};
  //const double sm_disc[] = {.62, .69, .74, .80, .95};
  d_vec default_discounts;
  default_discounts.assign(sm_disc,&sm_disc[11]);
  // Declare the supported options.
  po::options_description generic("Generic options");
  po::options_description dumping("Dumping options");
  po::options_description hidden("Hidden options");
  generic.add_options()
    ("help", "produce help message")
    ("debug,D", "Print debugging output")
    ("sum,s", "Check that probabilities sum to one")
    ("print-tree", "Print the context tree to the screen")
    ("fragment", po::value<int>()->default_value(1), "1: nofrag; 2: frag; 3:below")
    ("read-int32", "Read input data as 32 bit integers")
    ("test-file", po::value<string>(), "Test file")
    ("save-serialized-nodes", po::value<string>(), "File to contain serialized nodes")
    ("load-serialized-nodes", po::value<string>(), "File to contain serialized nodes")
    ("head",po::value<int>()->default_value(0), "If given, cuts input to this number of symbols")
    ("mode", po::value<int>()->default_value(1), "1: particle filter, 2: no fragment, 3: fragment")
    ("sampler", po::value<int>()->default_value(1), "1: gibbs, 2: remove-add")
    ("restaurant", po::value<int>()->default_value(1),
     "0:KN, 1: SimpleFull, 2: Histogram, 3: ReinstantiatingCompact, 4: StirlingCompact, 5: Switching, 6: PowerLaw, 7: Fractional")
    ("parameters", po::value<int>()->default_value(0),
     "0:Simple, 1: Gradient")
    ("burn-in",po::value<int>()->default_value(0), "Number of Gibbs iterations for burn in")
    ("samples,s",po::value<int>()->default_value(1), "Number of samples used for prediction")
    ("lag,l",po::value<int>()->default_value(0), "Lag for deleted prediction (0=off)")
    ("num-types", po::value<int>()->default_value(256), "Number of types") 
    ("alpha,a", po::value<double>()->default_value(5), "Concentration parameter") 
    ("disc,d", po::value<d_vec>()->default_value(default_discounts,"..."), "Discount parameter(s)") 
    ;

  dumping.add_options()
    ("dump-losses", "Dump per-symbol losses into a file named losses")
    ("prefix",po::value<std::string>()->default_value(""), "Prefix for output files")
    ;

  hidden.add_options()
    ("input-file", po::value<string>(), "input file")
    ;

  po::options_description cmdline_options;
  cmdline_options.add(generic).add(hidden).add(dumping);

  po::options_description visible_options;
  visible_options.add(generic).add(dumping);


  po::positional_options_description p;
  p.add("input-file", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
      options(cmdline_options).positional(p).run(), vm);
  po::notify(vm);

  string usage = "Usage: score_file [OPTIONS]... FILENAME\n";

  if (vm.count("help") || !vm.count("input-file")) {
    cout << usage << visible_options << "\n";
    exit(1);
  }

  num_types = vm["num-types"].as<int>();

  init_rng();

  if (vm.count("input-file")) {
    double score;
    score = score_file(vm);
    cout << "log-loss: " << score << endl;
    //exit(0); // force exit to avoid expensive cleanup
  }
  free_rng();
}
