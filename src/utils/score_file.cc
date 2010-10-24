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
//#include <boost/filesystem.hpp>
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
  return new SimpleParameters(vm["disc"].as<d_vec>(), vm["alpha"].as<double>());
}

IAddRemoveRestaurant* getRestaurant(po::variables_map& vm) {
  switch(vm["restaurant"].as<int>()) {
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

void runSampler(po::variables_map& vm, HPYPModel& model) {
  model.runGibbsSampler();
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

  boost::scoped_ptr<IParameters> parameters(getParameters(vm));
  boost::scoped_ptr<IAddRemoveRestaurant> restaurant(getRestaurant(vm));
  boost::scoped_ptr<INodeManager> nodeManager(
      getNodeManager(vm, restaurant->getFactory()));

  HPYPModel model(seq, *nodeManager, *restaurant, *parameters, num_types);

  //cout << "Discounts: " << iterableToString(m.parameters.discounts) << endl;
  //cout << "Concentration: " << m.parameters.alpha << endl;
  d_vec losses;
  if (vm.count("load-serialized-nodes")) {
    Serializer nodeSerializer(vm["load-serialized-nodes"].as<string>());
    nodeSerializer.loadNodesAndPayloads(*nodeManager, restaurant->getFactory());
  } else {
    losses = model.computeLosses(0,seq.size());
  }
  
  if (vm.count("debug")) {
    model.checkConsistency();
  }

  cout << "Training loss: " << mean(losses) << endl;

  if (vm.count("print-tree")) {
    cout << model.toString() << endl;
  }

  //if (vm.count("traverse-tree")) {
  //   typename Tree::DFSPathIterator it = contextTree.dfs_path_iterator();
  //    while(it.hasMore()) {
  //        cout << contextTree.pathToString(*it) << endl;
  //        it++;
  //    }
  //}

  if (vm.count("dump-losses")) {
    iterableToCSVFile(losses,"losses");
  }

  //if (vm.count("dump-stats")) {
  //    typename Model::CollectStatisticsVisitor csv;
  //    contextTree.visitDFS(csv);
  //    cout << csv.toString() << endl;
  //}



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
        runSampler(vm, model);
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
        //if (vm.count("per-sample-dump")) {
        //  if (vm.count("dump-stats")) {
        //    typename Model::CollectStatisticsVisitor csv;
        //    contextTree.visitDFS(csv);
        //    cout << csv.toString() << endl;
        //    cout << "Stirling stats: " << stirling_generator::statsToString() << endl;
        //  }
        if (vm.count("debug") && vm.count("print-tree")) {
          cout << model.toString() << endl;
        }
      }
    }
  }
  //         cout << "Using test mode: ";
  //         switch (mode) {
  //             case 1:
  //                 cout << "particle filter" << endl;
  //                 losses = m.computeLosses(start_pos,seq.size());
  //                 break;
  //             case 2:
  //                 cout << "Prediction (fast mode)" << endl;
  //                 losses = predict(vm,m, start_pos,seq);
  //                 log2_vec(losses);
  //                 mult_vec(losses,-1.);
  //                 break;
  // 
  //             case 3:
  //                 cout << "Prediction without fragmentation (safe mode)" << endl;
  //                 losses.clear();
  //                 for(unsigned int i=start_pos+1;i<seq.size();++i) {    
  //                     cout << i << endl;
  //                     d_vec predictive = m.predictiveDistribution(start_pos,i);
  //                     if (!closeTo(sum(predictive),1)) {
  //                         cerr << "Warning: Predictive distribution does not sum to 1: " << sum(predictive) << endl;
  //                     }
  //                     //cout << sum(predictive) << endl;
  //                     losses.push_back(-log2(predictive[seq[i]]));
  //                 }
  //                 break;
  //             case 4:
  //                 cout << "Prediction with fragmentation" << endl;
  //                 break;
  //             case 5:
  //                 cout << "Prediction without fragmentation with Gibbs sampling" << endl;
  //                 losses.clear();
  //                 losses.assign(seq.size()-start_pos-1,0);
  //                 int num_samples = vm["samples"].as<int>();
  //                 for (int s=0;s<num_samples;s++) {
  //                     cout << "Sampling iteration: " << s +1 << endl;
  //                     d_vec cur_losses = predict(vm,m,start_pos,seq);
  //                     current_sample_losses.push_back(prob2loss(cur_losses));
  //                     losses_f << current_sample_losses.back() << ", ";
  //                     losses_f.flush();
  //                     for (unsigned int i=0; i<losses.size();++i) {
  //                         losses[i] += cur_losses[i];
  //                     }
  //                     d_vec tmp = losses;
  //                     mult_vec(tmp,1/(double)(s+1));
  //                     log2_vec(tmp);
  //                     mult_vec(tmp,-1.);
  //                     cout << "Averaged loss after " << s+1 << " samples: " << mean(tmp) << ", Perplexity: " << pow(2,mean(tmp)) << endl;
  //                     runSampler(vm,m);
  //                 }
  //                 mult_vec(losses,1/(double)(num_samples));
  //                 log2_vec(losses);
  //                 mult_vec(losses,-1.);
  //                 break;
  //         } 
  //         losses_f << mean(losses);
  //     }

  if (vm.count("save-serialized-nodes")) {
    Serializer nodeSerializer(vm["save-serialized-nodes"].as<string>());
    nodeSerializer.saveNodesAndPayloads(*nodeManager, restaurant->getFactory());
  }

  return mean(losses);
}

//template<typename NodeManager, typename stirling_generator>
//void incrementally_build(po::variables_map& vm) {
//    typedef typename NodeManager::Payload Payload; 
//    typedef ContextTree<NodeManager> Tree;
//    typedef ContextTreeModel<Tree, HPYPPayloadManager<Payload> > Model;
//	
//    string filename = vm["input-file"].as<string>();
//    
//    fs::path input_path(filename);
//    if (!fs::exists(input_path)) {
//        cerr << "Input file not found!" << endl;
//        exit(1);
//    }
//    
//    seq_type seq;
//    pushFileToVec<e_type>(filename, seq, vm["head"].as<int>());
//    cout << "Sequence length: " << seq.size() << endl;
//    ostringstream out_fn;
//    out_fn << "stats_" << fs::basename(input_path) << ".csv";
//    ofstream stats_f(out_fn.str().c_str());
//    std::vector<double> losses;
//    NodeManager nodeManager;
//    Tree contextTree(nodeManager, seq);
//    Model m(contextTree, num_types);
//    m.parameters.discounts = vm["disc"].as<d_vec>();
//    m.parameters.alpha = vm["alpha"].as<double>();
//    int interval = vm["interval"].as<int>();
//    m.buildTree(1);
//    typename Model::CollectStatisticsVisitor csv;
//    stats_f << csv.header() << endl;
//    d_vec times;
//    for (unsigned int i=1;i<seq.size();i++) {
//        if (i%interval==0) {
//            cout << i << endl;
//            if (vm.count("sample")) {
//                for (int j=0;j<vm["sample"].as<int>();j++) {
//                    tic();
//                    runSampler(vm,m);
//                    times.push_back(toc());
//                    cout << "Time per sample: " << times.back() << endl;
//                    cout << j << endl;
//                }
//            }
//            if (vm.count("dump-stats")) {
//                contextTree.visitDFS(csv);
//                stats_f <<  i << ", " << iterableToString(csv.toVec()) << endl;
//            }
//        }
//        m.updateTree(i,i+1);
//    }
//    ofstream time_f("times.csv");
//    time_f << iterableToString(times);
//}

int main(int argc, char* argv[]) {
  //const double discs[] =  {0.05, 0.7, 0.8, 0.82, 0.84, 0.88, 0.91, 0.92, 0.93, 0.94, 0.95, 0.9999};
  //const double discs[] =  {0.05, 0.7, 0.8, 0.82, 0.84, 0.88, 0.91, 0.92, 0.93, 0.94, 0.95};
  const double sm_disc[] = {.62, .69, .74, .80, .95};
  //const double discs[] =  {0.05, 0.7, 0.8, 0.82, 0.84, 0.88, 0.91, 0.92, 0.93, 0.94, 0.95, 0.9999};
  d_vec default_discounts;
  //default_discounts.assign(discs,&discs[11]);
  default_discounts.assign(sm_disc,&sm_disc[5]);
  // Declare the supported options.
  po::options_description generic("Generic options");
  po::options_description dumping("Dumping options");
  po::options_description hidden("Hidden options");
  generic.add_options()
    ("help", "produce help message")
    ("debug,D", "Print debugging output")
    ("hash-table-tree,H", "Use hash table to store tree")
    ("print-tree", "Print the context tree to the screen")
    ("traverse-tree", "Traverse the tree in the order the sampling is performed to the screen")
    ("incremental", "Build tree incrementally")
    ("switching", "Use SwitchingNodeManager")
    ("fragment", po::value<int>()->default_value(1), "1: nofrag; 2: frag; 3:below")
    ("read-int32", "Read input data as 32 bit integers")
    ("test-file", po::value<string>(), "Test file")
    ("save-serialized-nodes", po::value<string>(), "File to contain serialized nodes")
    ("load-serialized-nodes", po::value<string>(), "File to contain serialized nodes")
    ("head",po::value<int>()->default_value(0), "If given, cuts input to this number of symbols")
    ("mode", po::value<int>()->default_value(1), "1: particle filter, 2: no fragment, 3: fragment")
    ("restaurant", po::value<int>()->default_value(1),
     "1: SimpleFull, 2: Histogram, 3: ReinstantiatingCompact, 4: StirlingCompact, 5: Switching")
    ("sample",po::value<int>()->default_value(0), "Number of Gibbs iteration (incremental mode)")
    ("burn-in",po::value<int>()->default_value(0), "Number of Gibbs iterations for burn in")
    ("samples,s",po::value<int>()->default_value(1), "Number of samples used for prediction")
    ("interval",po::value<int>()->default_value(10000), "Interval between sampling/stats runs")
    ("num-types", po::value<int>()->default_value(256), "Number of types") 
    ("alpha,a", po::value<double>()->default_value(5), "Concentration parameter") 
    ("disc,d", po::value<d_vec>()->default_value(default_discounts,"..."), "Discount parameter(s)") 
    ;

  dumping.add_options()
    ("dump-losses", "Dump per-symbol losses into a file named losses")
    ("dump-stats", "Dump tree statistics")
    ("per-sample-dump", "Dump statistics for each sample")
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
  //cerr << "Stirling stats: " << stirling_generator::statsToString() << endl;
  free_rng();
}
