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

#include "libplump/hpyp_model.h"

#include <cmath>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>

#include "libplump/utils.h"
#include "libplump/subseq.h"


namespace gatsby { namespace libplump {

HPYPModel::HPYPModel(seq_type& seq,
                     INodeManager& nodeManager,
                     const IAddRemoveRestaurant& restaurant,
                     IParameters& parameters,
                     int numTypes) 
    : seq(seq), 
      contextTree_(new ContextTree(nodeManager, seq)), 
      contextTree(*contextTree_), 
      restaurant(restaurant),
      parameters(parameters), 
      numTypes(numTypes) {
    baseProb = 1./((double) numTypes);
  }
 

void HPYPModel::insertRoot(e_type obs) {
  WrappedNodeList root_path = this->contextTree.findLongestSuffix(0,0);
  d_vec discount_path = this->parameters.getDiscounts(root_path);
  d_vec concentration_path = this->parameters.getConcentrations(root_path,
                                                                discount_path);
  d_vec prob_path = this->computeProbabilityPath(root_path,
                                                 discount_path,
                                                 concentration_path,
                                                 obs);
  this->updatePath(root_path,
                   prob_path,
                   discount_path,
                   concentration_path,
                   obs);
}


/**
 * Insert a context into the tree and handle a potentially occuring 
 * split.
 */
WrappedNodeList HPYPModel::insertContext(l_type start, l_type stop) {
  typedef ContextTree::InsertionResult InsertionResult;

  // insert context
  InsertionResult insertionResult = contextTree.insert(start, stop);

  // handle split if one occurred
  if (insertionResult.action != InsertionResult::INSERT_ACTION_NO_SPLIT) { 
    WrappedNode nodeA, nodeC;
    WrappedNodeList::iterator i = insertionResult.path.end();

    // move iterator to point to the split node which is either the 
    // last or second to last element on the path, depending on 
    // whether the inserted node was a suffix
    switch (insertionResult.action) {
      case InsertionResult::INSERT_ACTION_SPLIT :
        // split where inserted node is not a suffix of a node in 
        // the tree
        i--; i--; // split node is above the newly inserted node
        break;
      case InsertionResult::INSERT_ACTION_SPLIT_SUFFIX :
        // split where inserted node is a suffix -> inserted node
        // is the shorter node
        i--; // split node is the newly inserted node
        break;
      case InsertionResult::INSERT_ACTION_NO_SPLIT :
        break;
    }
    nodeC = *i;
    i--; // move iterator to parent
    nodeA = *i;
    this->handleSplit(nodeA, 
                      insertionResult.splitChild,
                      nodeC); 
  }
  return insertionResult.path;    
}


d_vec HPYPModel::insertContextAndObservation(l_type start, 
                                             l_type stop,
                                             e_type obs) {
  // insert context (and handle a potential split)
  WrappedNodeList path = insertContext(start, stop);

  // insert observation
  d_vec discountPath = this->parameters.getDiscounts(path);
  d_vec concentrationPath = this->parameters.getConcentrations(path,
                                                               discountPath);
  d_vec probabilityPath = this->computeProbabilityPath(path,
                                                       discountPath,
                                                       concentrationPath,
                                                       obs);
  this->updatePath(path, probabilityPath, discountPath,
                   concentrationPath, obs);
  return probabilityPath; 
}


d_vec HPYPModel::insertObservation(l_type start, l_type stop, e_type obs) {
  tracer << "HPYPModel::insertObservation(" << start << ", " << stop 
         << ", " << obs << ")" << std::endl;
  
  WrappedNodeList path = this->contextTree.findLongestSuffix(start,stop);
  
  tracer << "  HPYPModel::insertObservation: longest suffix path: " 
         << std::endl << this->contextTree.pathToString(path) << std::endl;
  
  d_vec discount_path = this->parameters.getDiscounts(path);
  d_vec concentration_path = this->parameters.getConcentrations(path,
                                                                discount_path);
  d_vec prob_path = this->computeProbabilityPath(path,
                                                 discount_path,
                                                 concentration_path,
                                                 obs);
  this->updatePath(path, prob_path, discount_path,
                   concentration_path, obs);
  return prob_path; 
}


void HPYPModel::removeObservation(
    l_type start, l_type stop, e_type obs, 
    const HPYPModel::PayloadDataPath& payloadDataPath) {
  tracer << "HPYPModel::removeObservation(" << start << ", " << stop 
         << ", " << obs << ")" << std::endl;

  WrappedNodeList path = this->contextTree.findLongestSuffix(start, stop);
  
  tracer << "  HPYPModel::removeObservation: longest suffix path: " 
         << std::endl << contextTree.pathToString(path) << std::endl;
  
  d_vec discountPath = this->parameters.getDiscounts(path);
  
  this->removeObservationFromPath(path, 
                                  discountPath,
                                  obs,
                                  payloadDataPath);
}


d_vec HPYPModel::computeLosses(l_type start, l_type stop) {
  d_vec losses;

  // deal with first symbol: add loss and insert customer
  losses.push_back(log2((double) this->numTypes));
  insertRoot(this->seq[start]);

  // start timer
  clock_t start_t,end_t;
  start_t = clock();

  for (l_type i=start+1; i < stop; i++) {
    d_vec prob_path = this->insertContextAndObservation(start,i,this->seq[i]);
    double prob = prob_path[prob_path.size()-2];
    losses.push_back(-log2(prob));
    
    if (i%10000==0) {
      end_t = clock();
      std::cerr << makeProgressBarString(i/(double)stop) << " " 
                << ((double)i*CLOCKS_PER_SEC)/(end_t-start_t) << " chars/sec" 
                <<  "\r";
    }

  }
  end_t = clock();

  std::cerr << makeProgressBarString(1) << " " 
            << ((double)stop*CLOCKS_PER_SEC)/(end_t-start_t) << " chars/sec" 
            <<  std::endl;
  
  return losses;
}


d_vec HPYPModel::predictSequence(l_type start, l_type stop, PredictMode mode) {
  d_vec probs;
  for (l_type i = start; i < stop; i++) {
    switch(mode) {
      case ABOVE:
        probs.push_back(this->predict(start, i, this->seq[i]));
        break;
      case FRAGMENT:
        probs.push_back(this->predictWithFragmentation(start, i, this->seq[i]));
        break;
      case BELOW:
        probs.push_back(this->predictBelow(start, i, this->seq[i]));
        break;
    }
  }
  return probs;
}


void HPYPModel::buildTree(l_type stop) {
  this->insertRoot(this->seq[0]);
  for (l_type i=1; i < stop; ++i) {
    this->insertContextAndObservation(0, i, this->seq[i]);
  }
}


void HPYPModel::updateTree(l_type start, l_type stop) {
  for (l_type i = start; i < stop; ++i) {
    this->insertContextAndObservation(0, i, this->seq[i]);
  }
}


double HPYPModel::predict(l_type start, l_type stop, e_type obs) {
  WrappedNodeList path = this->contextTree.findLongestSuffix(start,stop);
  d_vec discount_path = this->parameters.getDiscounts(path);
  d_vec concentration_path = this->parameters.getConcentrations(path,
                                                                discount_path);

  d_vec prob_path = this->computeProbabilityPath(path,
                                                 discount_path,
                                                 concentration_path,
                                                 obs);
  return prob_path.back();
}


/**
 * Predict prob; in the case of required fragmentation, predict form 
 * _below_ the split point!
 */
double HPYPModel::predictBelow(l_type start, l_type stop, e_type obs) {
  WrappedNodeList path = 
      this->contextTree.findLongestSuffixVirtual(start,stop).second;
  d_vec discount_path = this->parameters.getDiscounts(path);
  d_vec concentration_path = this->parameters.getConcentrations(path,
                                                                discount_path);

  d_vec prob_path = this->computeProbabilityPath(path,
                                                 discount_path,
                                                 concentration_path,
                                                 obs);
  return prob_path.back();
}


double HPYPModel::predictWithFragmentation(l_type start, 
                                           l_type stop,
                                           e_type obs) {
  std::pair<int, WrappedNodeList> path = contextTree.findLongestSuffixVirtual(
      start, stop);

  d_vec discountPath = parameters.getDiscounts(path.second);
  d_vec concentrationPath = parameters.getConcentrations(path.second,
                                                         discountPath);
  d_vec probabilityPath = this->computeProbabilityPath(path.second,
                                                       discountPath,
                                                       concentrationPath,
                                                       obs);

  double probability = 0;
  if (path.first != 0) {
    // create a new payload for the node we are predicting from
    // fragmentation -- last probability on the path needs to be recomputed
    // by creating a new, split node of length path.second. 
    void* splitNode = this->restaurant.getFactory().make();
    WrappedNodeList::iterator it = path.second.end();
    it--; it--; // one before last; parent of node we need to split
    int parentLength = it->end - it->start; 
    // path.first is length of parent after split
    double discountAfter = parameters.getDiscount(path.first, it->end - it->start);
    it++; // last node
    this->restaurant.updateAfterSplit(
        it->payload,
        splitNode,
        discountPath.back(),
        discountAfter,
        true); // update splitNode only
    double discountFragmented = this->parameters.getDiscount(parentLength,
                                                             path.first);
    double concentrationFragmented = this->parameters.getConcentration(
        discountFragmented, parentLength, path.first);
    probability = this->restaurant.computeProbability(
        splitNode, obs, probabilityPath[probabilityPath.size()-2],
        discountFragmented, concentrationFragmented);
    this->restaurant.getFactory().recycle(splitNode);
  } else {
    probability = probabilityPath.back();
  }
  return probability;
}


d_vec HPYPModel::predictiveDistribution(l_type start, l_type stop) {
  d_vec predictive;
  predictive.reserve(this->numTypes);
  WrappedNodeList path = this->contextTree.findLongestSuffix(start,stop);
  d_vec discount_path = this->parameters.getDiscounts(path);
  d_vec concentration_path = this->parameters.getConcentrations(path,
                                                                discount_path);

  for(int i = 0; i < this->numTypes; ++i) {
    d_vec prob_path = this->computeProbabilityPath(path,
                                                   discount_path,
                                                   concentration_path,
                                                   i);
    predictive.push_back(prob_path.back());
  }

  return predictive;
}


d_vec HPYPModel::predictiveDistributionWithMixing(l_type start, 
                                                  l_type stop, 
                                                  d_vec& mixingWeights) {
  d_vec predictive;
  predictive.reserve(this->numTypes);
  const WrappedNodeList path = this->contextTree.findLongestSuffix(start, stop);
  const d_vec discount_path = this->parameters.getDiscounts(path);
  const d_vec concentration_path 
      = this->parameters.getConcentrations(path, discount_path);

  for(int i = 0; i < this->numTypes; ++i) {
    d_vec prob_path = this->computeProbabilityPath(path,
                                                   discount_path,
                                                   concentration_path,
                                                   i);
    double prob = 0;
    double sum = 0;
    for (int j = 0; 
         j < std::min<int>(mixingWeights.size(), prob_path.size());
         ++j) {
      prob += mixingWeights[j]*prob_path[j];
      sum += mixingWeights[j];
    }
    predictive.push_back(prob + (1 - sum)*prob_path.back());
  }

  return predictive;
}
        

/**
 * Given a path from the root to some node (not a leaf), 
 * perform Gibbs sampling of the last node by repeatedly 
 * removing and adding customers, cus times for each type s.
 */
void HPYPModel::gibbsSamplePath(
    const WrappedNodeList& path, 
    const d_vec& discountPath, 
    const d_vec& concentrationPath, 
    const HPYPModel::PayloadDataPath& payloadDataPath,
    double baseProb) {
  assert(path.size() > 0);
  assert(path.size() == discountPath.size());
  assert(path.size() == concentrationPath.size());

  bool useAdditionalData = payloadDataPath.size() == path.size();
  void* main = path.back().payload;
  const IAddRemoveRestaurant& r = this->restaurant; // shortcut
  
  IHPYPBaseRestaurant::TypeVector types = r.getTypeVector(main);
  for(IHPYPBaseRestaurant::TypeVectorIterator it = types.begin(); 
      it != types.end(); ++it) { // for each type of customer 

    e_type type = *it;
    l_type cw = r.getC(main, type);
    
    if (cw == 1) {
      continue; // no point in reseating in a 1 customer restaurant
    }

    d_vec probabilityPath = this->computeProbabilityPath(path,
                                                         discountPath,
                                                         concentrationPath,
                                                         type);
    WrappedNodeList::const_iterator current;
    for (l_type i = 0; i < cw; ++i) { // for each customer of this type
      current = --path.end(); // set current to last restaurant in path
      
      // index into d/alpha vectors for current restaurant
      int j = discountPath.size() - 1;
      
      bool goUp = true;
      while(goUp && j != -1) {
        void* additionalData = NULL;
        if (useAdditionalData) {
          additionalData = payloadDataPath[j].get();
        }
        bool removed = r.removeCustomer(current->payload,
                                        type,
                                        discountPath[j],
                                        additionalData);
        if (removed) {
          --current;
          --j;
        } else {
          goUp = false;
        }
      }

      // recompute probabilities back down; need not recompute last probability
      while(j < (int)probabilityPath.size() - 1) {
        if (j == -1) {
          // can't recompute base distribution probabilities at prob_path[0]
          j = 0; // 
          ++current;
        }
        probabilityPath[j+1] = r.computeProbability(current->payload, 
                                                    type,
                                                    probabilityPath[j],
                                                    discountPath[j],
                                                    concentrationPath[j]);
        ++j; 
        ++current;
      }

      // set current and j to the last restaurant on path
      current = --path.end(); 
      j = discountPath.size() - 1; 
      goUp = true;
      while(goUp && j != -1) {
        void* additionalData = NULL;
        if (useAdditionalData) {
          additionalData = payloadDataPath[j].get();
        }
        bool inserted = r.addCustomer(current->payload, 
                                      type,
                                      probabilityPath[j],
                                      discountPath[j],
                                      concentrationPath[j], 
                                      additionalData);
        if (inserted) {
          current--;
          j--;
        } else {
          goUp = false;
        }
      }

    }

  }

}


boost::shared_ptr<void> HPYPModel::makeAdditionalDataPtr(void* payload, 
                                                         double discount, 
                                                         double concentration) {
  return boost::shared_ptr<void>(
      this->restaurant.createAdditionalData(payload,
        discount,
        concentration),
      boost::bind(&IAddRemoveRestaurant::freeAdditionalData, 
        &(this->restaurant),
        _1));
}


void HPYPModel::runGibbsSampler() {
  ContextTree::DFSPathIterator pathIterator = contextTree.getDFSPathIterator();
  d_vec discountPath = parameters.getDiscounts(*pathIterator);
  d_vec concentrationPath = parameters.getConcentrations(*pathIterator, 
                                                         discountPath);

  // initialize payloadDataPath; by using shared_ptr with the proper
  // destruction function, all clean-up should be automatic.
  HPYPModel::PayloadDataPath payloadDataPath;
  int j = 0;
  for (WrappedNodeList::const_iterator it = (*pathIterator).begin();
       it != (*pathIterator).end(); ++it) {
      payloadDataPath.push_back(this->makeAdditionalDataPtr(
            it->payload, discountPath[j], concentrationPath[j]));
    j++;
  }

  this->gibbsSamplePath(*pathIterator, discountPath, concentrationPath,
                        payloadDataPath, baseProb);

  size_t pathLength = (*pathIterator).size();
  
  while(pathIterator.hasMore()) { // loop over all paths in the tree
    ++pathIterator;
    if ((*pathIterator).size() == 0) {
      break;
    }

    if ((*pathIterator).size() == pathLength) {
      // sibling
      discountPath.pop_back();
      parameters.extendDiscounts(*pathIterator, discountPath);
      concentrationPath.pop_back();
      parameters.extendConcentrations(*pathIterator,
                                      discountPath, 
                                      concentrationPath);
      payloadDataPath.pop_back();
      payloadDataPath.push_back(
          this->makeAdditionalDataPtr((*pathIterator).back().payload, 
                                      discountPath.back(), 
                                      concentrationPath.back()));

    } else {
      if ((*pathIterator).size() == pathLength - 1) {
        // we went up -- just drop the last term
        discountPath.pop_back();
        concentrationPath.pop_back();
        payloadDataPath.pop_back();
      } else {
        // we went up one and then down some number of levels -- recompute
        discountPath.pop_back();
        concentrationPath.pop_back();
        parameters.extendDiscounts(*pathIterator, discountPath);
        parameters.extendConcentrations(*pathIterator,
                                        discountPath,
                                        concentrationPath);
        payloadDataPath.pop_back();
        WrappedNodeList::const_iterator it = (*pathIterator).begin();
        // move it to the first item not covered by the payloadDataPath
        for (size_t i = 0; i < payloadDataPath.size(); ++i) {
          ++it;
        }
        for (size_t i = payloadDataPath.size(); i < discountPath.size(); ++i) {
          payloadDataPath.push_back(
              this->makeAdditionalDataPtr(it->payload, 
                                          discountPath[i], 
                                          concentrationPath[i]));

          ++it;
        }
        assert(it == (*pathIterator).end());
      }
    }
    pathLength = (*pathIterator).size();
    
    tracer << (*pathIterator).size() << " " << discountPath.size() << " " 
           << concentrationPath.size() << " " <<  payloadDataPath.size() 
           << std::endl;
    
    
    this->gibbsSamplePath(*pathIterator, discountPath, concentrationPath,
                          payloadDataPath, baseProb);
  }
}


d_vec HPYPModel::computeProbabilityPath(const WrappedNodeList& path, 
                                        const d_vec& discount_path, 
                                        const d_vec& concentration_path,
                                        e_type obs) {
  d_vec out;
  out.reserve(path.size() + 1);
  
  double prob = this->baseProb; // base distribution
  out.push_back(prob);
  
  int j = 0;
  for(WrappedNodeList::const_iterator it = path.begin();
      it!= path.end();
      ++it) {
    prob = this->restaurant.computeProbability(it->payload, 
                                               obs,
                                               prob,
                                               discount_path[j],
                                               concentration_path[j]);
    out.push_back(prob);
    j++;
  } 
  return out;
}


void HPYPModel::updatePath(const WrappedNodeList& path, 
                           const d_vec& prob_path, 
                           const d_vec& discount_path, 
                           const d_vec& concentration_path, 
                           e_type obs) {

  unsigned int j=path.size()-1;
  for(WrappedNodeList::const_reverse_iterator it = path.rbegin(); 
      it != path.rend();
      ++it) {
    bool newTable = this->restaurant.addCustomer(it->payload,
                                                 obs,
                                                 prob_path[j],
                                                 discount_path[j],
                                                 concentration_path[j]);
    if (!newTable) {
      break;
    }
    j--;
  }
}
    

void HPYPModel::removeObservationFromPath(
    const WrappedNodeList& path,
    const d_vec& discountPath,
    e_type obs,
    const HPYPModel::PayloadDataPath& payloadDataPath) {
  int j = path.size()-1;

  for(WrappedNodeList::const_reverse_iterator it = path.rbegin();
    it != path.rend(); it++) {

    void* payloadData = NULL;
    if (payloadDataPath.size() == path.size()) {
      payloadData = payloadDataPath[j].get();
    }

    bool tableDeleted = this->restaurant.removeCustomer(it->payload,
                                                        obs,
                                                        discountPath[j],
                                                        payloadData);
    if (!tableDeleted)
      break;
    j--;
  }
}
        

void HPYPModel::handleSplit(const WrappedNode& nodeA,
                            const WrappedNode& nodeB, 
                            const WrappedNode& nodeC) {
    int lengthA = nodeA.end - nodeA.start;
    int lengthB = nodeB.end - nodeB.start;
    int lengthC = nodeC.end - nodeC.start;
    
    // parent context should be shorter than both its children 
    assert(lengthA < lengthB && lengthA < lengthC);
    // length of node that required splitting should be longer than the result
    assert(lengthC < lengthB);

    double discBBeforeSplit = this->parameters.getDiscount(lengthA, lengthB);
    double discBAfterSplit  = this->parameters.getDiscount(lengthC, lengthB);
    this->restaurant.updateAfterSplit(nodeB.payload, 
                                      nodeC.payload, 
                                      discBBeforeSplit,
                                      discBAfterSplit);
}


bool HPYPModel::checkConsistency(const WrappedNode& node, 
                      const std::list<WrappedNode>& children) const {
  bool consistent = this->restaurant.checkConsistency(node.payload);
  std::map<e_type, int> table_counts;
  for(std::list<WrappedNode>::const_iterator it = children.begin();
      it != children.end(); ++it) {
    IHPYPBaseRestaurant::TypeVector keys = 
        this->restaurant.getTypeVector(it->payload);
    for(IHPYPBaseRestaurant::TypeVectorIterator key_it = keys.begin();
        key_it != keys.end(); ++key_it) {
      table_counts[*key_it] += this->restaurant.getT(it->payload, *key_it);
    }
  }

  for(std::map<e_type, int>::iterator it = table_counts.begin();
      it != table_counts.end(); ++it) {
    consistent = (this->restaurant.getC(node.payload, (*it).first) 
                  >= (*it).second) && consistent;
    if (!consistent) {
      std::cerr << "Child table sum is: " << (*it).second 
                << ", parent customers is: " 
                << this->restaurant.getC(node.payload, (*it).first) 
                << std::endl;
    }
  }
  return consistent;
}


bool HPYPModel::checkConsistency() const {
  CheckConsistencyVisitor v(*this);
  this->contextTree.visitDFSWithChildren(v);
  return v.consistent;
}


HPYPModel::ToStringVisitor::ToStringVisitor(seq_type& seq, 
    const IHPYPBaseRestaurant& restaurant) 
    : outstream(), seq(seq), restaurant(restaurant) {}


void HPYPModel::ToStringVisitor::operator()(const WrappedNode& n) {
  for (l_type i = 0; i < n.depth; ++i) {
    this->outstream << " "; 
  }
  this->outstream << SubSeq::toString(n.start, n.end, seq);
  this->outstream << " " << this->restaurant.toString(n.payload);
  this->outstream << std::endl;
}


HPYPModel::CheckConsistencyVisitor::CheckConsistencyVisitor(
    const HPYPModel& model) : consistent(true), model(model) {}


void HPYPModel::CheckConsistencyVisitor::operator()(
    WrappedNode& n, std::list<WrappedNode>& children) {
  bool nodeConsistent = this->model.checkConsistency(n, children);
  if (!nodeConsistent) {
    std::cerr << "Node " << n.toString() << " not consistent!" << std::endl;
  }
  consistent = consistent && nodeConsistent;
}


std::string HPYPModel::toString() {
  HPYPModel::ToStringVisitor visitor(this->seq, this->restaurant);
  this->contextTree.visitDFS(visitor);
  return visitor.outstream.str();
}


}} // namespace gatsby::libplump
