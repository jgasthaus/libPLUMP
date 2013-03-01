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

#ifndef HPYP_MODEL_H_
#define HPYP_MODEL_H_


#include <vector>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "libplump/config.h"
#include "libplump/context_tree.h"
#include "libplump/node_manager_interface.h"
#include "libplump/hpyp_restaurant_interface.h"
#include "libplump/hpyp_parameters_interface.h"

namespace gatsby { namespace libplump {
  

class HPYPModel {
  public: 
    typedef std::vector<boost::shared_ptr<void> > PayloadDataPath;

    enum PredictMode {ABOVE, FRAGMENT, BELOW};
    
    /**
     * Construct a new HPYP model using the given nodeManager and restaurant.
     */
    HPYPModel(seq_type& seq,
              INodeManager& nodeManager, 
              const IAddRemoveRestaurant& restaurant,
              IParameters& parameters,
              int numTypes);

    ~HPYPModel() {}

    /**
     * Create the root node and insert the given observation into it.
     */
    void insertRoot(e_type obs);

    /**
     * Insert a context into the tree and handle a potentially occuring 
     * split.
     */
    WrappedNodeList insertContext(l_type start, l_type stop);

    /**
     * Insert a context into the tree, compute the probability of the
     * observation under the current state, and insert the observation
     * into the model.
     *
     * A basic particle filter based on HPYPModel thus consists of 
     * repeated calls to this method. 
     */
    d_vec insertContextAndObservation(l_type start, l_type stop, e_type obs);

    /**
     * Insert an observation into an existing context.
     */
    d_vec insertObservation(l_type start, l_type stop, e_type obs, 
                            WrappedNodeList* path = NULL);

    /**
     * Remove an observation of type obs from the given context.
     *
     * This assumes that the given context (start,stop) exists explicitly
     * in the tree. 
     */
    void removeObservation(l_type start,
                           l_type stop,
                           e_type obs,
                           const PayloadDataPath&, 
                           WrappedNodeList* path = NULL);

    /**
     * For each i in [start, stop), insert the context 
     * (start, i) into the tree, compute the probability of seq[i] and 
     * then insert seq[i] into the model. 
     */
    d_vec computeLosses(l_type start, l_type stop);
    
    
    d_vec computeLossesWithDeletion(l_type start, l_type stop, l_type lag);

    
    /**
     * Update model by calling removeCustomer followed by
     * addCustomer on all observations/
     */
    void removeAddSweep(l_type start, l_type stop);


    /**
     * Build the context tree for all contexts [start,i) for 
     * i in [start, stop].
     */
    void buildTree(l_type stop);

    /**
     * Update the context tree by inserting all contexts [start, i]
     * for i in [start, stop].
     */
    void updateTree(l_type start, l_type stop);

    /**
     * Compute the predictive probability of the given observation
     * in the context [start, stop). If the required context is not in tree, 
     * use the parent node of the required context for prediction (i.e. this
     * method does _not_ use fragmentation for prediction)..
     */
    double predict(l_type start, l_type stop, e_type obs);

    /**
     * Predict prob; in the case of required fragmentation, predict form 
     * _below_ the split point!
     */
    double predictBelow(l_type start, l_type stop, e_type obs);

    /**
     * Compute the predictive probability of the given observation
     * in the context [start, stop). If the required context is not in tree, 
     * produce a temporary restaurant using fragmentation and use that for
     * prediction.
     */
    double predictWithFragmentation(l_type start, l_type stop, e_type obs);


    /**
     * Compute predictive probability for a sequence of observations, 
     * i.e. for each i in [start, stop), compute p(x_i|x_{start:i}).
     */
    d_vec predictSequence(l_type start, l_type stop,
                          PredictMode mode = ABOVE);

    /**
     * Compute the entire predictive distribution in the given context.
     */
    d_vec predictiveDistribution(l_type start, l_type stop);

    /**
     * Compute the entire predictive distribution in the given context by
     * mixing the distributions in all contexts up to the root with
     * the given weights.
     */
    d_vec predictiveDistributionWithMixing(l_type start, 
                                           l_type stop,
                                           d_vec& mixingWeights);

    /**
     * Run one iteration of Gibbs sampling in the model.
     */
    void runGibbsSampler();

    /**
     * Return a string representation of entire model.
     */
    std::string toString();

    /**
     * Check the consistency of the current seating arrangement by checking
     * that for every (parent, children, type) triple the sum of tables of that
     * type in the child restaurants matches the number of customers in the
     * parent restaurant.
     */
    bool checkConsistency() const;


  private:

    /** 
     * Compute the predictive probability for a symbol along a path, starting
     * with the base distribution at the root. At position i it contains the
     * parentProb for node i.
     */
    d_vec computeProbabilityPath(const WrappedNodeList& path, 
                                 const d_vec& discount_path, 
                                 const d_vec& concentration_path,
                                 e_type obs);


    /**
     * Insert a customer of type obs into the last node in path, then
     * recursively insert customers up the path if a new table was created by
     * the last insertion.
     */
    void updatePath(const WrappedNodeList& path, 
                    const d_vec& prob_path, 
                    const d_vec& discount_path, 
                    const d_vec& concentration_path, 
                    e_type obs);
    


    /**
     * Remove and observation from the last node of the path, then recursively
     * remove customers up the path if a table was deleted by the last 
     * removal.
     *
     * @param payloadCache must be either NULL or an array of void*
     *        of the same length as path. These void* are then passed
     *        to the restaurant when removeCustomer is called. This can be used
     *        to e.g. cache stirling numbers along a path. 
     *
     */
    void removeObservationFromPath(const WrappedNodeList& path, 
                                   const d_vec& discountPath, 
                                   e_type obs,
                                   const PayloadDataPath& payloadDataPath);
    

    /**
     * Handle the split of a node that occurred while inserting a new node.
     *
     * Assume we are inserting a node X, and the longest common suffix
     * of X with a node in the tree lies within a node B, whose parent
     * in the tree is a node A. Before insertion we have A -> B and X.
     * After insertion, we have A -> C -> (B, X) or A -> C -> B (where
     * C=X),  i.e. either a new, shorter node C becomes the parent of both
     * B and X or X becomes the parent of B (if X is a suffix of B). 
     * @param nodeA      The parent of the node where the split occurred
     * @param nodeB      The old node that was split
     * @param nodeC      The new, shorter node created during the split
     */
    void handleSplit(const WrappedNode& nodeA,
                     const WrappedNode& nodeB, 
                     const WrappedNode& nodeC);


    void gibbsSamplePath(const WrappedNodeList& path, 
                         const d_vec& discountPath, 
                         const d_vec& concentrationPath, 
                         const PayloadDataPath& payloadDataPath,
                         double baseProb);
    
    boost::shared_ptr<void> makeAdditionalDataPtr(void* payload, 
                                                  double discount, 
                                                  double concentration);
    

    bool checkConsistency(const WrappedNode& node, 
                          const std::list<WrappedNode>& children) const;
    
    
    class ToStringVisitor {
      public:
        ToStringVisitor(seq_type& seq, const IHPYPBaseRestaurant& restaurant);
        void operator()(const WrappedNode& n);
        
        std::ostringstream outstream;
      private:
        seq_type& seq;
        const IHPYPBaseRestaurant& restaurant;
    };
    
    class CheckConsistencyVisitor {
      public:
        CheckConsistencyVisitor(const HPYPModel& model);
        void operator()(WrappedNode& n, std::list<WrappedNode>& children);

        bool consistent;
      private:
        const HPYPModel& model;
    };


    seq_type& seq;
    boost::scoped_ptr<ContextTree> contextTree_;
    ContextTree& contextTree;
    const IAddRemoveRestaurant& restaurant;
    IParameters& parameters;
    int numTypes;
    double baseProb;


};




// class CollectStatisticsVisitor {
//   public:
// 
//     unsigned int contexts, context_symbol, tables_no_single, tables_single,
//                  tables_one_cust, customers, tables, tables_per_symbol;
// 
//     CollectStatisticsVisitor() : 
//       contexts(0), context_symbol(0), tables_no_single(0), tables_single(0), tables_one_cust(0), customers(0), tables(0), tables_per_symbol(0) {}
//     
//     void operator()(typename ContextTree::WrappedNode_t node) {
//       contexts++;
//       customers += node.payload.storage->getC();
//       tables += node.payload.storage->getT();
//       typedef std::list<e_type> key_t;
//       key_t keys = node.payload.storage->keys();
//       for (typename key_t::iterator it = keys.begin(); it!=keys.end(); ++it) {
//         context_symbol++;
//         l_type tw = node.payload.storage->getT(*it);
//         tables_per_symbol += tw;
//         if (tw == 1) {
//           tables_single++;
//           if (node.payload.storage->getC(*it) == 1) {
//             tables_one_cust++;
//           }
//         } else {
//           tables_no_single++;
//         }
//       }
//     }
// 
//     std::string toString() {
//       std::ostringstream out;
//       out << "contexts: " << contexts << ", customers: " << customers << ", context/symbol: " << context_symbol;
//       out << ", tables_single: " << tables_single << ", tables_no_single: " << tables_no_single << ", tables: " << tables;
//       out << ", tables2:" << tables_per_symbol;
// 
//       return out.str();
//     }
// 
//     static std::string header() {
//       return "i, contexts , customers, context_symbol, tables_single, tables_no_single, tables_one_cust, tables";
//     }
// 
// 
// 
//     std::vector<unsigned int> toVec() {
//       std::vector<unsigned int> out;
//       out.push_back(contexts);
//       out.push_back(customers);
//       out.push_back(context_symbol);
//       out.push_back(tables_single);
//       out.push_back(tables_no_single);
//       out.push_back(tables_one_cust);
//       out.push_back(tables);
//       return out;
//     }
// };
    

}} // namespace gatsby::libplump

#endif
