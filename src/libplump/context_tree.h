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

#ifndef CONTEXT_TREE_H_
#define CONTEXT_TREE_H_

#include <list>
#include <string>
#include <sstream>
#include <stack>
#include "libplump/config.h"
#include "libplump/node_manager.h"

namespace gatsby { namespace libplump {

/**
 * Wrapper context tree node used at the interface of the ContextTree class.
 * It contains the information stored in a node of the context tree, namely
 * its start and end position in the underlying string, as well as a pointer
 * to some associated payloads.
 *
 * WrappedNodes and lists of wrapped Nodes (WrappedNodeList) are used in the
 * interface of the ContextTree.
 */ 
class WrappedNode {
 public:
  l_type start, end, depth;

  // non-owning pointer
  void* payload;

  WrappedNode() : start(0), end(0), depth(0), payload(NULL) {}

  WrappedNode(const WrappedNode& other) {
    *this = other;
  }

  /**
   * Constructs a valid WrappedNode from the given information.
   *
   * @param start Start position in input string
   * @param end   End position
   * @param payload Payload for this node (by value!)
   * @param depth Depth of this node in the tree
   */
  WrappedNode(l_type start, l_type end, void* payload, l_type depth) : 
    start(start), end(end), depth(depth), payload(payload) {}

  WrappedNode& operator=(const WrappedNode& other) {
    this->start   = other.start;
  	this->end     = other.end;
  	this->depth   = other.depth;
  	this->payload = other.payload;
  	return *this;
  }

  std::string toString() {
    std::ostringstream outstream;
    outstream << "WrappedNode[" << start << ":" << end  << "]";
    return outstream.str();
  }
};


/**
 * List of WrappedNodes; mainly used for representing pathes through the
 * context tree, e.g. as returned by ContextTree::findLongestSuffix.
 */
typedef std::list<WrappedNode> WrappedNodeList;

/**
 * The ContextTree class implements basic operations on context trees 
 * (reverse prefix trees) over arbitrary sequences. These basic operations
 * include inserting a new node into the tree (insert()) and finding the node
 * with the longest common suffix with a given string (findLongestSuffix)
 */ 
class ContextTree {
  public:
    typedef INodeManager::NodeId NodeId;
    class DFSPathIterator; // defined further below
    class InsertionResult; // defined further below


    ContextTree(INodeManager& nodeManager, seq_type& seq);

    ~ContextTree();
   
    /**
     * Insert the given subsequence into the tree.
     */
    InsertionResult insert(l_type start, l_type end); 

    /**
     * Find the path to the node which is the longest suffix of the given
     * subsequence within the tree.
     */
    WrappedNodeList findLongestSuffix (l_type start, l_type end) const;

    /**
     * Find the path to the node which is the longest suffix of the given
     * subsequence within the tree.
     */
    std::pair<int, WrappedNodeList> findLongestSuffixVirtual (l_type start,
                                                              l_type end) const;
    
    DFSPathIterator getDFSPathIterator() const;

    std::string pathToString(WrappedNodeList path, 
                             bool printString = false, 
                             bool printPayload = true) const;

    std::string wrappedNodeToString(WrappedNode wrappedNode,
                                    bool printString = false,
                                    bool printPayload = true) const;

    std::string toString() const;


    /////////////////////// NESTED CLASSES BELOW //////////////////////////////

    /**
     * Result type of the insert() method. Contains the path to the 
     * inserted node and information about a possible split event that
     * occurred during insertion.
     */
    struct InsertionResult {
      // Possible operations performed during insertion:
      //   NO_SPLIT     no split was necessary
      //   SPLIT        split occurred and node was inserted as child
      //   SPLIT_SUFFIX split occurred and node was suffix of split node
      enum InsertAction {INSERT_ACTION_NO_SPLIT, 
                         INSERT_ACTION_SPLIT, 
                         INSERT_ACTION_SPLIT_SUFFIX};

      // the path to the inserted node (inclusive)
      WrappedNodeList path;

      // the action performed during insertion
      InsertAction action;

      // the node in which a split occurred that is no longer on the path
      WrappedNode splitChild;
    };



    class ToStringVisitor {
      public:
        ToStringVisitor(seq_type& seq);
        void operator()(WrappedNode n);
        
        std::ostringstream outstream;
      private:
        seq_type& seq;
    };

    template<typename Visitor>
    void visitDFS(Visitor& visitor) const;

    template<typename Visitor>
    void visitDFSWithChildren(Visitor& visitor) const;

    class DFSPathIterator {
      public:
        DFSPathIterator(NodeId root, const INodeManager& nm, 
                        const ContextTree& ct);
        WrappedNodeList& operator*();
        DFSPathIterator& operator ++(); // prefix version
        DFSPathIterator operator ++(int); // postfix version
        bool hasMore();
      
      private:
        class IteratorState {
          public:
            IteratorState(const INodeManager& nm, NodeId node);
            bool moreChildren();
            NodeId pop();
          private:
            friend class ContextTree::DFSPathIterator;
            INodeManager::ChildMap* childMap;
            INodeManager::ChildMap::iterator current;
            NodeId node;
            bool haveChildren;
        };

        NodeId root;
        const INodeManager& nm;
        const ContextTree& ct;
        WrappedNodeList currentPath;
        std::stack<IteratorState> iteratorStateStack;
    };


  private: 
    INodeManager& nm;
    seq_type& seq;
    INodeManager::NodeId root;

    /**
     * Determine the first position where the subsequence delimited by start 
     * and end and the subsequence s differ, 
     * i.e. the length of the longest common suffix of this sequence and s.
     * The argument offset sets an offset from where the comparison should be
     * started, e.g. start=1 makes the assumption that the longest common suffix
     * is at least of length one.
     */
    l_type suffixUntil(
        l_type thisStart,  l_type thisEnd, 
        l_type otherStart, l_type otherEnd,
        l_type offset) const; 

    /**
     * Same as above, but does not assume that 
     * thisEnd-thisStart < otherEnd - otherStart.
     */
    l_type suffixUntilCheck(
        l_type thisStart,  l_type thisEnd, 
        l_type otherStart, l_type otherEnd,
        l_type offset) const;
    
    WrappedNode wrap(NodeId node, l_type depth) const;



};

template<typename Visitor>
void ContextTree::visitDFSWithChildren(Visitor& visitor) const {
    typedef std::pair<NodeId, int> mypair;
    std::stack<mypair> mystack;
    mystack.push(mypair(root,0));
    while(!mystack.empty()) {
      mypair current = mystack.top();
      mystack.pop();
      WrappedNode n(nm.getStart(current.first), 
                    nm.getEnd(current.first),
                    nm.getPayload(current.first),
                    current.second);
      INodeManager::ChildMap children = nm.getChildren(current.first);
      std::list<WrappedNode> childList;
      for (typename INodeManager::ChildMapIterator it = children.begin();
           it != children.end(); it++) {
        mypair child;
        child.first = (*it).second;
        child.second = current.second + 1;
        mystack.push(child);
        childList.push_back(WrappedNode(
              nm.getStart(child.first),
              nm.getEnd(child.first),
              nm.getPayload(child.first),
              child.second));
      }
      visitor(n, childList);
    }
  }

template<typename Visitor>
void ContextTree::visitDFS(Visitor& visitor) const {
  typedef std::pair<NodeId, int> mypair;
  std::stack<mypair> mystack;
  mystack.push(mypair(root,0));
  while(!mystack.empty()) {
    mypair current = mystack.top();
    mystack.pop();
    WrappedNode n(nm.getStart(current.first), 
        nm.getEnd(current.first), 
        nm.getPayload(current.first), 
        current.second);
    visitor(n);
    INodeManager::ChildMap children = nm.getChildren(current.first);
    for (typename INodeManager::ChildMapIterator it = children.begin(); 
        it != children.end(); it++) {
      mypair child;
      child.first = (*it).second;
      child.second = current.second + 1;
      mystack.push(child);
    }
  }
}

}} // namespace gatsby::libplump

#endif
