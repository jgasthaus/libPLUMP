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

#include "libplump/context_tree.h"

#include <sstream>
#include "libplump/subseq.h"
#include "libplump/utils.h"

namespace gatsby { namespace libplump {

ContextTree::ContextTree(INodeManager& nodeManager, seq_type& seq) 
    : nm(nodeManager), seq(seq) {
  root = nm.getRoot();
}


ContextTree::~ContextTree() {
  this->nm.destroyNodeRecursive(this->root);
}


l_type ContextTree::suffixUntil(l_type thisStart,
                                l_type thisEnd,
                                l_type otherStart,
                                l_type otherEnd,
                                l_type offset) const {
  //l_type l = ((tend-start) < s.length()) ? (end-start) : s.length();
  l_type l = thisEnd - thisStart; // assume this seq is always shorter
  l_type i;
  for (i = offset; i<l; i++) {
    if (seq[thisEnd - 1 - i] != seq[otherEnd - 1 - i])
      break;
  }
  return i;
}


l_type ContextTree::suffixUntilCheck(l_type thisStart,
                                     l_type thisEnd,
                                     l_type otherStart,
                                     l_type otherEnd,
                                     l_type offset) const {
  l_type l = std::min(thisEnd-thisStart, otherEnd-otherStart);
  l_type i;
  for (i = offset; i<l; i++) {
    if (seq[thisEnd - 1 - i] != seq[otherEnd - 1 - i])
      break;
  }
  return i;
}


ContextTree::InsertionResult ContextTree::insert(l_type start, l_type end) 
{
  tracer << "ContextTree::insert(" << start << ", " << end << ")" << std::endl;
  ContextTree::InsertionResult result;
  l_type offset = 0;
  WrappedNodeList& path = result.path;
  NodeId current = root;
  NodeId parent;
  e_type parentKey = 0;
  l_type depth = 0;
  while (true) {
    l_type curStart = nm.getStart(current);
    l_type curEnd = nm.getEnd(current);
    l_type curLength = curEnd - curStart;

    // wrap the current node and put it in the list
    WrappedNode n(curStart, curEnd, nm.getPayload(current), depth);
    path.push_back(n);

    // determine the length l of the longest common suffix with s
    l_type longestSuffixLen = suffixUntilCheck(curStart, curEnd,
                                               start, end, offset); 
    if (curLength == longestSuffixLen) { // if context is fully consumed
      if (curLength == end - start) { // context is already in the tree!
        result.action = InsertionResult::INSERT_ACTION_NO_SPLIT;
        break;
      }
      // determine key for the child pointer
      e_type key = seq[end - 1 - curLength]; 
      NodeId child = nm.getChild(current,key);
      if (child != NULL) {
        offset = curLength;
        parent = current;
        current = child;
        parentKey = key;
        depth++;
        continue;
      } else {
        child = nm.setChild(current, key, start, end);
        // The payload of the parent may have changed!
        path.back().payload = nm.getPayload(current);
        WrappedNode wrapped_child(start, end, nm.getPayload(child), depth + 1);
        path.push_back(wrapped_child);
        result.action = InsertionResult::INSERT_ACTION_NO_SPLIT;
        break;
      }
    } else { // context not fully consumed: split node and insert
      // We are trying to insert a context X into the tree and
      // have reached a node B such that C is a common suffix
      // of B and C, and this suffix is shorter than B.
      // We now need to insert a node for C into the tree 
      // and make B a child of C. Also, if C is a proper 
      // suffix of X, we need to insert X as a child of C.

      // start position for the common suffix C
      l_type shorterStart = curEnd - longestSuffixLen;

      // key for B in C
      e_type oldNodeKey = seq[curEnd - 1 - longestSuffixLen];
      NodeId newParent = nm.insertBetween(parent,
                                          parentKey,
                                          shorterStart,
                                          curEnd,
                                          oldNodeKey);

      // return node B
      result.splitChild = path.back();
      // pop node B off the path
      path.pop_back(); 

      // push node C onto the path
      path.push_back(WrappedNode(shorterStart,
                                 curEnd,
                                 nm.getPayload(newParent),
                                 depth));

      // if C is a proper suffix of X, insert X as a child of C
      // and push C onto the path
      if (end - start > curEnd - shorterStart) {
        e_type childKey = seq[end - longestSuffixLen - 1];
        NodeId child = nm.setChild(newParent, childKey, start, end);
        path.push_back(WrappedNode(start, end, nm.getPayload(child),depth+1));
        result.action = InsertionResult::INSERT_ACTION_SPLIT;
      } else {   
        result.action = InsertionResult::INSERT_ACTION_SPLIT_SUFFIX;
      }
      break;
    }
  }
  return result;
}


/**
 * Find the path to the node which is the longest suffix of the given
 * subsequence within the tree.
 */
WrappedNodeList ContextTree::findLongestSuffix (l_type start, l_type end) const {
  l_type offset = 0;
  WrappedNodeList path;
  NodeId current = root;
  NodeId parent;
  e_type parentKey = 0;
  l_type depth = 0;
  bool done = false;
  while (!done) {
    l_type curStart = nm.getStart(current);
    l_type curEnd = nm.getEnd(current);
    l_type curLength = curEnd - curStart;
    // wrap the current node and put it in the list
    WrappedNode n(curStart, curEnd, nm.getPayload(current), depth);
    path.push_back(n);
    // determine the length l of the longest common suffix with s
    l_type longestSuffixLen = suffixUntilCheck(curStart, curEnd, 
                                               start, end, offset);
    if (curLength == longestSuffixLen) { // if context is fully consumed
      if (curLength == end - start) { // no more input to consume
        done = true;
      } else {
        // determine key for the child pointer
        e_type key = seq[end - 1 - curLength];
        current = nm.getChild(current,key);
        if (current == NULL) {
          done = true;
          break;
        }
        offset = curLength;
        parent = current;
        parentKey = key;
        depth++;
      }
    } else {
      // drop last node from path as it would have to be split
      path.pop_back();
      done = true;
    }
  }
  return path;
}


/**
 * Find the path to the node which is the longest suffix of the given
 * subsequence within the tree.
 */
std::pair<int,WrappedNodeList> ContextTree::findLongestSuffixVirtual(
    l_type start, l_type end) const {
  std::pair<int,WrappedNodeList> ret;
  ret.first = 0;
  l_type offset = 0;
  WrappedNodeList& path = ret.second;
  NodeId current = root;
  NodeId parent;
  e_type parentKey = 0;
  l_type depth = 0;
  bool done = false;
  while (!done) {
    l_type curStart = nm.getStart(current);
    l_type curEnd = nm.getEnd(current);
    l_type curLength = curEnd - curStart;
    // wrap the current node and put it in the list
    WrappedNode n(curStart, curEnd, nm.getPayload(current), depth);
    path.push_back(n);
    // determine the length l of the longest common suffix with s
    l_type longestSuffixLen = suffixUntilCheck(curStart, curEnd, 
                                               start, end, offset);
    if (curLength == longestSuffixLen) { // if context is fully consumed
      if (curLength == end - start) { // no more input to consume
        done = true;
      } else {
        // determine key for the child pointer
        e_type key = seq[end - 1 - curLength];
        current = nm.getChild(current,key);
        if (current == NULL) {
          done = true;
          break;
        }
        offset = curLength;
        parent = current;
        parentKey = key;
        depth++;
      }
    } else {
      ret.first = longestSuffixLen; // length of 
      done = true;
    }
  }
  return ret;
}


WrappedNode ContextTree::wrap(NodeId node, l_type depth) const {
  return WrappedNode(nm.getStart(node),
                     nm.getEnd(node),
                     nm.getPayload(node),
                     depth); 
}


ContextTree::DFSPathIterator ContextTree::getDFSPathIterator() const {
  return ContextTree::DFSPathIterator(root,nm,*this);
}


std::string ContextTree::pathToString(WrappedNodeList path,
                                      bool printString,
                                      bool printPayload) const {
  std::ostringstream outstream;
  for(WrappedNodeList::const_iterator it=path.begin();it!=path.end();it++) {
    outstream << "start: " << it->start << ", end: " << it->end;
    if (printString) {
      outstream << ", string: " << SubSeq::toString(it->start, it->end, seq);
    }
    //if (printPayload) {
    //    outstream << ", payload: " << it->payload.toString();
    //}
    outstream<< std::endl;
  }
  return outstream.str();
}


std::string ContextTree::wrappedNodeToString(WrappedNode wrappedNode,
                                             bool printString,
                                             bool printPayload) const {
  std::ostringstream outstream;
  outstream << "start: " << wrappedNode.start << ", end: " << wrappedNode.end;
  if (printString) {
    outstream << ", string: " 
              << SubSeq::toString(wrappedNode.start, wrappedNode.end, seq);
  }
  if (printPayload) {
    //outstream << ", payload: " << wrappedNode.payload.toString(); 
    //TODO: do we want this feature?
  }
  outstream<< std::endl;
  return outstream.str();
}


std::string ContextTree::toString() const {
  ContextTree::ToStringVisitor visitor(seq);
  visitDFS(visitor);
  return visitor.outstream.str();
}



////////////////////////////////////////////////////////////////////////////////
////////////   ContextTree::ToStringVisitor   //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
ContextTree::ToStringVisitor::ToStringVisitor(seq_type& seq) : seq(seq) {}
void ContextTree::ToStringVisitor::operator()(WrappedNode n) {
  for (l_type i = 0; i < n.depth; ++i) {
    this->outstream << " "; 
  }
  this->outstream << SubSeq::toString(n.start, n.end, seq);
  this->outstream << std::endl;
}
       
////////////////////////////////////////////////////////////////////////////////
////////////   ContextTree::DFSPathIterator   //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
ContextTree::DFSPathIterator::DFSPathIterator(NodeId root, 
                                            const INodeManager& nm, 
                                            const ContextTree& ct) 
    : root(root), nm(nm), ct(ct) {
  IteratorState r(this->nm, this->root);
  this->iteratorStateStack.push(r);
  this->currentPath.push_back(ct.wrap(this->root, 0));
  int j = 1;
  while (this->iteratorStateStack.top().moreChildren()) {
    IteratorState c(this->nm,this->iteratorStateStack.top().pop());
    this->iteratorStateStack.push(c);
    this->currentPath.push_back(ct.wrap(c.node,j));
    ++j;
  }
}

WrappedNodeList& ContextTree::DFSPathIterator::operator*() {
  return this->currentPath;
}

// prefix version
ContextTree::DFSPathIterator& ContextTree::DFSPathIterator::operator ++() {
  if (this->iteratorStateStack.empty()) {
    return *this;
  }
  if (!this->iteratorStateStack.top().moreChildren()) {
    // no more children, back up
    this->iteratorStateStack.pop();
    this->currentPath.pop_back();
  }

  if (this->iteratorStateStack.empty())
    return *this;

  while(this->iteratorStateStack.top().moreChildren()) {
    IteratorState c(this->nm,this->iteratorStateStack.top().pop());
    this->iteratorStateStack.push(c);
    this->currentPath.push_back(ct.wrap(c.node,0));
  }
  return *this;
}

// postfix version
ContextTree::DFSPathIterator ContextTree::DFSPathIterator::operator ++(int) {
  ContextTree::DFSPathIterator old(*this);
  this->operator++();
  return old;
}

bool ContextTree::DFSPathIterator::hasMore() {
  return !this->iteratorStateStack.empty();
}


////////////////////////////////////////////////////////////////////////////////
////////////   ContextTree::DFSPathIterator::IteratorState   ///////////////////
////////////////////////////////////////////////////////////////////////////////
ContextTree::DFSPathIterator::IteratorState::IteratorState(
    const INodeManager& nm, NodeId node) : node(node), haveChildren(false) {
  INodeManager::ChildMap& children = nm.getChildren(node);
  haveChildren = !children.empty();
  if (haveChildren) {
    this->childMap = &children;
    this->current = this->childMap->begin();
  }
}

bool ContextTree::DFSPathIterator::IteratorState::moreChildren() {
  if (haveChildren) {
    return (this->current != this->childMap->end());
  } else {
    return false;
  }
}

ContextTree::NodeId ContextTree::DFSPathIterator::IteratorState::pop() {
  ContextTree::NodeId ret = (*(this->current)).second;
  ++(this->current);
  return ret;
}


}} // namespace gatsby::libplump
