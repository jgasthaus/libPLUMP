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

#ifndef NODE_MANAGER_H_
#define NODE_MANAGER_H_

#include "libplump/config.h"
#include "libplump/utils.h"
#include "libplump/node_manager_interface.h"
#include "libplump/mini_map.h"
#include "libplump/pool.h"
#include "libplump/serialization.h"

namespace gatsby { namespace libplump {

/**
 * A SimpleNodeManager implements the NodeManager functionality by 
 * storing a map from keys to pointers to children in each node. 
 * It manages memory using a Boost pool for efficient creation/deletion of
 * nodes.
 * It treats all nodes the same (i.e. all have Payloads).
 */
class SimpleNodeManager : public INodeManager {

  public:

    SimpleNodeManager(const IPayloadFactory& payloadFactory) 
      :  payloadFactory(payloadFactory), root(createNode(0,0)){} 

    ~SimpleNodeManager() {
      this->destroyNode(this->root);
    }

    NodeId getRoot() const {
      return root;
    }

    /**
     * Get the child of this node with the given key.
     */
    NodeId getChild(NodeId node, e_type key) const {
      ChildMapIterator child = static_cast<Node*>(node)->children.find(key);
      if (child != static_cast<Node*>(node)->children.end()) {
        return (*child).second;
      } else {
        return NULL;
      } 
    }

    /**
     * Create a new child of node with the given key, start and end positions 
     * and returns a handle to the newly created node.
     *
     * TODO: This method does not check if a child with this key
     *       already exists. If it does, we just drop that old pointer and 
     *       have ourselves a nice memory leak.
     */
    NodeId setChild(NodeId node, e_type key, l_type start, l_type end,
                    void* payload = NULL) {
      NodeId newNode = createNode(start, end, payload);
      static_cast<Node*>(node)->children[key] = newNode;
      return newNode;
    }

    /**
     * Insert a node between two other nodes at the given key position.
     *
     * If parent already has a old_key-child, we create a new child of parent 
     * a make the old child a child of the new child with key new_key.
     */
    NodeId insertBetween(NodeId parent, e_type oldKey,
                         l_type newStart, l_type newEnd,  e_type newKey) {
      // we can safely assume that parent is an inner node
      NodeId oldChild = static_cast<Node*>(parent)->children[oldKey];
      // create intermediate node
      Node* newParent = createNode(newStart, newEnd); 
      // make old node child of intermediate
      newParent->children[newKey] = oldChild; 
      // make intermediate node child of parent
      static_cast<Node*>(parent)->children[oldKey] = newParent;
      return newParent;
    }

    /**
     * Remove the child with key key from node's child list.
     * XXX: Not implemented -- and not needed?
     */
    void removeChild(NodeId node, e_type key) {

    }

    /**
     * Get a reference to the payload associated with the
     * given node.
     */
    void* getPayload(NodeId node) const {
      return static_cast<Node*>(node)->payload;
    }
    
    void setPayload(NodeId node, void* payload) {
      this->payloadFactory.recycle(static_cast<Node*>(node)->payload);
      static_cast<Node*>(node)->payload = payload;
    }

    l_type getStart(NodeId node) const {
      return static_cast<Node*>(node)->start;
    }

    l_type getEnd(NodeId node) const {
      return static_cast<Node*>(node)->end;
    }

    ChildMap& getChildren(NodeId node) const {
      return static_cast<Node*>(node)->children;
    }
    
    
    /**
     * Destroy the node with the given handle, i.e.
     * call its destructor and free the allocated memory.
     */
    void destroyNode(NodeId node) {
      if (node != NULL) {
        this->payloadFactory.recycle(static_cast<Node*>(node)->payload);
        delete static_cast<Node*>(node);
      }
    }


    /**
     * Destroy the given node and all its children.
     */
    void destroyNodeRecursive(NodeId node) {
      if (node != NULL) {
        for (ChildMapIterator it = static_cast<Node*>(node)->children.begin();
             it != static_cast<Node*>(node)->children.end(); ++it) {
          this->destroyNodeRecursive((*it).second);
        }
        this->destroyNode(static_cast<Node*>(node));
        if (node == this->root) {
          // if we have destroyed the root, create a new one
          this->root = createNode(0,0);
        }
      }
    }



  private:
    class Node : public PoolObject<Node> {
      public:
        l_type start;
        l_type end;
        void* payload;
        ChildMap children;

        Node(l_type start, l_type end, void* payload) 
          : start(start), end(end), payload(payload), children() {}

        Node() : start(0), end(0), payload(NULL), children() {}
    };


    /**
     * Create a new, initally unreferenced node and return its handle.
     *
     * The caller takes ownership of the returned object and should 
     * call destroyNode or destroyNodeRecursive when done using it.
     */
    Node* createNode(l_type start, l_type end, void* payload = NULL) {
      //NodeId node(new (node_pool.malloc()) Node());
      Node* node = new Node();
      node->start = start;
      node->end = end;
      if (payload != NULL) {
        node->payload = payload;
      } else {
        node->payload = this->payloadFactory.make();
      }
      return node;
    }
    

    const IPayloadFactory& payloadFactory;
    Node* root;

    DISALLOW_COPY_AND_ASSIGN(SimpleNodeManager);
};

// 
// template<typename _OtherNodeManager>
// class DebugProxyNodeManager {
// 
//     private:
// 
// 		_OtherNodeManager nm;
//         unsigned int num_getChild, 
//                      num_setChild,
//                      num_insertBetween,
//                      num_removeChild,
//                      num_getPayload, 
//                      num_getStart,
//                      num_getEnd,
//                      num_getChildren;
// 
// 
//     public:
//         typedef typename _OtherNodeManager::Payload Payload;
// 
//         typedef typename _OtherNodeManager::NodeId NodeId;
//         typedef typename _OtherNodeManager::child_map_t child_map_t;
// 
//         static const NodeId NO_CHILD;
// 
//         DebugProxyNodeManager() {
//             num_getChild = num_setChild = num_insertBetween = num_removeChild = num_getPayload = num_getStart = num_getEnd = num_getChildren = 0;
//         }
// 
// 
//        // NodeId createNode(l_type start, l_type end) {
//        //     NodeId node = nm.createNode(start, end);
//        // 	std::cerr << node << " = createNode(" << start << ", " << end << ")" << std::endl;
//        //     return node;
//        // }
// 
//        // void destroyNode(NodeId& node) {
//        // 	std::cerr << "destroyNode(" << node << ")" << std::endl;
//        // 	nm.destroyNode(node);
//        // }
//         
//         NodeId getRoot() const {
//             std::cerr << "getRoot()" << std::endl;
//             return nm.getRoot();
//         }
// 
//         std::pair<bool,NodeId> getChild(NodeId& node, e_type key) {
//             num_getChild++;
//         	std::cerr << "getChild(" << node << ", " << key << ")";
//             std::pair<bool, NodeId> child = nm.getChild(node,key);
//         	std::cerr << ", returned: " << child.first << ", " << child.second << std::endl;
//             return child;
//         }
// 
//         NodeId setChild(NodeId& node, e_type key, l_type start, l_type end){
//             num_setChild++;
//         	std::cerr << "setChild(" << node << ", " << key << ", " << start<< ", " << end << ")";
//             NodeId new_child = nm.setChild(node,key,start,end);
//             std::cerr << ", returned: " << new_child << std::endl;
//         	return new_child;
//         }
//         
//         NodeId insertBetween(NodeId& parent, e_type old_key, l_type new_start, l_type new_end,  e_type new_key) {
//             num_insertBetween++;
//             std::cerr << "insertBetween(" << parent << ", " << old_key << ", " << new_start << ", " << new_end << ", " << new_key << ")";
//             NodeId new_node = nm.insertBetween(parent,old_key,new_start,new_end,new_key);
//             std::cerr << ", returned: " << new_node << std::endl;
//             return new_node;
//         }
// 
//         void removeChild(NodeId& node, e_type key) {
//             num_removeChild++;
//         	std::cerr << "removeChild(" << node << ", " << key << ")" << std::endl;
//         	nm.removeChild(node,key);
//         }
// 
//         Payload getPayload(NodeId& node){
//             num_getPayload++;
//         	std::cerr << "getPayload(" << node << ")" << std::endl;
//             return nm.getPayload(node);
//         }
// 
//         l_type getStart(NodeId& node) {
//             num_getStart++;
//         	std::cerr << "getStart(" << node << ")";
//             l_type start = nm.getStart(node);
//             std::cerr << ", returned: " << start << std::endl;
//             return start;
//         }
// 
//         l_type getEnd(NodeId& node) {
//             num_getEnd++;
//         	std::cerr << "getEnd(" << node << ")";
//         	l_type end = nm.getEnd(node);
//         	std::cerr << ", returned: " << end << std::endl;
//             return end;
//         }
// 
//         std::pair<bool, child_map_t*> getChildren(NodeId& node) {
//             num_getChildren++;
//         	std::cerr << "getChildren(" << node << ")" << std::endl;
//             return nm.getChildren(node);
//        }
// 
//        void printStats() {
//            std::cerr <<   "num_getChild       = " <<  num_getChild     << std::endl; 
//            std::cerr <<   "num_setChild       = " <<  num_setChild     << std::endl; 
//            std::cerr <<   "num_insertBetween  = " <<  num_insertBetween<< std::endl;
//            std::cerr <<   "num_removeChild    = " <<  num_removeChild  << std::endl; 
//            std::cerr <<   "num_getPayload     = " <<  num_getPayload   << std::endl; 
//            std::cerr <<   "num_getStart       = " <<  num_getStart     << std::endl; 
//            std::cerr <<   "num_getEnd         = " <<  num_getEnd       << std::endl; 
//            std::cerr <<   "num_getChildren    = " <<  num_getChildren  << std::endl;
//        }
// 
// };

}} // namespace gatsby::libplump
#endif
