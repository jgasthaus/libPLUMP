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

#ifndef NODE_MANAGER_H_
#define NODE_MANAGER_H_

#include "libplump/config.h"
#include "libplump/utils.h"
#include "libplump/mini_map.h"
#include "libplump/pool.h"
#include "libplump/serialization.h"


namespace gatsby { namespace libplump {

/**
 * The node manager should abstract simple operations on nodes 
 * such as finding the children of a node or returning the payload 
 * associated with a node. For "heavy" nodes, this can just redirect to the 
 * appropriate node functions, but more clever storage schemes may be possible 
 * by using this abstraction. 
 */
class INodeManager {
  public:
    
    /**
     * Type of the unique identifier for each node managed by the NodeManager.
     */
    typedef void* NodeId;

    /**
     * Type of a map of edge label/NodeId pairs
     */
    typedef MapType<e_type, NodeId>::Type ChildMap;  

    /**
     * ChildMap iterator.
     */
    typedef ChildMap::iterator ChildMapIterator;

    virtual ~INodeManager() {}

    /**
     * Get the root node.
     *
     * The root node is created when the NodeManager is created.
     */
    virtual NodeId getRoot() const = 0;

    /**
     * Get the child of the given node with the given key.
     * 
     * Returns NULL if no such child exists.
     */
    virtual NodeId getChild(NodeId node, e_type key) const = 0;

    /**
     * Create a new child of node with the given key, start and end positions 
     * and returns a handle to the newly created node.
     */
    virtual NodeId setChild(NodeId node,
                            e_type key,
                            l_type start,
                            l_type end, 
                            void* payload = NULL) = 0;

    /**
     * Insert a node between two other nodes at the given key position.
     *
     * If parent already has a old_key-child, we create a new child of parent 
     * a make the old child a child of the new child with key new_key.
     */
    virtual NodeId insertBetween(NodeId parent,
                                 e_type old_key,
                                 l_type new_start,
                                 l_type new_end,
                                 e_type new_key) = 0;

    /**
     * Remove the child with key key from node's child list.
     */
    virtual void removeChild(NodeId node, e_type key) = 0;

    /**
     * Get the payload associated with the given node
     */
    virtual void* getPayload(NodeId node) const = 0;
    
    /**
     * Set the payload associated with the given node
     *
     * Note that this destroys the payload that was 
     * previously assigned to this node.
     */
    virtual void setPayload(NodeId node, void* payload) = 0;

    /**
     * Get the start position of the given node in the underlying string.
     */
    virtual l_type getStart(NodeId node) const = 0;

    /**
     * Get the end position (exclusive) of the given node in the underlying
     * string.
     */
    virtual l_type getEnd(NodeId node) const = 0;

    /**
     * Get all children of the given  node.
     */
    virtual ChildMap& getChildren(NodeId node) const = 0;

    /**
     * Destroy the given node and free memory.
     */
    virtual void destroyNode(NodeId node) = 0;

    /**
     * Destroy the given node and all its children.
     */
    virtual void destroyNodeRecursive(NodeId node) = 0;
};



/**
 * IPayloadFactory is an interface used by NodeManagers to obtain storage
 * space for payloads. The void pointers returned by the make function 
 * should be stored by the NodeManager together with the corresponding node, 
 * and returned to the user when getPayload is called. The user can then pass
 * this pointer to the construtor of the actual Payload.
 * 
 * Note that the pointer returned by make() is an aliasing pointer. The 
 * class implementing this interface is responsible for free-ing the memory
 * when it is destroyed. It can also free (or reuse) the memory when a call 
 * to recycle() is made.
 */
class IPayloadFactory {
  public:
    virtual void* make() const = 0;
    virtual void recycle(void*) const = 0;
    virtual void save(void*, OutArchive&) const = 0;
    virtual void* load(InArchive&) const = 0;
};



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
