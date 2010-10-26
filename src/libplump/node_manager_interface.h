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

#ifndef NODE_MANAGER_INTERFACE_H_
#define NODE_MANAGER_INTERFACE_H_

#include "libplump/config.h"
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

}} // namespace gatsby::libplump

#endif
