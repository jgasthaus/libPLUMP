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

#include "libplump/serialization.h"

#include <vector>
#include <stack>
#include <map>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/serialization.hpp>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <iostream>
#include <fstream>


#include "libplump/node_manager.h"

namespace io = boost::iostreams;

namespace gatsby { namespace libplump {

struct ContextTreeEdge {
  ContextTreeEdge(INodeManager::NodeId id, e_type key, size_t parent) 
    : id(id), key(key), parent(parent) {}

  INodeManager::NodeId id;
  e_type key;
  size_t parent;
};

class SerializedNode {
  public:
    SerializedNode() {}
    SerializedNode(l_type start, l_type end, e_type key, size_t parent) 
      : start(start), end(end), key(key), parent(parent) {}
    l_type start;
    l_type end;
    e_type key;
    size_t parent;

  private: 
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & start;
        ar & end;
        ar & key;
        ar & parent;
    }
};


std::vector<void*> getSerializedNodeVector(const INodeManager& nm, 
                                           std::vector<SerializedNode>& nodes){
  std::vector<void*> payloads;
  std::stack<ContextTreeEdge> stack;
  nodes.push_back(SerializedNode(0,0,0,0));
  payloads.push_back(nm.getPayload(nm.getRoot()));
  INodeManager::ChildMap& rootChildren = nm.getChildren(nm.getRoot());
  for (INodeManager::ChildMapIterator it = rootChildren.begin();
       it != rootChildren.end();
       ++it) {
    stack.push(ContextTreeEdge((*it).second, (*it).first, 0));
  }
  while(!stack.empty()) {
    ContextTreeEdge current = stack.top();
    stack.pop();
    nodes.push_back(SerializedNode(nm.getStart(current.id),
                                   nm.getEnd(current.id),
                                   current.key,
                                   current.parent));
    payloads.push_back(nm.getPayload(current.id));
    for (INodeManager::ChildMapIterator it = nm.getChildren(current.id).begin();
         it != nm.getChildren(current.id).end();
         ++it) {
      stack.push(ContextTreeEdge( (*it).second, (*it).first, nodes.size()-1));
    }
  }
  return payloads;
}


std::vector<INodeManager::NodeId> insertSerializedNodeVector(
    INodeManager& nm, std::vector<SerializedNode>& nodes,
    std::vector<void*> payloads) {
  std::vector<INodeManager::NodeId> handles(nodes.size(), NULL);
  handles[0] = nm.getRoot();
  nm.setPayload(handles[0], payloads[0]);
  for (size_t i = 1; i < nodes.size(); ++i) {
    handles[i] = nm.setChild(handles[nodes[i].parent],
                             nodes[i].key,
                             nodes[i].start,
                             nodes[i].end, 
                             payloads[i]);
  }
  return handles;
}


std::vector<void*> Serializer::saveNodesToArchive(
    const INodeManager& nm, OutArchive& oa) {

  // serialize to file
  std::vector<SerializedNode> nodes;
  std::vector<void*> payloads = getSerializedNodeVector(nm, nodes);
  oa << nodes;
  return payloads;
}

void Serializer::saveNodesAndPayloadsToArchive(
    const INodeManager& nm, const IPayloadFactory& factory, OutArchive& oa) {
  std::vector<SerializedNode> nodes;
  std::vector<void*> payloads = getSerializedNodeVector(nm, nodes);
  savePayloadsToArchive(factory, payloads, oa);
  oa << nodes;
}

void Serializer::loadNodesAndPayloadsFromArchive(
    INodeManager& nm, const IPayloadFactory& factory, InArchive& ia) {
  std::vector<void*> payloads = loadPayloadsFromArchive(factory, ia);
  std::vector<SerializedNode> nodes;
  ia >> nodes;
  insertSerializedNodeVector(nm, nodes, payloads);
}

void Serializer::loadNodesFromArchive(INodeManager& nm,
                                      std::vector<void*> payloads,
                                      InArchive& ia) {
  std::vector<SerializedNode> nodes;
  ia >> nodes;
  insertSerializedNodeVector(nm, nodes, payloads);
}


std::vector<void*> Serializer::loadPayloadsFromArchive(
    const IPayloadFactory& factory, InArchive& ia) {
  size_t numPayloads; 
  ia >> numPayloads;
  std::vector<void*> payloads(numPayloads, NULL);
  for (size_t i = 0; i < numPayloads; ++i) {
    payloads[i] = factory.load(ia);
  }
  return payloads;
}

void Serializer::savePayloadsToArchive(
    const IPayloadFactory& factory, std::vector<void*> voidPayloads, 
    OutArchive& oa) {
  size_t size = voidPayloads.size();
  oa << size;
  for (size_t i = 0; i < voidPayloads.size(); ++i) {
    factory.save(voidPayloads[i], oa);
  }
}


void Serializer::saveNodesAndPayloads(const INodeManager& nodeManager,
                                      const IPayloadFactory& factory) {
  std::ofstream archiveFileStream(this->filename.c_str(),
                                  std::ios::out | std::ios::binary);
  io::filtering_streambuf<io::output> out;
  out.push(io::bzip2_compressor());
  out.push(archiveFileStream);
  std::ostream nodeStream(&out);

  OutArchive nodeArchive(nodeStream);
  Serializer::saveNodesAndPayloadsToArchive(nodeManager, factory, nodeArchive);
}


void Serializer::loadNodesAndPayloads(INodeManager& nodeManager,
                                      const IPayloadFactory& factory) {
  
  std::ifstream archiveFileStream(this->filename.c_str(),
                                  std::ios::in | std::ios::binary);
  io::filtering_streambuf<io::input> in;
  in.push(io::bzip2_decompressor());
  in.push(archiveFileStream);
  std::istream nodeStream(&in);
  
  InArchive nodeArchive(nodeStream);

  Serializer::loadNodesAndPayloadsFromArchive(nodeManager,
                                              factory,
                                              nodeArchive);
}

}}
