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

#ifndef SERIALIZATION_H_
#define SERIALIZATION_H_

#include "libplump/config.h"

namespace boost { namespace archive {
  // forward declaration
class binary_iarchive;
class binary_oarchive;
}}

namespace gatsby { namespace libplump {

typedef boost::archive::binary_iarchive InArchive;
typedef boost::archive::binary_oarchive OutArchive;

class INodeManager;
class IPayloadFactory;


class Serializer {
  public:
    Serializer(std::string filename) 
        : filename(filename) {}

    
    void saveNodesAndPayloads(const INodeManager& nm,
                              const IPayloadFactory& factory);

    void loadNodesAndPayloads(INodeManager& nm,
                              const IPayloadFactory& factory);

  private:
    static void loadNodesAndPayloadsFromArchive(INodeManager& nm,
                                                const IPayloadFactory& factory,
                                                InArchive& ia);

    static std::vector<void*> saveNodesToArchive(const INodeManager& nodeManager,
                                                 OutArchive& oa);

    static void loadNodesFromArchive(INodeManager& nodeManager,
                                     std::vector<void*> payloads,
                                     InArchive& ia);

    static std::vector<void*> loadPayloadsFromArchive(
        const IPayloadFactory& factory, InArchive& ia);


    static void savePayloadsToArchive(const IPayloadFactory& factory, 
                                      std::vector<void*> voidPayloads,
                                      OutArchive& oa);
    
    static void saveNodesAndPayloadsToArchive(
        const INodeManager& nm, const IPayloadFactory& factory,
        OutArchive& oa);

    std::string filename;
};

}} // namespace gatsby::libplump

#endif
