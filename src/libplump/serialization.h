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
