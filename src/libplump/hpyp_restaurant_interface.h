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

#ifndef HPYP_RESTAURANT_INTERFACE_H
#define HPYP_RESTAURANT_INTERFACE_H


#include "libplump/config.h"
#include "libplump/node_manager.h" // for IPayloadFactory

namespace gatsby { namespace libplump {

////////////////////////////////////////////////////////////////////////////////
////////////////////////   INTERFACES   ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class IHPYPBaseRestaurant {
  public:
    typedef std::vector<e_type> TypeVector;
    typedef TypeVector::iterator TypeVectorIterator;
    virtual ~IHPYPBaseRestaurant() {}
    virtual l_type getC(void* payloadPtr, e_type type) const = 0;
    virtual l_type getC(void* payloadPtr) const = 0;
    virtual l_type getT(void* payloadPtr, e_type type) const = 0;
    virtual l_type getT(void* payloadPtr) const = 0;
    virtual double computeProbability(void*  payloadPtr,
                                      e_type type, 
                                      double parentProbability,
                                      double discount, 
                                      double concentration) const = 0;
    virtual TypeVector getTypeVector(void* payloadPtr) const = 0;
    virtual const IPayloadFactory& getFactory() const = 0;
    virtual void updateAfterSplit(void* longerPayload, 
                                  void* shorterPayload, 
                                  double discountBeforeSplit, 
                                  double discountAfterSplit,
                                  bool parentOnly = false) const = 0;
    virtual std::string toString(void* payloadPtr) const = 0;
    virtual bool checkConsistency(void* payloadPtr) const = 0;
};


class IAddRestaurant : public IHPYPBaseRestaurant {
  public:
    virtual ~IAddRestaurant() {}
    virtual double addCustomer(void*  payloadPtr, 
                             e_type type, 
                             double parentProbability, 
                             double discount, 
                             double concentration,
                             void*  additionalData = NULL,
                             double count = 1) const = 0;
};


class IAddRemoveRestaurant : public IAddRestaurant {
  public:
    virtual ~IAddRemoveRestaurant() {}
    virtual double removeCustomer(void* payloadPtr, 
                                e_type type,
                                double discount,
                                void* additionalData,
                                double count = 1) const = 0;

    virtual void* createAdditionalData(void* payloadPtr, 
                                       double discount, 
                                       double concentration) const = 0;

    /**
     * Free the memory allocated for additionalData. 
     *
     * This should be called for every piece of additionalData 
     * created using createAdditionalData.
     */
    virtual void freeAdditionalData(void* additionalData) const = 0;
};

}} // namespace gatsby::libplump

#endif
