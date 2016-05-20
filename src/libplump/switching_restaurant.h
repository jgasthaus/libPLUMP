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

#ifndef SWITCHING_RESTAURANT_H
#define SWITCHING_RESTAURANT_H

#include <boost/scoped_ptr.hpp>

#include "libplump/config.h"
#include "libplump/hpyp_restaurant_interface.h"
#include "libplump/serialization.h"

namespace gatsby { namespace libplump {

class SwitchingRestaurant : public IAddRemoveRestaurant {
  public:
    SwitchingRestaurant(IAddRemoveRestaurant* switchedRestaurant,
                        int numSlots);

    ~SwitchingRestaurant() {}

    l_type getC(void* payloadPtr, e_type type) const;
    l_type getC(void* payloadPtr) const;
    l_type getT(void* payloadPtr, e_type type) const;
    l_type getT(void* payloadPtr) const;
    
    double computeProbability(void*  payloadPtr,
                              e_type type, 
                              double parentProbability,
                              double discount, 
                              double concentration) const;

    TypeVector getTypeVector(void* payloadPtr) const;
    
    const IPayloadFactory& getFactory() const;
    
    void updateAfterSplit(void* longerPayload, 
                          void* shorterPayload, 
                          double discountBeforeSplit, 
                          double discountAfterSplit, 
                          bool parentOnly = false) const;

    std::string toString(void* payloadPtr) const;
    
    bool checkConsistency(void* payloadPtr) const;

    double addCustomer(void*  payloadPtr, 
                     e_type type, 
                     double parentProbability, 
                     double discount, 
                     double concentration,
                     void*  additionalData = NULL,
                     double count = 1) const;


    double removeCustomer(void* payloadPtr, 
                        e_type type,
                        double discount,
                        void* additionalData,
                        double count = 1) const;

    void* createAdditionalData(void* payloadPtr, 
                               double discount, 
                               double concentration) const;

    void freeAdditionalData(void* additionalData) const;

    bool selectSlot(int slot);

  private:
    void* getCurrent(void* payloadPtr) const {
      return ((Payload*)payloadPtr)->payloads[this->currentSlot];
    }
    

    struct Payload {
      std::vector<void*> payloads;
    };

    class PayloadFactory : public IPayloadFactory {
      public:
        PayloadFactory(const SwitchingRestaurant& switchingRestaurant)
            : switchingRestaurant(switchingRestaurant) {}
        
        void* make() const;
        void recycle(void* payloadPtr) const;
        void save(void* payloadPtr, OutArchive& oa) const;
        void* load(InArchive& ia) const;

      private: 
        const SwitchingRestaurant& switchingRestaurant;
    };


    const PayloadFactory payloadFactory;
    boost::scoped_ptr<IAddRemoveRestaurant> switchedRestaurant;
    int numSlots;
    int currentSlot;
};


}} // namespace gatsby::libplump
#endif

