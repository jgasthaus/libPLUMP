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

