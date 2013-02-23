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
    virtual bool removeCustomer(void* payloadPtr, 
                                e_type type,
                                double discount,
                                void* additionalData) const = 0;

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
