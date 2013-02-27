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

#include "libplump/switching_restaurant.h"

#include <vector>
#include <sstream>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

namespace gatsby { namespace libplump {

////////////////////////////////////////////////////////////////////////////////
//////////////////////   class SwitchingRestaurant   ///////////////////////////
////////////////////////////////////////////////////////////////////////////////

SwitchingRestaurant::SwitchingRestaurant(
    IAddRemoveRestaurant* switchedRestaurant, int numSlots) 
    : payloadFactory(*this), switchedRestaurant(switchedRestaurant), 
      numSlots(numSlots), currentSlot(0) {}


l_type SwitchingRestaurant::getC(void* payloadPtr, e_type type) const {
  return this->switchedRestaurant->getC(getCurrent(payloadPtr), type);
}


l_type SwitchingRestaurant::getC(void* payloadPtr) const {
  return this->switchedRestaurant->getC(getCurrent(payloadPtr));
}


l_type SwitchingRestaurant::getT(void* payloadPtr, e_type type) const {
  return this->switchedRestaurant->getT(getCurrent(payloadPtr), type);
}


l_type SwitchingRestaurant::getT(void* payloadPtr) const {
  return this->switchedRestaurant->getT(getCurrent(payloadPtr));
}


double SwitchingRestaurant::computeProbability(void*  payloadPtr,
                          e_type type, 
                          double parentProbability,
                          double discount, 
                          double concentration) const {
  return this->switchedRestaurant->computeProbability(getCurrent(payloadPtr),
                                                     type,
                                                     parentProbability,
                                                     discount,
                                                     concentration);
}


IHPYPBaseRestaurant::TypeVector SwitchingRestaurant::getTypeVector(
    void* payloadPtr) const {
  return this->switchedRestaurant->getTypeVector(getCurrent(payloadPtr));
}


const IPayloadFactory& SwitchingRestaurant::getFactory() const {
  return this->payloadFactory;
}


void SwitchingRestaurant::updateAfterSplit(void* longerPayload, 
                                           void* shorterPayload, 
                                           double discountBeforeSplit,
                                           double discountAfterSplit,
                                           bool parentOnly) const {
  for (int i = 0; i < this->numSlots; ++i) {
    this->switchedRestaurant->updateAfterSplit(
        ((Payload*)longerPayload)->payloads[i],
        ((Payload*)shorterPayload)->payloads[i], 
        discountBeforeSplit,
        discountAfterSplit,
        parentOnly);
  }
}


std::string SwitchingRestaurant::toString(void* payloadPtr) const {
  std::ostringstream out;
  Payload* p = (Payload*)payloadPtr;
  out << "[";
  for (int i = 0; i < this->numSlots; ++i) {
    out << i << ":" 
        << this->switchedRestaurant->toString(p->payloads[i])
        << ", ";
  }
  out << "]";
  return out.str();
}


bool SwitchingRestaurant::checkConsistency(void* payloadPtr) const {
  bool consistent = true;
  for (int i = 0; i < this->numSlots; ++i) {
    consistent = consistent && this->switchedRestaurant->checkConsistency(
        ((Payload*)payloadPtr)->payloads[i]);
  }
  return consistent;
}


double SwitchingRestaurant::addCustomer(void*  payloadPtr, 
                 e_type type, 
                 double parentProbability, 
                 double discount, 
                 double concentration,
                 void*  additionalData,
                 double count) const {
  return this->switchedRestaurant->addCustomer(getCurrent(payloadPtr),
                                              type, 
                                              parentProbability, 
                                              discount, 
                                              concentration,
                                              additionalData,
                                              count);
}


double SwitchingRestaurant::removeCustomer(void* payloadPtr, 
                    e_type type,
                    double discount,
                    void* additionalData,
                    double count) const {
  return this->switchedRestaurant->removeCustomer(getCurrent(payloadPtr),
                                                 type,
                                                 discount,
                                                 additionalData,
                                                 count);
}


void* SwitchingRestaurant::createAdditionalData(void* payloadPtr, 
                           double discount, 
                           double concentration) const {
  return this->switchedRestaurant->createAdditionalData(getCurrent(payloadPtr),
                                                       discount,
                                                       concentration);

}


void SwitchingRestaurant::freeAdditionalData(void* additionalData) const {
  return this->switchedRestaurant->freeAdditionalData(additionalData);
}


bool SwitchingRestaurant::selectSlot(int slot) {
  if (slot >= 0 && slot < this->numSlots) {
    this->currentSlot = slot;
    return true;
  } else {
    return false;
  }
}


void* SwitchingRestaurant::PayloadFactory::make() const {
  Payload* p = new Payload();
  for (int i = 0; i < this->switchingRestaurant.numSlots; ++i) {
    p->payloads.push_back(
        this->switchingRestaurant.switchedRestaurant->getFactory().make());
  }
  return p;
}


void SwitchingRestaurant::PayloadFactory::recycle(void* payloadPtr) const {
  for (int i = 0; i < this->switchingRestaurant.numSlots; ++i) {
    this->switchingRestaurant.switchedRestaurant->getFactory().recycle(
        ((Payload*)payloadPtr)->payloads[i]);
  }

  delete (Payload*)payloadPtr;
}

void SwitchingRestaurant::PayloadFactory::save(void* payloadPtr, 
                                               OutArchive& oa) const {

  size_t size = ((Payload*)payloadPtr)->payloads.size();
  oa << size;
  for (size_t i = 0; i < size; ++i) {
    switchingRestaurant.switchedRestaurant->getFactory().save(
        ((Payload*)payloadPtr)->payloads[i], oa);
  }
}

void* SwitchingRestaurant::PayloadFactory::load(InArchive& ia) const {
  Payload* p = new Payload();
  size_t size;
  ia >> size;
  p->payloads.resize(size, NULL);
  for (size_t i = 0; i < size; ++i) {
    p->payloads[i] = 
        switchingRestaurant.switchedRestaurant->getFactory().load(ia);
  }
  return p;
}

}} // namespace gatsby::libplump
