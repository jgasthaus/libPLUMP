#include "libplump/hpyp_restaurants.h"

#include <sstream>
#include <boost/scoped_ptr.hpp>

// serialization stuff
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/serialization.hpp>

#include "libplump/random.h"
#include "libplump/pyp_sample.h"
#include "libplump/stirling.h"


namespace gatsby { namespace libplump {

////////////////////////////////////////////////////////////////////////////////
//////////////////////   class SimpleFullRestaurant   //////////////////////////
////////////////////////////////////////////////////////////////////////////////

l_type SimpleFullRestaurant::getC(void* payloadPtr, e_type type) const {
  Payload& payload = *((Payload*)payloadPtr);
  Payload::TableMap::iterator it = payload.tableMap.find(type);
  if (it != payload.tableMap.end()) {
    return it->second.first; //cw
  } else {
    return 0;
  }
}


l_type SimpleFullRestaurant::getC(void* payloadPtr) const {
  return ((Payload*)payloadPtr)->sumCustomers;
}


l_type SimpleFullRestaurant::getT(void* payloadPtr, e_type type) const {
  Payload& payload = *((Payload*)payloadPtr);
  Payload::TableMap::iterator it = payload.tableMap.find(type);
  if (it != payload.tableMap.end()) {
    return it->second.second.size();
  } else {
    return 0;
  }
}


l_type SimpleFullRestaurant::getT(void* payloadPtr) const {
  return ((Payload*)payloadPtr)->sumTables;
}


double SimpleFullRestaurant::computeProbability(void*  payloadPtr,
                                                e_type type, 
                                                double parentProbability,
                                                double discount, 
                                                double concentration) const {
  Payload& payload = *((Payload*)payloadPtr);
  int cw = 0;
  int tw = 0;
  Payload::TableMap::iterator it = payload.tableMap.find(type);
  if (it != payload.tableMap.end()) {
    cw = (*it).second.first;
    tw = (*it).second.second.size();
  }
  return computeHPYPPredictive(cw, // cw
                               tw, // tw
                               payload.sumCustomers, // c
                               payload.sumTables, // t
                               parentProbability,
                               discount,
                               concentration);
}


IHPYPBaseRestaurant::TypeVector SimpleFullRestaurant::getTypeVector(
    void* payloadPtr) const {
  Payload& payload = *((Payload*)payloadPtr);
  IHPYPBaseRestaurant::TypeVector typeVector;
  typeVector.reserve(payload.tableMap.size());
  for (Payload::TableMap::iterator it = payload.tableMap.begin();
       it != payload.tableMap.end(); ++it) {
    typeVector.push_back(it->first); 
  }
  return typeVector;
}


const IPayloadFactory& SimpleFullRestaurant::getFactory() const {
  return this->payloadFactory;
}


void SimpleFullRestaurant::updateAfterSplit(void* longerPayloadPtr, 
                                            void* shorterPayloadPtr, 
                                            double discountBeforeSplit, 
                                            double discountAfterSplit,
                                            bool parentOnly) const {
  Payload& payload = *((Payload*)longerPayloadPtr);
  Payload& newParent = *((Payload*)shorterPayloadPtr);
    
  // make sure the parent is empty
  assert(newParent.sumCustomers == 0);
  assert(newParent.sumTables == 0);
  assert(newParent.tableMap.size() == 0);

  for(Payload::TableMap::iterator it = payload.tableMap.begin();
      it != payload.tableMap.end(); ++it) {
    e_type type = it->first;
    Payload::Arrangement& arrangement = payload.tableMap[type];
    Payload::Arrangement& parentArrangement = newParent.tableMap[type];

    if (arrangement.first == 1) { // just one customer -- can't split
      // seat customer at his own table in parent
      ++newParent.sumCustomers; // c
      ++newParent.sumTables; // t
      parentArrangement.first = 1; // cw
      parentArrangement.second.push_back(1); // cwk
    } else {
      // The correct thing to do is the following: 
      // For each type s in this restaurant
      //      for each table with cwk customers
      //          sample a partition cwkj from a PYP(-d1d2,d2)
      //          make the resulting table the tables in this restaurant
      //          for each k, create a table in the parent restaurant and 
      //          seat |cskj| customers on that table

      // make copy of old seating arrangement
      std::vector<l_type> oldTables = arrangement.second;

      // aliases
      std::vector<l_type>& newTables = arrangement.second;
      std::vector<l_type>& parentTables = parentArrangement.second;

      if (!parentOnly) {
        // remove old seating arrangement from this restaurant
        payload.sumCustomers -= arrangement.first;
        payload.sumTables -= oldTables.size();
        arrangement.first = 0;
        newTables.clear();
      }

      for (l_type k = 0; k < (l_type)oldTables.size(); ++k) { // for each table
        std::vector<int> frag = sample_crp_c(discountAfterSplit,
                                             -discountBeforeSplit,
                                             oldTables[k]);
        // add table to parent with cwk = frag.size()
        parentTables.push_back(frag.size());
        parentArrangement.first += frag.size();
        newParent.sumCustomers += frag.size();
        newParent.sumTables += 1;

        if (!parentOnly) {
          // add split tables to this restaurant
          for (int j = 0; j < (int)frag.size(); ++j) {
            int cwkj = frag[j];
            newTables.push_back(cwkj);
            payload.sumTables += 1;
            arrangement.first += cwkj;
            payload.sumCustomers += cwkj;
          }
        }
      }
    }
  }
}


bool SimpleFullRestaurant::addCustomer(void*  payloadPtr, 
                                       e_type type, 
                                       double parentProbability, 
                                       double discount, 
                                       double concentration,
                                       void*  additionalData) const {
  assert(additionalData == NULL); 

  Payload& payload = *((Payload*)payloadPtr);
  Payload::Arrangement& arrangement = payload.tableMap[type];
  std::vector<l_type>& tables = arrangement.second;

  ++payload.sumCustomers; // c
  ++arrangement.first; // cw 

  if (arrangement.first == 1) {
    // first customer sits at first table
    tables.push_back(1);
    ++payload.sumTables;
    return true;
  }
 
  // probs for old tables: \propto cwk - d
  d_vec tableProbs(tables.size() + 1, 0);
  for(int i = 0; i < (int)tables.size(); ++i) {
    tableProbs[i] = tables[i] - discount;
  }
  // prob for new table: \propto (alpha + d*t)*P0
  // this can be 0 for the first customer if concentration=0, but that is ok
  tableProbs[tables.size()] = 
    (concentration + discount*payload.sumTables)*parentProbability;
  
  // choose table for customer to sit at
  int table = sample_unnormalized_pdf(tableProbs);
  assert(table <= (int)tables.size());

  if(table == (int)tables.size()) {
    // sit at new table
    tables.push_back(1);
    ++payload.sumTables;
    return true;
  } else {
    // existing table
    ++tables[table];
    return false;
  }
}


bool SimpleFullRestaurant::removeCustomer(void* payloadPtr, 
                                          e_type type,
                                          double discount,
                                          void* additionalData) const {
  assert(additionalData == NULL);

  Payload& payload = *((Payload*)payloadPtr);
  assert(payload.tableMap.count(type) == 1);
  Payload::Arrangement& arrangement = payload.tableMap[type];
  std::vector<l_type>& tables = arrangement.second;

  assert(payload.sumCustomers > 0);
  assert(payload.sumTables > 0);
  assert(arrangement.first > 0);
  assert(tables.size() > 0);
  
  --payload.sumCustomers; // c
  --arrangement.first; // cw
  
  d_vec tableProbs(tables.begin(), tables.end()); // cast to double
  // chose a table to delete the customer from; prob proportional to table size
  int table = sample_unnormalized_pdf(tableProbs);
  assert(table < (int)tables.size());
  
  // remove customer from table
  --tables[table];
  if (tables[table] == 0) { // if table became empty
      tables.erase(tables.begin() + table); // drop from table list
      --payload.sumTables;
      return true;
  } else {
      return false;
  }
}


void* SimpleFullRestaurant::createAdditionalData(void* payloadPtr, 
                                                 double discount, 
                                                 double concentration) const {
  return NULL; // we don't need additional data
}


void SimpleFullRestaurant::freeAdditionalData(void* additionalData) const {
  // nothing to be done
}


std::string SimpleFullRestaurant::toString(void* payloadPtr) const {
  Payload& payload = *((Payload*)payloadPtr);
  std::ostringstream out;

  out << "[";
  for(Payload::TableMap::iterator it = payload.tableMap.begin();
      it != payload.tableMap.end(); ++it) {
    out << it->first << ":(" << it->second.first 
        << "/" << it->second.second.size() << "), ";
  }
  out << "]";
  return out.str();
}


bool SimpleFullRestaurant::checkConsistency(void* payloadPtr) const {
  Payload& payload = *((Payload*)payloadPtr);
  bool consistent = true; 

  int sumCustomers = 0;
  int sumTables = 0;

  for (Payload::TableMap::iterator it = payload.tableMap.begin();
       it != payload.tableMap.end(); ++it) {
    int cw = sum(it->second.second); 
    if (it->second.first != cw) {
      consistent = false;
      tracer << "sum_k(cwk) [" << cw << "] != cw [" << it->second.first << "]"
             << std::endl;
    }
    sumCustomers += cw;
    sumTables += it->second.second.size();
  }

  consistent =    (sumCustomers == payload.sumCustomers) 
               && (sumTables == payload.sumTables)
               && consistent;
  if (!consistent) {
    tracer << "Restaurant internally inconsistent!" 
           << " " << sumCustomers << "!=" << payload.sumCustomers
           << ", " << sumTables << "!=" << payload.sumTables 
           << std::endl;
  }
  return consistent;
}
    
void* SimpleFullRestaurant::newPayloadFromOther(
    const IHPYPBaseRestaurant& fromRestaurant,
    void* otherPayload,
    double discount) const {
  Payload* payload = new Payload();
  IHPYPBaseRestaurant::TypeVector types = 
      fromRestaurant.getTypeVector(otherPayload);
  
  for (IHPYPBaseRestaurant::TypeVectorIterator it = types.begin();
       it != types.end(); ++it) {
    e_type type = *it;
    int cw = fromRestaurant.getC(otherPayload, type);
    int tw = fromRestaurant.getT(otherPayload, type);
    payload->tableMap[type].first = cw;
    std::vector<int> cwk = sample_crp_ct(discount, cw, tw);
    assert((int)cwk.size() == tw);

    payload->tableMap[type].second = cwk;
    payload->sumCustomers += cw;
    payload->sumTables += tw;
  }
  
  return payload;
}


void SimpleFullRestaurant::Payload::serialize(InArchive & ar,
                                              const unsigned int version) {
  ar >> tableMap;
  ar >> sumCustomers;
  ar >> sumTables;
}


void SimpleFullRestaurant::Payload::serialize(OutArchive & ar,
                                              const unsigned int version) {
  ar << tableMap;
  ar << sumCustomers;
  ar << sumTables;
}


void SimpleFullRestaurant::PayloadFactory::save(
    void* payloadPtr, OutArchive& oa) const {
  oa << *((Payload*)payloadPtr);
}


void* SimpleFullRestaurant::PayloadFactory::load(InArchive& ia) const {
  Payload* p = new Payload();
  ia >> *p;
  return p;
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////   class HistogramRestaurant   ///////////////////////////
////////////////////////////////////////////////////////////////////////////////

l_type HistogramRestaurant::getC(void* payloadPtr, e_type type) const {
  Payload& payload = *((Payload*)payloadPtr);
  Payload::TableMap::iterator it = payload.tableMap.find(type);
  if (it != payload.tableMap.end()) {
    return it->second.cw; //cw
  } else {
    return 0;
  }
}


l_type HistogramRestaurant::getC(void* payloadPtr) const {
  return ((Payload*)payloadPtr)->sumCustomers;
}


l_type HistogramRestaurant::getT(void* payloadPtr, e_type type) const {
  Payload& payload = *((Payload*)payloadPtr);
  Payload::TableMap::iterator it = payload.tableMap.find(type);
  if (it != payload.tableMap.end()) {
    return it->second.tw; //tw
  } else {
    return 0;
  }
}


l_type HistogramRestaurant::getT(void* payloadPtr) const {
  return ((Payload*)payloadPtr)->sumTables;
}


double HistogramRestaurant::computeProbability(void*  payloadPtr,
                                               e_type type, 
                                               double parentProbability,
                                               double discount, 
                                               double concentration) const {
  Payload& payload = *((Payload*)payloadPtr);
  int cw = 0;
  int tw = 0;
  Payload::TableMap::iterator it = payload.tableMap.find(type);
  if (it != payload.tableMap.end()) {
    cw = (*it).second.cw;
    tw = (*it).second.tw;
  }
  return computeHPYPPredictive(cw, // cw
                               tw, // tw
                               payload.sumCustomers, // c
                               payload.sumTables, // t
                               parentProbability,
                               discount,
                               concentration);
}


IHPYPBaseRestaurant::TypeVector HistogramRestaurant::getTypeVector(
    void* payloadPtr) const {
  Payload& payload = *((Payload*)payloadPtr);
  IHPYPBaseRestaurant::TypeVector typeVector;
  typeVector.reserve(payload.tableMap.size());
  for (Payload::TableMap::iterator it = payload.tableMap.begin();
       it != payload.tableMap.end(); ++it) {
    typeVector.push_back(it->first); 
  }
  return typeVector;
}


const IPayloadFactory& HistogramRestaurant::getFactory() const {
  return this->payloadFactory;
}


void HistogramRestaurant::updateAfterSplit(void* longerPayloadPtr, 
                                           void* shorterPayloadPtr, 
                                           double discountBeforeSplit, 
                                           double discountAfterSplit,
                                           bool parentOnly) const {
  Payload& payload = *((Payload*)longerPayloadPtr);
  Payload& newParent = *((Payload*)shorterPayloadPtr);
    
  // make sure the parent is empty
  assert(newParent.sumCustomers == 0);
  assert(newParent.sumTables == 0);
  assert(newParent.tableMap.size() == 0);

  for(Payload::TableMap::iterator it = payload.tableMap.begin();
      it != payload.tableMap.end(); ++it) {
    e_type type = it->first;
    Payload::Arrangement& arrangement = payload.tableMap[type];
    Payload::Arrangement& parentArrangement = newParent.tableMap[type];

    if (arrangement.cw == 1) { // just one customer -- can't split
      // seat customer at his own table in parent
      ++newParent.sumCustomers; // c
      ++newParent.sumTables; // t
      parentArrangement.cw = 1; // cw
      parentArrangement.tw = 1; // cw
      parentArrangement.histogram[1] = 1; // cwk
    } else {
      // The correct thing to do is the following: 
      // For each type s in this restaurant
      //      for each table with cwk customers
      //          sample a partition cwkj from a PYP(-d1d2,d2)
      //          make the resulting table the tables in this restaurant
      //          for each k, create a table in the parent restaurant and 
      //          seat |cskj| customers on that table

      // make copy of old seating arrangement
      Payload::Arrangement oldArrangement = arrangement;

      if (!parentOnly) {
        // remove old seating arrangement from this restaurant
        payload.sumCustomers -= arrangement.cw;
        payload.sumTables -= arrangement.tw;
        arrangement.cw = 0;
        arrangement.tw = 0;
        arrangement.histogram.clear();
      }

      for (Payload::Histogram::iterator it = oldArrangement.histogram.begin();
           it != oldArrangement.histogram.end();
           ++it) { // for all table sizes
        for (l_type k = 0; k < (*it).second; ++k) { // all tables of this size
          std::vector<int> frag = sample_crp_c(discountAfterSplit,
              -discountBeforeSplit,
              (*it).first);
          // add table to parent with cwk = frag.size()
          parentArrangement.histogram[frag.size()] += 1;
          parentArrangement.cw += frag.size();
          newParent.sumCustomers += frag.size();
          parentArrangement.tw += 1;
          newParent.sumTables += 1;
          
          if (!parentOnly) {
            // add split tables to this restaurant
            for (int j = 0; j < (int)frag.size(); ++j) {
              int cwkj = frag[j];
              arrangement.histogram[cwkj] += 1;
              arrangement.cw += cwkj;
              arrangement.tw += 1;
              payload.sumTables += 1;
              payload.sumCustomers += cwkj;
            }
          }
        }
      }
    }
  }
}


bool HistogramRestaurant::addCustomer(void*  payloadPtr, 
                                      e_type type, 
                                      double parentProbability, 
                                      double discount, 
                                      double concentration,
                                      void* additionalData) const {
  assert(additionalData == NULL); 
  tracer << "HistogramRestaurant::addCustomer(" << type << "," 
         << parentProbability << "," << discount << "," << concentration 
         << ", " << additionalData
         << ")" << std::endl;

  Payload& payload = *((Payload*)payloadPtr);
  Payload::Arrangement& arrangement = payload.tableMap[type];

  payload.sumCustomers += 1; // c
  arrangement.cw += 1; // cw 
  if (arrangement.cw == 1) {
    // first customer sits at the first table
    // this special case is needed as otherwise things will break when alpha=0
    // and the restaurant has 0 customers in the singleton bucket
    arrangement.histogram[1] += 1;
    arrangement.tw += 1;
    payload.sumTables += 1;
    return true;
  }


  int numBuckets = arrangement.histogram.size(); 
  d_vec tableProbs(numBuckets + 1, 0);
  std::vector<int> assignment(numBuckets,0);
  int i = 0;
  for(Payload::Histogram::iterator it = arrangement.histogram.begin();
      it != arrangement.histogram.end();
      ++it) {
    // prob for joining a table of size k: \propto (k - d)t[k]
    tableProbs[i] = ((*it).first - discount) * (*it).second;
    assignment[i] = (*it).first;
    ++i;
  }
  // prob for new table: \propto (alpha + d*t)*P0
  // this can be 0 for the first customer if concentration=0, but that is ok
  tableProbs[numBuckets] = 
    (concentration + discount*payload.sumTables)*parentProbability;
  
  // choose table for customer to sit at
  int sample = sample_unnormalized_pdf(tableProbs);
  assert(sample <= (int)numBuckets);

  if(sample == (int)numBuckets) {
    // sit at new table
    arrangement.histogram[1] += 1;
    arrangement.tw += 1;
    payload.sumTables += 1;
    return true;
  } else {
    // existing table
    arrangement.histogram[assignment[sample]] -= 1;
    if (arrangement.histogram[assignment[sample]] == 0) {
      // delete empty bucket from histogram
      arrangement.histogram.erase(assignment[sample]);
    }
    arrangement.histogram[assignment[sample]+1] += 1;
    return false;
  }
}


bool HistogramRestaurant::removeCustomer(void* payloadPtr, 
                                         e_type type,
                                         double discount,
                                         void* additionalData) const {
  assert(additionalData == NULL); 
  //assert(this->checkConsistency(payloadPtr));

  Payload& payload = *((Payload*)payloadPtr);
  Payload::Arrangement& arrangement = payload.tableMap[type];

  arrangement.cw -= 1; // cw 
  payload.sumCustomers -= 1; // c

  assert(arrangement.cw >= 0);
  assert(payload.sumCustomers >= 0);

  int numBuckets = arrangement.histogram.size(); 
  int singletonBucket = -1; // invalid bucket
  d_vec tableProbs(numBuckets, 0);
  std::vector<int> assignment(numBuckets,0);

  int i = 0;
  for(Payload::Histogram::iterator it = arrangement.histogram.begin();
      it != arrangement.histogram.end();
      ++it) {
    // prob for choosing a bucket k*t[k]
    tableProbs[i] = (*it).first * (*it).second;
    assignment[i] = (*it).first;
    if ((*it).first == 1) {
      singletonBucket = i;
    }
    ++i;
  }
  
  // choose table for customer to sit at
  int sample = sample_unnormalized_pdf(tableProbs);
  assert(sample <= (int)numBuckets);
  assert(tableProbs[sample] > 0);

  if (sample == singletonBucket) {
    tracer << "HistogramRestaurant: deleting from singleton bucket" << std::endl;
    assert(arrangement.histogram[1] > 0);
    // singleton -> drop table
    arrangement.histogram[1] -= 1;
    if (arrangement.histogram[1] == 0) {
      // delete empty bucket from histogram
      //arrangement.histogram.erase(1);
    }
    arrangement.tw -= 1;
    payload.sumTables -= 1;

    assert(arrangement.tw >= 0);
    assert(payload.sumTables >= 0);
    
    return true;
  } else {
    // non-singleton bucket
    arrangement.histogram[assignment[sample]] -= 1;
    assert(arrangement.histogram[assignment[sample]] >= 0);
    if (arrangement.histogram[assignment[sample]] == 0) {
      // delete empty bucket from histogram
      arrangement.histogram.erase(assignment[sample]);
    }
    arrangement.histogram[assignment[sample]-1] += 1;
    return false;
  }
}


void* HistogramRestaurant::createAdditionalData(void* payloadPtr, 
                                                double discount, 
                                                double concentration) const {
  return NULL;
}


void HistogramRestaurant::freeAdditionalData(void* additionalData) const {
  // nothing to be done
}


std::string HistogramRestaurant::toString(void* payloadPtr) const {
  Payload& payload = *((Payload*)payloadPtr);
  std::ostringstream out;

  out << "[";
  for (Payload::TableMap::iterator it = payload.tableMap.begin();
      it != payload.tableMap.end(); ++it) {
    out << it->first << ":(" << it->second.cw 
        << "/" << it->second.tw << " {";
    for (Payload::Histogram::iterator hit = (*it).second.histogram.begin();
         hit != (*it).second.histogram.end();
         ++hit) {
      out << (*hit).first << ":" << (*hit).second << ", ";
    }
    out << "}, ";
  }
  out << "]";
  return out.str();
  
}


bool HistogramRestaurant::checkConsistency(void* payloadPtr) const {
  Payload& payload = *((Payload*)payloadPtr);
  bool consistent = true; 

  int sumCustomers = 0;
  int sumTables = 0;

  for (Payload::TableMap::iterator it = payload.tableMap.begin();
       it != payload.tableMap.end(); ++it) {
    int cw = 0;
    int tw = 0;
    for (Payload::Histogram::iterator hit = (*it).second.histogram.begin();
         hit != (*it).second.histogram.end();
         ++hit) {
      cw += (*hit).first * (*hit).second;
      tw += (*hit).second;
    }

    if (it->second.cw != cw || it->second.tw != tw) {
      consistent = false;
      std::cout << "sum_k(cwk) [" << cw << "] != cw [" << it->second.cw << "]"
                << std::endl;
    }
    sumCustomers += cw;
    sumTables += tw;
  }

  consistent =    (sumCustomers == payload.sumCustomers) 
               && (sumTables == payload.sumTables)
               && consistent;
  if (!consistent) {
    std::cout << "Restaurant internally inconsistent!" 
              << " " << sumCustomers << "!=" << payload.sumCustomers
              << ", " << sumTables << "!=" << payload.sumTables 
              << std::endl;
  }
  return consistent;
}

/**
 * Construct a HistogramRestaurant::Payload from any other type of
 * payload by resampling a full seating arrangement for each type.
 */
void* HistogramRestaurant::newPayloadFromOther(
    const IHPYPBaseRestaurant& fromRestaurant, 
    void* otherPayload, 
    double discount) const {
  Payload* payload = new Payload();
  IHPYPBaseRestaurant::TypeVector types = 
      fromRestaurant.getTypeVector(otherPayload);
  
  for (IHPYPBaseRestaurant::TypeVectorIterator it = types.begin();
       it != types.end(); ++it) {
    e_type type = *it;
    int cw = fromRestaurant.getC(otherPayload, type);
    int tw = fromRestaurant.getT(otherPayload, type);
    Payload::Arrangement& arrangement = payload->tableMap[type];
    arrangement.cw = cw;
    arrangement.tw = tw;
    std::vector<int> cwk = sample_crp_ct(discount, cw, tw);
    assert((int)cwk.size() == tw);
    for (int i = 0; i < (int)cwk.size(); ++i) {
      arrangement.histogram[cwk[i]] += 1;
    }

    payload->sumCustomers += cw;
    payload->sumTables += tw;
  }
  
  return payload;
}

void HistogramRestaurant::Payload::Arrangement::serialize(
    InArchive & ar, const unsigned int version) {
  ar >> cw;
  ar >> tw;
  ar >> histogram;
}


void HistogramRestaurant::Payload::Arrangement::serialize(
    OutArchive & ar, const unsigned int version) {
  ar << cw;
  ar << tw;
  ar << histogram;
}


void HistogramRestaurant::Payload::serialize(
    InArchive & ar, const unsigned int version) {
  ar >> tableMap;
  ar >> sumCustomers;
  ar >> sumTables;
}


void HistogramRestaurant::Payload::serialize(
    OutArchive & ar, const unsigned int version) {
  ar << tableMap;
  ar << sumCustomers;
  ar << sumTables;
}


void HistogramRestaurant::PayloadFactory::save(
    void* payloadPtr, OutArchive& oa) const {
  oa << *((Payload*)payloadPtr);
}


void* HistogramRestaurant::PayloadFactory::load(InArchive& ia) const {
  Payload* p = new Payload();
  ia >> *p;
  return p;
}


////////////////////////////////////////////////////////////////////////////////
//////////////////////   class BaseCompactRestaurant   /////////////////////////
////////////////////////////////////////////////////////////////////////////////

l_type BaseCompactRestaurant::getC(void* payloadPtr, e_type type) const {
  Payload& payload = *((Payload*)payloadPtr);
  Payload::TableMap::iterator it = payload.tableMap.find(type);
  if (it != payload.tableMap.end()) {
    return (*it).second.first; // cw
  } else {
    return 0;
  }
}


l_type BaseCompactRestaurant::getC(void* payloadPtr) const {
  return ((Payload*)payloadPtr)->sumCustomers;
}


l_type BaseCompactRestaurant::getT(void* payloadPtr, e_type type) const {
  Payload& payload = *((Payload*)payloadPtr);
  Payload::TableMap::iterator it = payload.tableMap.find(type);
  if (it != payload.tableMap.end()) {
    return (*it).second.second; // tw
  } else {
    return 0;
  }
}


l_type BaseCompactRestaurant::getT(void* payloadPtr) const {
  return ((Payload*)payloadPtr)->sumTables;
}


double BaseCompactRestaurant::computeProbability(void*  payloadPtr,
                                                 e_type type, 
                                                 double parentProbability,
                                                 double discount, 
                                                 double concentration) const {
  Payload& payload = *((Payload*)payloadPtr);
  int cw = 0;
  int tw = 0;
  Payload::TableMap::iterator it = payload.tableMap.find(type);
  if (it != payload.tableMap.end()) {
    cw = (*it).second.first;
    tw = (*it).second.second;
  }
  return computeHPYPPredictive(cw, // cw
                               tw, // tw
                               payload.sumCustomers, // c
                               payload.sumTables, // t
                               parentProbability,
                               discount,
                               concentration);
}


IHPYPBaseRestaurant::TypeVector BaseCompactRestaurant::getTypeVector(
    void* payloadPtr) const {
  Payload& payload = *((Payload*)payloadPtr);
  IHPYPBaseRestaurant::TypeVector typeVector;
  typeVector.reserve(payload.tableMap.size());
  for (Payload::TableMap::iterator it = payload.tableMap.begin();
       it != payload.tableMap.end(); ++it) {
    typeVector.push_back((*it).first); 
  }
  return typeVector;
}


const IPayloadFactory& BaseCompactRestaurant::getFactory() const {
  return this->payloadFactory;
}


void BaseCompactRestaurant::updateAfterSplit(void* longerPayloadPtr, 
                                             void* shorterPayloadPtr, 
                                             double discountBeforeSplit, 
                                             double discountAfterSplit,
                                             bool parentOnly) const {
  tracer << "BaseCompactPayload::updateAfterSplit("
         << this->toString(longerPayloadPtr)
         << ", " << this->toString(shorterPayloadPtr)
         << ", " << discountBeforeSplit
         << ", " << discountAfterSplit
         << ")"  << std::endl;
  
  Payload& payload = *((Payload*)longerPayloadPtr);
  Payload& newParent = *((Payload*)shorterPayloadPtr);
    
  // make sure the parent is empty
  assert(newParent.sumCustomers == 0);
  assert(newParent.sumTables == 0);
  assert(newParent.tableMap.size() == 0);

  for(Payload::TableMap::iterator it = payload.tableMap.begin();
      it != payload.tableMap.end(); ++it) {
    e_type type = (*it).first;
    Payload::Arrangement& arrangement = payload.tableMap[type];
    Payload::Arrangement& parentArrangement = newParent.tableMap[type];
  
    if (arrangement.first == 1) { // just one customer -- can't split
      // seat customer at his own table in parent
      ++newParent.sumCustomers; // c
      ++newParent.sumTables; // t
      parentArrangement.first = 1; // cw
      parentArrangement.second = 1; // tw
    } else {
      // in order to split the restaurant we will re-instantiate a full 
      // seating arrangement
      std::vector<int> cwk = sample_crp_ct(discountBeforeSplit, 
                                           arrangement.first,
                                           arrangement.second);
      
      int totalTables = 0;
      for (int k = 0; k < (int)cwk.size(); ++k) { // for each table
        std::vector<int> frag = sample_crp_c(discountAfterSplit,
                                             -discountBeforeSplit,cwk[k]);
        totalTables += frag.size();
      }

      // should have at least as many tables after split
      assert(totalTables >= arrangement.second);
      
      // parent has totalTables customers around tw tables
      parentArrangement.first = totalTables;
      parentArrangement.second = arrangement.second;
      newParent.sumCustomers += parentArrangement.first;
      newParent.sumTables    += parentArrangement.second;
      
      if (!parentOnly) {
        // bottom has unchanged # of customers at totalTables tables
        payload.sumTables -= arrangement.second;
        arrangement.second = totalTables;
        payload.sumTables += totalTables;
      }
    }
  }
}
  

bool BaseCompactRestaurant::addCustomer(void*  payloadPtr, 
                                        e_type type, 
                                        double parentProbability, 
                                        double discount, 
                                        double concentration, 
                                        void*  additionalData) const {
  tracer << "BaseCompactRestaurant::addCustomer(" << type << "," 
         << parentProbability << "," << discount << "," << concentration 
         << ", " << additionalData
         << ")" << std::endl;

  Payload& payload = *((Payload*)payloadPtr);
  Payload::Arrangement& arrangement = payload.tableMap[type];

  // add customer, incC returns incremented count
  l_type cw = arrangement.first;
  arrangement.first += 1; // inc(cw)
  payload.sumCustomers += 1; // inc(c)

  if (cw == 0) {
    // first customer always creates a table
    arrangement.second += 1; // inc(tw)
    payload.sumTables += 1; // inc(t)
    return true;
  } else {
    double incTProb =   (concentration + discount*payload.sumTables)
      * parentProbability;
    incTProb = incTProb/(incTProb + cw - arrangement.second * discount);
    if (coin(incTProb)) {
      arrangement.second += 1;
      payload.sumTables += 1;
      return true;
    } else {
      return false;
    }
  }
}


std::string BaseCompactRestaurant::toString(void* payloadPtr) const {
  Payload& payload = *((Payload*)payloadPtr);
  std::ostringstream out;

  out << "[";
  for(Payload::TableMap::iterator it = payload.tableMap.begin();
      it != payload.tableMap.end(); ++it) {
    out << (*it).first << ":(" << (*it).second.first 
        << "/" << (*it).second.second << "), ";
  }
  out << "]";
  return out.str();
}


bool BaseCompactRestaurant::checkConsistency(void* payloadPtr) const {
  Payload& payload = *((Payload*)payloadPtr);
  bool consistent = true; 

  int sumCustomers = 0;
  int sumTables = 0;

  for (Payload::TableMap::iterator it = payload.tableMap.begin();
       it != payload.tableMap.end(); ++it) {
    sumCustomers += (*it).second.first;
    sumTables    += (*it).second.second;
  }

  consistent =    (sumCustomers == payload.sumCustomers) 
               && (sumTables == payload.sumTables)
               && consistent;
  if (!consistent) {
    tracer << "Restaurant internally inconsistent!" 
           << " " << sumCustomers << "!=" << payload.sumCustomers
           << ", " << sumTables << "!=" << payload.sumTables 
           << std::endl;
  }
  return consistent;
}


void BaseCompactRestaurant::Payload::serialize(
    InArchive & ar, const unsigned int version) {
  ar >> tableMap;
  ar >> sumCustomers;
  ar >> sumTables;
}


void BaseCompactRestaurant::Payload::serialize(
    OutArchive & ar, const unsigned int version) {
  ar << tableMap;
  ar << sumCustomers;
  ar << sumTables;
}


void BaseCompactRestaurant::PayloadFactory::save(
    void* payloadPtr, OutArchive& oa) const {
  oa << *((Payload*)payloadPtr);
}


void* BaseCompactRestaurant::PayloadFactory::load(InArchive& ia) const {
  Payload* p = new Payload();
  ia >> *p;
  return p;
}
////////////////////////////////////////////////////////////////////////////////
//////////////////////   class ReinstantiatingCompactRestaurant   //////////////
////////////////////////////////////////////////////////////////////////////////

bool ReinstantiatingCompactRestaurant::removeCustomer(
    void* payloadPtr, e_type type, double discount,
    void* additionalData) const {
  Payload& payload = *((Payload*)payloadPtr);
  Payload::Arrangement& arrangement = payload.tableMap[type];

  bool removedTable;
  if (additionalData != NULL) {
    removedTable = this->fullRestaurant.removeCustomer(additionalData,
                                                       type,
                                                       discount,
                                                       NULL);
  } else {
    std::cerr << "Additional data MUST be provided for now!" << std::endl;
    exit(1);
  }
  arrangement.first -= 1;
  payload.sumCustomers -= 1;
  if (removedTable) {
    arrangement.second -= 1;
    payload.sumTables -= 1;
  }

  return removedTable;
}


void* ReinstantiatingCompactRestaurant::createAdditionalData(
    void* payloadPtr, double discount, double concentration) const {
  return this->fullRestaurant.newPayloadFromOther(*this, payloadPtr, discount);
}


void ReinstantiatingCompactRestaurant::freeAdditionalData(
    void* additionalData) const {
  this->fullRestaurant.getFactory().recycle(additionalData);
}


bool ReinstantiatingCompactRestaurant::addCustomer(
    void*  payloadPtr, e_type type, double parentProbability, double discount,
    double concentration, void*  additionalData) const {
  if (additionalData != NULL) {
    // need to stay in sync with full restaurant
    Payload& payload = *((Payload*)payloadPtr);
    Payload::Arrangement& arrangement = payload.tableMap[type];
    arrangement.first += 1; // inc(cw)
    payload.sumCustomers += 1; // inc(c)
    if (this->fullRestaurant.addCustomer(additionalData,
                                         type,
                                         parentProbability,
                                         discount, 
                                         concentration, 
                                         NULL)) {
    arrangement.second += 1; // inc(tw)
    payload.sumTables += 1; // inc(t)
    return true;
    } else {
      return false;
    }
  } else {
    return BaseCompactRestaurant::addCustomer(payloadPtr, 
                                              type,
                                              parentProbability, 
                                              discount,
                                              concentration,
                                              additionalData);
  }
}



////////////////////////////////////////////////////////////////////////////////
//////////////////////   class StirlingCompactRestaurant ///////////////////////
////////////////////////////////////////////////////////////////////////////////

bool StirlingCompactRestaurant::removeCustomer(void* payloadPtr, 
                                               e_type type,
                                               double discount,
                                               void* additionalData) const {
  Payload& payload = *((Payload*)payloadPtr);
  Payload::Arrangement& arrangement = payload.tableMap[type];
 
  boost::scoped_ptr<stirling_generator_full_log> additionalDataDeleter;
  if (additionalData == NULL) {
    // this will delete the additionalData at the end of the function
    additionalDataDeleter.reset(new stirling_generator_full_log(
        discount, arrangement.first, arrangement.second));
    additionalData = additionalDataDeleter.get();
  }

  double decTProb = ((stirling_generator_full_log*)additionalData)->ratio(
      arrangement.first, arrangement.second);

  arrangement.first -= 1;
  payload.sumCustomers -= 1;

  if (arrangement.first == 0 || coin(decTProb)) {
    arrangement.second -= 1;
    payload.sumTables -= 1;
    return true;
  } else {
    return false;
  }
}


void* StirlingCompactRestaurant::createAdditionalData(
    void* payloadPtr, double discount, double concentration) const {
    return new stirling_generator_full_log(discount, 1, 1);
}


void StirlingCompactRestaurant::freeAdditionalData(void* additionalData) const {
  delete (stirling_generator_full_log*)additionalData;

}



////////////////////////////////////////////////////////////////////////////////
//////////////////////   class KneserKeyRestaurant   ///////////////////////////
////////////////////////////////////////////////////////////////////////////////

l_type KneserNeyRestaurant::getC(void* payloadPtr, e_type type) const {
  Payload& payload = *((Payload*)payloadPtr);
  Payload::TableMap::iterator it = payload.tableMap.find(type);
  if (it != payload.tableMap.end()) {
    return (*it).second; // cw
  } else {
    return 0;
  }
}


l_type KneserNeyRestaurant::getC(void* payloadPtr) const {
  return ((Payload*)payloadPtr)->sumCustomers;
}


l_type KneserNeyRestaurant::getT(void* payloadPtr, e_type type) const {
  Payload& payload = *((Payload*)payloadPtr);
  Payload::TableMap::iterator it = payload.tableMap.find(type);
  if (it != payload.tableMap.end()) {
    return 1; // if we have a customer of this type, we have exactly one table
  } else {
    return 0;
  }
}


l_type KneserNeyRestaurant::getT(void* payloadPtr) const {
  return ((Payload*)payloadPtr)->tableMap.size();
}


double KneserNeyRestaurant::computeProbability(void*  payloadPtr,
                                               e_type type, 
                                               double parentProbability,
                                               double discount, 
                                               double concentration) const {
  Payload& payload = *((Payload*)payloadPtr);
  
  if (payload.sumCustomers == 0) {
    return parentProbability;
  }

  Payload::TableMap::iterator it = payload.tableMap.find(type);
  int cw = 0;
  int tw = 0;
  if (it != payload.tableMap.end()) {
    cw = (*it).second;
    tw = 1;
  }

  return computeHPYPPredictive(cw, // cw
                               tw, // tw
                               payload.sumCustomers, // c
                               payload.tableMap.size(), // t
                               parentProbability,
                               discount,
                               concentration);
}


IHPYPBaseRestaurant::TypeVector KneserNeyRestaurant::getTypeVector(
    void* payloadPtr) const {
  Payload& payload = *((Payload*)payloadPtr);
  IHPYPBaseRestaurant::TypeVector typeVector;
  typeVector.reserve(payload.tableMap.size());
  for (Payload::TableMap::iterator it = payload.tableMap.begin();
       it != payload.tableMap.end(); ++it) {
    typeVector.push_back((*it).first); 
  }
  return typeVector;
}


const IPayloadFactory& KneserNeyRestaurant::getFactory() const {
  return this->payloadFactory;
}


void KneserNeyRestaurant::updateAfterSplit(void* longerPayloadPtr, 
                                           void* shorterPayloadPtr, 
                                           double discountBeforeSplit, 
                                           double discountAfterSplit,
                                           bool parentOnly) const {
  tracer << "KneserNeyRestaurant::updateAfterSplit("
         << this->toString(longerPayloadPtr)
         << ", " << this->toString(shorterPayloadPtr)
         << ", " << discountBeforeSplit
         << ", " << discountAfterSplit
         << ")"  << std::endl;
  
  Payload& payload = *((Payload*)longerPayloadPtr);
  Payload& newParent = *((Payload*)shorterPayloadPtr);
    
  // make sure the parent is empty
  assert(newParent.sumCustomers == 0);
  assert(newParent.tableMap.size() == 0);

  for(Payload::TableMap::iterator it = payload.tableMap.begin();
      it != payload.tableMap.end(); ++it) {
    e_type type = (*it).first;
    newParent.tableMap[type] = 1;
  }
  newParent.sumCustomers = payload.tableMap.size();
}
  

bool KneserNeyRestaurant::addCustomer(void*  payloadPtr, 
                                        e_type type, 
                                        double parentProbability, 
                                        double discount, 
                                        double concentration, 
                                        void*  additionalData) const {
  tracer << "KneserNeyRestaurant::addCustomer(" << type << "," 
         << parentProbability << "," << discount << "," << concentration 
         << ", " << additionalData
         << ")" << std::endl;

  Payload& payload = *((Payload*)payloadPtr);
  int& cw = payload.tableMap[type];
  cw += 1;
  payload.sumCustomers += 1;
  return (cw == 1); // true if we created a new table
}


bool KneserNeyRestaurant::removeCustomer(
    void* payloadPtr, e_type type, double discount,
    void* additionalData) const {
  Payload& payload = *((Payload*)payloadPtr);
  Payload::TableMap::iterator it = payload.tableMap.find(type);
  int& cw = (*it).second;

  cw -= 1;
  payload.sumCustomers -= 1;
  if (cw == 0) {
    payload.tableMap.erase(it);
    return true;
  } else {
    return false;
  }
}


std::string KneserNeyRestaurant::toString(void* payloadPtr) const {
  Payload& payload = *((Payload*)payloadPtr);
  std::ostringstream out;

  out << "[";
  for(Payload::TableMap::iterator it = payload.tableMap.begin();
      it != payload.tableMap.end(); ++it) {
    out << (*it).first << ":" << (*it).second << ", ";
  }
  out << "]";
  return out.str();
}


bool KneserNeyRestaurant::checkConsistency(void* payloadPtr) const {
  Payload& payload = *((Payload*)payloadPtr);
  bool consistent = true; 

  int sumCustomers = 0;

  for (Payload::TableMap::iterator it = payload.tableMap.begin();
       it != payload.tableMap.end(); ++it) {
    sumCustomers += (*it).second;
  }

  consistent = sumCustomers == payload.sumCustomers;
  if (!consistent) {
    tracer << "Restaurant internally inconsistent!" 
           << " " << sumCustomers << "!=" << payload.sumCustomers
           << std::endl;
  }
  return consistent;
}


void KneserNeyRestaurant::Payload::serialize(
    InArchive & ar, const unsigned int version) {
  ar >> tableMap;
  ar >> sumCustomers;
}


void KneserNeyRestaurant::Payload::serialize(
    OutArchive & ar, const unsigned int version) {
  ar << tableMap;
  ar << sumCustomers;
}


void KneserNeyRestaurant::PayloadFactory::save(
    void* payloadPtr, OutArchive& oa) const {
  oa << *((Payload*)payloadPtr);
}


void* KneserNeyRestaurant::PayloadFactory::load(InArchive& ia) const {
  Payload* p = new Payload();
  ia >> *p;
  return p;
}


void* KneserNeyRestaurant::createAdditionalData(
    void* payloadPtr, double discount, double concentration) const {
  return NULL;
}


void KneserNeyRestaurant::freeAdditionalData(
    void* additionalData) const {
  // do nothing
}

    


double ExpectedTablesCompactRestaurant::computeProbability(void*  payloadPtr,
                                                 e_type type, 
                                                 double parentProbability,
                                                 double discount, 
                                                 double concentration) const {
  Payload& payload = *((Payload*)payloadPtr);
  int cw = 0;
  int tw = 0;
  Payload::TableMap::iterator it = payload.tableMap.find(type);
  if (it != payload.tableMap.end()) {
    cw = (*it).second.first;
    tw = (*it).second.second;
  }
  //double expectedNumberOfTables = PYPExpectedNumberOfTables(c);
  return computeHPYPPredictive(cw, // cw
                               tw * parentProbability,
                               payload.sumCustomers, // c
                               payload.sumTables, // t
                               parentProbability,
                               discount,
                               concentration);
}


double PowerLawRestaurant::computeProbability(void*  payloadPtr,
                                               e_type type, 
                                               double parentProbability,
                                               double discount, 
                                               double concentration) const {
  Payload& payload = *((Payload*)payloadPtr);
  
  if (payload.sumCustomers == 0) {
    return parentProbability;
  }

  Payload::TableMap::iterator it = payload.tableMap.find(type);
  int cw = 0;
  double tw = 0;
  double t = 0;
  const double POWER = 0.5;
  if (it != payload.tableMap.end()) {
    cw = (*it).second;
    tw = pow(cw, POWER); //TODO: should parameterize this
  }
  for (Payload::TableMap::iterator it = payload.tableMap.begin(); it != payload.tableMap.end(); ++it) {
    t += pow((*it).second, POWER);
  }

  return computeHPYPPredictiveDouble(cw, // cw
                                     tw, // tw
                                     payload.sumCustomers, // c
                                     t,
                                     parentProbability,
                                     discount,
                                     concentration);
}

bool FractionalRestaurant::addCustomer(void*  payloadPtr, 
                                       e_type type, 
                                       double parentProbability, 
                                       double discount, 
                                       double concentration, 
                                       void*  additionalData) const {
  tracer << "FractionalRestaurant::addCustomer(" << type << "," 
         << parentProbability << "," << discount << "," << concentration 
         << ", " << additionalData
         << ")" << std::endl;

  Payload& payload = *((Payload*)payloadPtr);
  int& cw = payload.tableMap[type];
  cw += 1;
  payload.sumCustomers += 1;
  return (cw == 1); // true if we created a new table
}

double FractionalRestaurant::computeProbability(void*  payloadPtr,
                                                e_type type, 
                                                double parentProbability,
                                                double discount, 
                                                double concentration) const {
  Payload& payload = *((Payload*)payloadPtr);
  
  if (payload.sumCustomers == 0) {
    return parentProbability;
  }

  Payload::TableMap::iterator it = payload.tableMap.find(type);
  l_type cw = 0;
  l_type tw = 0;
  if (it != payload.tableMap.end()) {
    cw = (*it).second;
    tw = 1;
  }

  return computeHPYPPredictive(cw, // cw
                               tw, // tw
                               payload.sumCustomers, // c
                               payload.tableMap.size(), // t
                               parentProbability,
                               discount,
                               concentration);
}

}} // namespace gatsby::libplump
