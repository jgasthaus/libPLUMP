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

#ifndef MINI_MAP_H_
#define MINI_MAP_H_

#include <algorithm>
#include <ostream>
#include <sstream>
#include <iterator>
#include <string>
#include <vector>
#include <boost/scoped_array.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>

namespace gatsby { namespace libplump {

/**
 * A map container that stores (key,value) pairs in a pair of 
 * sorted arrays, which are grown by doubling their size as needed.
 *
 * The benefits of this type of container over the usual 
 * red-black tree implementation of e.g. std::map are that:
 *   - small constant overhead (2 pointers to the arrays and one integer
 *     storing the current number of elements in the container.
 * 
 * Like for std::map, lookup is logarithmic, but insertion is linear.
 *
 * Only a subset of the operations required by the STL associative container
 * concept are supported, but the ones that are supported should behave in the
 * required way.
 */
template<class Key, class T, typename size_type = unsigned int>
class MiniMap {
  public:
    typedef Key key_type;
    typedef T mapped_type;
    typedef std::pair<const Key,T> value_type;

    class const_iterator; // defined below
    
    // alias const_iterator for compatability with STL map
    typedef const_iterator iterator;

    MiniMap() : keys(new Key[1]), values(new T[1]), _size(0) {}

    MiniMap(const MiniMap& other) : 
        keys(new Key[other.capacity()]),
        values(new T[other.capacity()]) {
          std::copy(other.keys.get(), 
                    other.keys.get() + other._size, 
                    this->keys.get());
          std::copy(other.values.get(), 
                    other.values.get() + other._size,
                    this->values.get());
          this->_size = other._size;
    }

    MiniMap& operator=(const MiniMap& other) {
      if (this != &other) {
        size_type cap = other.capacity();

        Key* newKeys = 0;
        T* newValues = 0;

        try {
          newKeys = new Key[cap];
          newValues = new T[cap];

          std::copy(other.keys.get(), other.keys.get() + other._size, newKeys);
          std::copy(other.values.get(), 
                    other.values.get() + other._size,
                    newValues);

        } 
        catch (...) {
          delete[] newKeys;
          delete[] newValues;
          throw;
        }

        this->keys.reset(newKeys);
        this->values.reset(newValues);
        _size = other._size;
      }
      return *this;
    }

    ~MiniMap() {}

    void clear() {
      this->keys.reset(new Key[1]);
      this->values.reset(new T[1]);
      this->_size = 0;
    }


    size_type count(const Key& x) const {
      if (this->find(x) != this->end()) {
        return 1;
      } else {
        return 0;
      }
    }

    bool empty() const {
      return _size == 0;
    }

    size_type size() const {
      return _size;
    }

    const_iterator find(const Key& key) const {
      Key* pos = std::lower_bound(keys.get(),&keys[_size],key);
      if (pos != &keys[_size] && *pos == key) {
        return const_iterator(keys.get(),values.get(),pos);
      } else {
        return this->end();
      }
    }


    T& operator[](const Key& key) {
      Key* pos = std::lower_bound(keys.get(),&keys[_size],key);
      size_type offset = pos - keys.get();
      if (pos != &keys[_size] && *pos == key) {
        return values[offset];
      } else {
        return insert(offset, key, T());
      }
    }

    const_iterator insert(const_iterator position, const value_type& x) {
      if (*(position.pos) == x.first) {
        this->values[position.pos - position.keys] = x.second;
        return position; // iterator unchanged
      } else {
        Key* pos;
        if (x.first > *(position.pos)) {
          pos = std::lower_bound(position.pos,&keys[_size],x.first);
        } else {
          pos = std::lower_bound(keys.get(),&keys[_size],x.first);
        }
        int offset = pos - keys.get();
        if (pos != &keys[_size] && *pos == x.first) {
          values[offset] = x.second;
        } else {
          insert(offset, x.first, x.second);
        }
        return const_iterator(keys.get(), values.get(), keys.get() + offset); 
      }
    }


    const_iterator begin() const {
      return const_iterator(keys.get(),values.get(),keys.get());
    }


    const_iterator end() const {
      return const_iterator(keys.get(),values.get(),&keys[_size]);
    }


    std::string toString() const {
      std::ostringstream os;
      os << "{";
      for (const_iterator i = begin(); i != end(); ++i) {
        os << ((*i).first);
        os << ":" <<  (*i).second;
        os << ", ";
      }
      os << "}";
      return os.str();
    }
    

    class const_iterator : public std::iterator<std::output_iterator_tag,value_type> {
      public:

        const_iterator() : keys(NULL), values(NULL), pos(NULL) {} 

        const_iterator(Key* keys, T* values, Key* pos) : keys(keys), values(values), pos(pos) {}

        value_type operator *() const {
          return value_type(*pos, values[pos - keys]);
        }


        const_iterator& operator ++(){
          ++pos;
          return *this;
        }

        const_iterator operator ++(int){
          ++pos;
          return *this;
        }

        bool operator ==(const const_iterator& other) {
          return pos == other.pos;
        }

        bool operator !=(const const_iterator& other) {
          return pos != other.pos;
        }
      
      private:
        Key* keys;
        T* values;
        Key* pos;
        friend class MiniMap;
    };
  
  private:
    boost::scoped_array<Key> keys;
    boost::scoped_array<T> values;
    size_type _size;

    T& insert(size_type offset, const Key& key, const T& value) {
      ensure_capacity();
      for(size_type i=_size;i>offset;i--) {
        keys[i] = keys[i-1];
        values[i] = values[i-1];
      }
      keys[offset] = key;
      values[offset] = value;
      _size++;
      return values[offset];
    }

    /**
     * Ensure that the allocated arrays can hold at least one additional
     * element.
     *
     * If the storage arrays need to be resized, the are re-allocating to twice
     * the size of the old ones.
     */
    void ensure_capacity() {
      size_type c = capacity();
      assert(c >= _size); // we should never be over capacity
      if(c == _size) {
        size_type new_capacity = c * 2;
        Key* newKeys = new Key[new_capacity];
        T* newValues = new T[new_capacity];
        std::copy(this->keys.get(), this->keys.get() + this->_size, newKeys);
        std::copy(this->values.get(), 
                  this->values.get() + this->_size,
                  newValues);
        this->keys.reset(newKeys);
        this->values.reset(newValues);
      }
    }

    /** Returns the capacity of this map, i.e. the number of elements that the
     * currently allocated arrays can hold.
     *
     * As we always double the size of the arrays when we resize, the capacity
     * is equal to the smallest power of 2 that is larger than the number of
     * elements in the map.
     */
    size_type capacity() const {
      size_type powof2 = 1;
      // Double powof2 until >= val
      while (powof2 < _size) powof2 <<= 1;
      return powof2;
    }
    
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const {
      // copy data into STL vectors to simplify serialization
      std::vector<Key> keyVec(this->keys.get(),
                              this->keys.get() + this->_size);
      std::vector<T> valueVec(this->values.get(),
                              this->values.get() + this->_size);
      ar << _size;
      ar << keyVec;
      ar << valueVec;
    }
    

    template<class Archive>
    void load(Archive & ar, const unsigned int version) {
      std::vector<Key> keyVec;
      std::vector<T> valueVec;
      ar >> this->_size;
      ar >> keyVec;
      ar >> valueVec;
      this->keys.reset(new Key[this->capacity()]);
      std::copy(keyVec.begin(), keyVec.end(), this->keys.get());
      this->values.reset(new T[this->capacity()]);
      std::copy(valueVec.begin(), valueVec.end(), this->values.get());
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()
};


template<class Key, class T>
std::ostream& operator<<(std::ostream& stream, MiniMap<Key,T> map) {
  stream << map.toString();
  return stream;
}

}} // namespace gatsby::libplump


#endif /* MINI_MAP_H_ */
