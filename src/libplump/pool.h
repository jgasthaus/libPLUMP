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
 **/

#ifndef POOL_H_
#define POOL_H_

#include <boost/pool/poolfwd.hpp>
#include <boost/pool/pool.hpp>

namespace gatsby { namespace libplump {

template <class T>
  class PoolObject {
    public:
      static void* operator new(size_t size) {
        return memPool.malloc();
      }

      static void operator delete(void *p) {
        memPool.free(p);
      }

    private:
      static boost::pool<> memPool;
  };

template <class T>
boost::pool<> PoolObject<T>::memPool(sizeof(T));

}} // namespace gatsby::libplump

#endif
