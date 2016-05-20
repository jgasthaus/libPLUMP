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
