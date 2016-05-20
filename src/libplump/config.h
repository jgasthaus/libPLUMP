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

#ifndef CONFIG_H_
#define CONFIG_H_

#include <vector>
#include <iostream> // for cerr
#include <map>
#include <stdint.h>

#include "libplump/mini_map.h"

namespace gatsby { namespace libplump {

typedef int32_t e_type;
typedef int32_t l_type;
typedef std::vector<e_type> seq_type;

/*
 * We want an easy way to switch the map type that is used throughout to 
 * test the different performance characteristics.
 *
 * This should only be typedef'ed to a type that is compatible to the STL map
 * container.
 */
template<typename K, typename V>
struct MapType {
    typedef MiniMap<K,V> Type;
    //typedef std::map<K,V> Type;
};


}} // namespace gatsby::libplump

#endif /* CONFIG_H_ */
