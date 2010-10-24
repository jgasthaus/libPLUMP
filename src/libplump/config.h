/**
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
