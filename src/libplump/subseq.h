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

#ifndef SUBSEQ_H_
#define SUBSEQ_H_

#include <string>
#include <vector>
#include <sstream>
#include "config.h"

namespace gatsby { namespace libplump {

/*
 * Class that represents a subsequence of a sequence, i.e. its start and
 * end position within the sequence. Supplies the helper methods
 * getRelativeReversed to read a subsequence backwards, as well as 
 * suffixUntil that determines the longest common suffix between
 * two subsequences. Note: no reference to the underlying sequence is stored,
 * so it must be passed to all methods that require access to it.
 */
class SubSeq {
public:
	l_type start;
	l_type end;
    
    SubSeq() {
    	this->start = 0;
    	this->end = 0;
    }
    
    SubSeq(SubSeq& other) {
    	this->start = other.start;
    	this->end = other.end;
    }
    
    SubSeq(l_type start, l_type end) {
    	this->start = start;
    	this->end = end;
    }
    
    l_type length() {
	    return end - start;
    }

    e_type getRelativeReversed(seq_type* seq, l_type offset) {
	    return (*seq)[end - 1 - offset];
    }

    /*
     * Determine the first position where this subsequence and s differ, 
     * i.e. the length of the longest common suffix of this sequence and s.
     * The argument start sets an offset from where the comparison should be
     * started, e.g. start=1 makes the assumption that the longest common suffix
     * is at least of length one.
     */
    l_type suffixUntil(seq_type* seq, SubSeq& s, l_type start) {
    	l_type l = (length() < s.length()) ? length() : s.length();
    	//int l = this.length(); // assume this seq is always shorter
    	l_type i = start;
    	for (; i<l; i++) {
    		if (getRelativeReversed(seq,i) != s.getRelativeReversed(seq,i))
    			break;
    	}
    	return i;
    }

    static std::string toString(l_type start, l_type end, const seq_type& seq) {
    	std::ostringstream out;
    	for (l_type i=0;i<end-start;i++) {
    		out << seq[start + i];
    	}
    	return out.str();
    }

};


}} // namespace gatsby::libplump
#endif /* SUBSEQ_H_ */
