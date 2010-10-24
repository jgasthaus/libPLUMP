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
