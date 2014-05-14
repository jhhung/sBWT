#ifndef ABIT_HPP_
#define ABIT_HPP_
#include <map>
#include <locale>
#include <algorithm>
#include <bitset>
#include <iostream>
#include <deque>
#include <string>
#include <set>
#include "boost/dynamic_bitset.hpp"
#include "boost/utility/binary.hpp"
#include "../constant_def.hpp"
#include <cctype>

template <class T = std::string>
class ABSequence
	: public T
{
public:
	ABSequence(T &sequence)
		:T (sequence)
	{
		std::cerr << "Seqence Size : " << this->size() << std::endl;
	}
	ABSequence()
		:T ()
	{}
	std::string& getContent()
	{
		return *this;
	}
};

#endif
