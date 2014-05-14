#ifndef DIFFERENCE_COVER_HPP_
#define DIFFERENCE_COVER_HPP_
#include <fstream>
#include <map>
#include <list>
#include <random>
#include <algorithm>
#include <functional>
#include <boost/ref.hpp>
#include <utility>
#include <cstdlib>
#include <queue>
#include <tuple>
#include <ctime>

#include "compression/abit.hpp"
#include "mkq_sort.hpp"
#include "split_sort.hpp"
#include "lssort.hpp"


//#define INTTYPE uint64_t

template <typename T>
void FreeAll( T & t ) {
		T tmp;
		t.swap( tmp );
}

template<typename SEQTYPE, typename SORTTYPE>
class PSEUDO_DCS_FOR_SBWT
{
public:
	PSEUDO_DCS_FOR_SBWT (SEQTYPE &sequence, INTTYPE dcs_size=512)
	{
		std::cerr << "No DCS..." << std::endl;
	}
	inline bool compare(INTTYPE i, INTTYPE j)	{return true;}
	inline void release()
	{}
};

#endif
