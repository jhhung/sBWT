#ifndef MKQ_SORT2_HPP_
#define MKQ_SORT2_HPP_
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
#include "bucket_sort.hpp"
#include "constant_def.hpp"

//#define INTTYPE uint64_t

//
// Record_rank
class Record_rank_disable
{
public:
	std::vector<std::pair<INTTYPE,INTTYPE>> &_iD_same_rank;
	Record_rank_disable(std::vector<std::pair<INTTYPE,INTTYPE>> &id_same_rank)
		:_iD_same_rank(id_same_rank)
	{}
	inline void record_rank(INTTYPE a, INTTYPE b)
	{
		return;
	}
};

class Record_rank_enable
{
public:
	std::vector<std::pair<INTTYPE,INTTYPE>> &_iD_same_rank;
	Record_rank_enable(std::vector<std::pair<INTTYPE,INTTYPE>> &id_same_rank)
		:_iD_same_rank(id_same_rank)
	{}
	inline void record_rank(INTTYPE a, INTTYPE b)
	{
		_iD_same_rank.push_back({a, b});
		return;
	}
};

//
// COMPARE with dcs or without
template <typename SEQTYPE>
class Compare_DCS
{
public:
	SEQTYPE &_seq;
	std::function<bool(INTTYPE,INTTYPE)> &_dcs_compare;
	INTTYPE _limit_depth;
	
	Compare_DCS(SEQTYPE &s, std::function<bool(INTTYPE,INTTYPE)> &dp, INTTYPE ld)
		:_seq(s), _dcs_compare(dp), _limit_depth(ld)
	{};
	
	inline long default_compare(INTTYPE a,INTTYPE b, INTTYPE depth){		
		if( depth >= _limit_depth)
		{
			if( _dcs_compare(a-depth, b-depth) )
				return -1;
			else 
				return +1;
		}
		return _seq[ a ] - _seq[ b ];
	};
	
	template<class Iter>
	inline bool final_sort(Iter begin, Iter end)
	{
		std::sort( begin, end, 
			[this](typename Iter::value_type a, typename Iter::value_type b){
				return this->_dcs_compare(a,b);
			}
		);
	}
	
};

template <typename SEQTYPE>
class Compare_default
{
public:
	SEQTYPE &_seq;
	std::function<bool(INTTYPE,INTTYPE)> &_dcs_compare;
	INTTYPE _limit_depth;
	
	Compare_default(SEQTYPE &s, std::function<bool(INTTYPE,INTTYPE)> &dp, INTTYPE ld)
		:_seq(s), _dcs_compare(dp), _limit_depth(ld)
	{};
	//sadfasdf
	inline long default_compare(INTTYPE a,INTTYPE b, INTTYPE depth){		
		if( depth >= _limit_depth)
		{
			return 0;
		}
		return _seq[ a ] - _seq[ b ];
	};
	
	template<class Iter>
	inline bool final_sort(Iter begin, Iter end)
	{
		return true;
	}
	
};

// SORT_SMALL_N
template <typename SEQTYPE, typename VECTORTYPE>
class Sort_small_n_enable
{
public:
	SEQTYPE &_seq;
	std::function<INTTYPE(typename VECTORTYPE::value_type&)> &_getTableV;
	std::function<bool(INTTYPE,INTTYPE)> &_dcs_compare;
	INTTYPE _limit_depth;
	
	Sort_small_n_enable(SEQTYPE &s, std::function<INTTYPE(typename VECTORTYPE::value_type&)> &gtv, INTTYPE ld, std::function<bool(INTTYPE,INTTYPE)> &dp)
		:_seq(s), _getTableV(gtv), _limit_depth(ld), _dcs_compare(dp)
	{}
	
	inline bool sort_small_n(typename VECTORTYPE::iterator begin,typename VECTORTYPE::iterator end, INTTYPE depth, INTTYPE n){
		if(n>6)	return false;
		
		std::sort(begin, end,
			[&]( typename VECTORTYPE::value_type A, typename VECTORTYPE::value_type B){
				INTTYPE a(_getTableV(A)), b(_getTableV(B));
				for(INTTYPE i(depth); i < _limit_depth; ++i)
				{
					//if( a+i >= _seq.size() )
					//	return true;
					//else if( b+i >= _seq.size() )
					//	return false;
					
					if( _seq[(a+i)] > _seq[(b+i)] )
						return false;
					else if( _seq[(a+i)] < _seq[(b+i)] )
						return true;
				}
				return _dcs_compare(a,b);
			}
		);
		return true;//
	};
	inline bool sort_small_n2(typename VECTORTYPE::iterator begin, typename VECTORTYPE::iterator end, INTTYPE depth, INTTYPE n){
		if(n>6)	return false;
		
		std::sort(begin, end,
			[&]( typename VECTORTYPE::value_type A, typename VECTORTYPE::value_type B){
				INTTYPE a(A), b(B);
				for(INTTYPE i(depth); i < _limit_depth; ++i)
				{
					if( _seq[(a+i)] > _seq[(b+i)] )
						return false;
					else if( _seq[(a+i)] < _seq[(b+i)] )
						return true;
				}
				return _dcs_compare(a,b);
			}
		);
		return true;
	};
	
	
};

template <typename SEQTYPE, typename VECTORTYPE>
class Sort_small_n_disable
{
public:
	SEQTYPE &_seq;
	std::function<INTTYPE(typename VECTORTYPE::value_type&)> &_getTableV;
	std::function<bool(INTTYPE,INTTYPE)> &_dcs_compare;
	INTTYPE _limit_depth;
	
	Sort_small_n_disable(SEQTYPE &s, std::function<INTTYPE(typename VECTORTYPE::value_type&)> &gtv, INTTYPE ld, std::function<bool(INTTYPE,INTTYPE)> &dp)
		:_seq(s), _getTableV(gtv), _limit_depth(ld), _dcs_compare(dp)
	{}
	inline bool sort_small_n(typename VECTORTYPE::iterator begin, typename VECTORTYPE::iterator end, INTTYPE depth, INTTYPE n){
		return false;
	};
	inline bool sort_small_n2(typename VECTORTYPE::iterator begin, typename VECTORTYPE::iterator end, INTTYPE depth, INTTYPE n){
		return false;
	};
};




template	<	typename SEQTYPE, 
						typename VECTORTYPE,
						typename RECORD_RANK,
						template	<	typename > class COMPARE,
						template	<	typename, typename > class SORT_SMALL_N,
						template	<	typename,
												typename, 
												typename, 
												template <typename> class,
												template <typename, typename> class,
												typename ...
											> class SORT_BIG_N
					>
class Multikey_quicksort
	:public RECORD_RANK,
					COMPARE <SEQTYPE>,
					SORT_SMALL_N <SEQTYPE, VECTORTYPE>
{
private:
	SEQTYPE &seq;
	VECTORTYPE &vec;
	std::function<bool(INTTYPE,INTTYPE)> dcs_compare;
	SORT_BIG_N<SEQTYPE,VECTORTYPE,RECORD_RANK,COMPARE,SORT_SMALL_N> bucket_sort;
	INTTYPE limit_depth, max_depth;

public:
	std::vector<std::pair<INTTYPE,INTTYPE>> &iD_same_rank, tmp_same_rank;
	std::function<INTTYPE(typename VECTORTYPE::value_type& )> getTableV;
	std::function<long(INTTYPE,INTTYPE,INTTYPE)> default_compare2;
	
	
	Multikey_quicksort(SEQTYPE &sequence, VECTORTYPE &d, INTTYPE limit)
		:	vec(d),
			seq(sequence),
			limit_depth(limit),
			RECORD_RANK(iD_same_rank),
			COMPARE <SEQTYPE> (sequence, dcs_compare, limit),
			SORT_SMALL_N<SEQTYPE, VECTORTYPE>(sequence, getTableV, limit, dcs_compare),
			bucket_sort(sequence, vec, limit, iD_same_rank),
			iD_same_rank(tmp_same_rank)
	{}
	
	Multikey_quicksort(SEQTYPE &sequence, VECTORTYPE &d, INTTYPE limit, std::vector<std::pair<INTTYPE,INTTYPE>> &id_same_rank)
		:	vec(d),
			seq(sequence),
			limit_depth(limit),
			RECORD_RANK(id_same_rank),
			COMPARE <SEQTYPE> (sequence, dcs_compare, limit),
			SORT_SMALL_N<SEQTYPE, VECTORTYPE>(sequence, getTableV, limit, dcs_compare),
			bucket_sort(sequence, vec, limit, id_same_rank),
			iD_same_rank(id_same_rank)
	{
		getTableV = [](typename VECTORTYPE::value_type a){
			return a;
		};
	}
	
	Multikey_quicksort(SEQTYPE &sequence, VECTORTYPE &d, INTTYPE limit, std::function<bool(INTTYPE,INTTYPE)> &compare_func)
		:	vec(d),
			seq(sequence),
			limit_depth(limit),
			RECORD_RANK(iD_same_rank),
			dcs_compare(compare_func),
			COMPARE <SEQTYPE> (sequence, dcs_compare, limit),
			SORT_SMALL_N <SEQTYPE, VECTORTYPE>(sequence, getTableV, limit, compare_func),
			bucket_sort(sequence, vec, limit, compare_func),
			iD_same_rank(tmp_same_rank)
	{
		getTableV = [](typename VECTORTYPE::value_type a){
			return a;
		};
	}
	
	
	
	inline void sort(INTTYPE start, INTTYPE size, INTTYPE depth)
	{
		mkq_sort(vec.begin()-start,size,depth);
		bucket_sort.release();
	}
	
	inline void mkq_sort(typename VECTORTYPE::iterator x, INTTYPE n, INTTYPE depth)
	{
			
			if (n <= 1)
				return;
			
			if(bucket_sort.sort(x-vec.begin(), n, depth) )
				return;
			
			if(depth == (limit_depth+1))
			{
				this->record_rank(x-vec.begin(), x-vec.begin()+n-1 );
				this->final_sort(x, x+n);
				return;
			}
				
			if(this->sort_small_n(x, x+n, depth, n))
				return;
			
			INTTYPE a, b, c, d, r;
			long r2;
			INTTYPE v_idx(0);
			tx = x;
			
			INTTYPE pa, pb, pc, pd, pl, pm, pn, t;
			pl = 0;
			pm = 0 + (n/2);
			pn = 0 + (n-1);
			if (n > 30) { // On big arrays, pseudomedian of 9
					d = (n/8);
					pl = med3(pl, pl+d, pl+2*d, depth);
					pm = med3(pm-d, pm, pm+d, depth);
					pn = med3(pn-2*d, pn-d, pn, depth);
			}
			pm = med3(pl, pm, pn, depth);
			mswap(pl, pm);
			
			
			v_idx = real_idx(0,depth);
			a = b = 1;
			c = d = n-1;
			
			while(true)
			{
					while (b <= c && ( r2 = this->default_compare( real_idx(b,depth) , v_idx, depth) ) <= 0) {
							if (r2 == 0) { 
								mswap(a, b); 
								a++; 
							}
							b++;
					}
					while (b <= c && ( r2 = this->default_compare( real_idx(c,depth), v_idx, depth) ) >= 0) {
							if (r2 == 0) {
								mswap(c, d); 
								d--; 
							}
							c--;
					}
					if (b > c) break;
					mswap(b, c);
					b++;
					c--;
			}
			r = std::min(a, b-a);		 
			vecmswap(0, b-r, r);
			
			r = std::min(d-c, n-d-1); 
			vecmswap(b, n-r, r);

			r = b-a;
			mkq_sort(x, r, depth);
			
			//if(depth == limit_depth && (a + n-d-1) > 1){
			//	this->record_rank(x+r-vec.begin(), x+r-vec.begin()+(a + n-d-1)-1 );
			//}
			mkq_sort(x + r, a + n-d-1, depth+1);
			
			r = d-c;
			mkq_sort(x + n-r, r, depth);
			
	}

private:
	typename VECTORTYPE::iterator tx;
	inline INTTYPE med3(INTTYPE a, INTTYPE b, INTTYPE c, INTTYPE depth)
	{	 
			int va, vb, vc;
			if ((va=i2c(a,depth)) == (vb=i2c(b,depth)))
					return a;
			if ((vc=i2c(c,depth)) == va || vc == vb)
					return c;			 
			return va < vb ?
						(vb < vc ? b : (va < vc ? c : a ) )
					: (vb > vc ? b : (va < vc ? a : c ) );
	}

	inline INTTYPE real_idx(INTTYPE i, INTTYPE depth)
	{
		return getTableV(vec[tx-vec.begin()+i])+depth;
	}
	inline char i2c(INTTYPE i, INTTYPE depth)
	{
		return seq[ getTableV(vec[tx-vec.begin()+i])+depth ];
	}
	inline void mswap(INTTYPE a, INTTYPE b)
	{
		std::swap(vec[tx-vec.begin()+a],vec[tx-vec.begin()+b]);
	}
	inline void vecmswap(INTTYPE i, INTTYPE j, INTTYPE n)
	{
		std::swap_ranges(tx+i, tx+i+n, tx+j);
	}
	
};


















#endif
