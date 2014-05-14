/// @file Searcher.hpp
#ifndef SEARCHER_HPP_
#define SEARCHER_HPP_
#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <ctime>
#include <atomic>

#include "boost/serialization/vector.hpp"
#include "boost/serialization/utility.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"
#include "boost/unordered_map.hpp"

#include "../compression/abit.hpp"
#include "../compression/jbit.hpp"
#include "../constant_def.hpp"
#include "../format/sam.hpp"

#include "job_distributer.hpp"
		
//template < Aligner_types AlignerType, ParallelTypes ParallelType, class AlignerTableType>
template< Aligner_types AlignerType
				 ,ParallelTypes ParallelType
				 ,class AlignerTableType
				 //,class QueryParserType
				 ,class FormatType
				 ,typename SearchReturnType
				 ,class... ARGS
				>
class Searcher_impl
{};



//template <class AlignerTableType>
//class Searcher_impl< AlignerTableType, Aligner_types::BWT_Aligner>
template< ParallelTypes ParallelType
				 ,class AlignerTableType
				 //,class QueryParserType
				 ,class FormatType
				 ,typename SearchReturnType
				 ,class... ARGS
				>
class Searcher_impl< Aligner_types::BWT_Aligner
							 			,ParallelType
							 			,AlignerTableType
							 			//,QueryParserType
							 			,FormatType
							 			,SearchReturnType
							 			,ARGS...
							 		 >

{
protected:
	AlignerTableType &abwt_table_;
	const INTTYPE &c_functions_interval;
	
	std::array<INTTYPE,256> mtable_;
	
	Searcher_impl(AlignerTableType &table, const INTTYPE &c_f_interval)
		: abwt_table_ (table)
		, c_functions_interval(c_f_interval)
	{
		mtable_['A'] = 0;
		mtable_['C'] = 1;
		mtable_['G'] = 2;
		mtable_['T'] = 3;
	}
	
	inline void init_exact_match( std::pair<INTTYPE, INTTYPE> &start_end_pos_, char c )
	{
		start_end_pos_.first = abwt_table_.c_function[ c ];
		start_end_pos_.second = abwt_table_.c_function[ c+1 ];
	}
	inline void init_exact_match( std::pair<INTTYPE, INTTYPE> &start_end_pos_, std::string tmp_str8 )
	{
		INTTYPE tmp_cs(0);
		for(INTTYPE i(0); i < c_functions_interval; ++i)
		{
			//tmp_cs += mtable_[ tmp_str8[ i ] ] << ((11-i)<<1); //slower.....@@
			//tmp_cs += mtable_[ tmp_str8[ i ] ] << ((c_functions_interval-1-i)<<1); //slower.....@@
			tmp_cs += mtable_[ tmp_str8[ i ] ] * std::pow(4,(c_functions_interval-1-i)) ;
			//std::cerr << tmp_str8[ i ] << " * " << std::pow(mtable_.size(),(7-i)) << " = " << tmp_c8 << std::endl;
		}
		//std::cerr << "tmp_cs " << tmp_cs << " size " << abwt_table_.c_functions.size() << std::endl;
		start_end_pos_.first = abwt_table_.c_functions[ tmp_cs ];
		start_end_pos_.second = abwt_table_.c_functions[ tmp_cs+1 ];
	}
	inline void init_exact_match( std::pair<INTTYPE, INTTYPE> &start_end_pos_, INTTYPE tmp_cs )
	{
		start_end_pos_.first = abwt_table_.c_functions[ tmp_cs ];
		start_end_pos_.second = abwt_table_.c_functions[ tmp_cs+1 ];
	}
	inline void exec_exact_match( std::pair<INTTYPE, INTTYPE> &start_end_pos_, char c )
	{
		start_end_pos_.first = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( start_end_pos_.first, c );
		start_end_pos_.second = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( start_end_pos_.second, c );
		//INTTYPE tmp_a = start_end_pos_.first;
		//INTTYPE tmp_b = start_end_pos_.second;
		
		//INTTYPE a = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( start_end_pos_.first, c );
		//INTTYPE b = abwt_table_.c_function[ c ] + abwt_table_.get_occ_using_jbwt( start_end_pos_.second, c );
		
		//start_end_pos_.first = abwt_table_.c_function[ c ] + abwt_table_.get_occ( start_end_pos_.first, c );
		//start_end_pos_.second = abwt_table_.c_function[ c ] + abwt_table_.get_occ( start_end_pos_.second, c );

		
	}
	
	inline INTTYPE find_nearest_mark (std::pair<INTTYPE, INTTYPE> &start_end_pos_, INTTYPE trace_point)
	{
		
		//INTTYPE tmp = trace_point;
		INTTYPE traceback_count(0);
		INTTYPE flocation = abwt_table_.fbwt[trace_point];
		
		for(INTTYPE i(0); i<flocation;++i)
		{
			//trace_point = abwt_table_.back_tracking(trace_point);
			trace_point = abwt_table_.back_tracking_using_jbwt(trace_point);
			++traceback_count;
		}
		
		auto pp = std::lower_bound( 
			abwt_table_.location_table.begin(), 
			abwt_table_.location_table.end(), 
			std::pair<INTTYPE, INTTYPE>{trace_point, 0}, 
				[]( const std::pair<INTTYPE, INTTYPE>& a, const std::pair<INTTYPE, INTTYPE>& b)
				{
					return a.first < b.first;
				} 
		);
		traceback_count += pp->second;
		return traceback_count;
		
		while(1)//trace to the nearest upstream location_table
		{
			
			//std::cerr << "last trace_point : " << trace_point << " flocation : " << flocation << std::endl;
			
			auto pp = std::lower_bound( 
				abwt_table_.location_table.begin(), 
				abwt_table_.location_table.end(), 
				std::pair<INTTYPE, INTTYPE>{trace_point, 0}, 
					[]( const std::pair<INTTYPE, INTTYPE>& a, const std::pair<INTTYPE, INTTYPE>& b)
					{
						return a.first < b.first;
					} 
			);
			if ( pp->first == trace_point )
			{
				//std::cerr << "traceback_count : " << traceback_count << std::endl;

				//traceback_count += hit.first->second;
				traceback_count += pp->second;
				break;
			}
			else//no hit yet
			{
				//trace_point = abwt_table_.back_tracking(trace_point);
				trace_point = abwt_table_.back_tracking_using_jbwt(trace_point);
				traceback_count++;

			}

		}
		return traceback_count;
	}
		
	inline void push_result(std::vector<INTTYPE>& result, INTTYPE pos_f, INTTYPE pos_e)
	{
		if (pos_e - pos_f > 100)
			return;
		for (int i = pos_f; i < pos_e; i++)
		{
			result.push_back( find_nearest_mark(i) );
		}
	}
	
public:
	static std::mutex mux_;
	void load_table(std::string prefix_name)
	{
	std::lock_guard <std::mutex> lk (mux_);
		if(abwt_table_.is_table_loaded == 0)
		{
			abwt_table_.is_table_loaded=1;
			if (prefix_name.back () != '.') {
				prefix_name += '.';
			}
			abwt_table_.readTable(prefix_name + "t_table.bwt");
			abwt_table_.readNPosLen (prefix_name + "NposLen.z");
			abwt_table_.readChrStartPos (prefix_name + "chrStart");
			abwt_table_.readChrLen (prefix_name + "chrLen");
		}
	}	
};

template< 
	ParallelTypes ParallelType
 ,class AlignerTableType
 //,class QueryParserType
 ,class FormatType
 ,typename SearchReturnType
 ,class... ARGS
>

std::mutex 
Searcher_impl< 
	Aligner_types::BWT_Aligner
	,ParallelType
	,AlignerTableType
	//,QueryParserType
	,FormatType
	,SearchReturnType
	,ARGS...
>::mux_;

//search_policy
//template<class AlignerTableType, int AlignerType, int SearcherPolicy ,typename... ARGS>
template< Aligner_types AlignerType
				 ,Searcher_types SearcherPolicy
				 ,ParallelTypes ParallelType
				 ,class AlignerTableType
				 //,class QueryParserType
				 //,class FormatType
				 ,typename SearchReturnType
				 ,typename... ARGS
				>
class Searcher
{};

/**
 * @struct 
 * 
 * @brief for sbwt search
 *  
 * @tparam
 *  
 */
//template <class AlignerTableType, typename... ARGS >
//class Searcher <AlignerTableType, Aligner_types::BWT_Aligner, Searcher_types::Exact_match, ARGS...>

template< 
	ParallelTypes ParallelType
	,class AlignerTableType
	//,class QueryParserType
	,class FormatType
	,typename SearchReturnType
	,class... ARGS
>
class Searcher< 
	Aligner_types::BWT_Aligner
	,Searcher_types::SBWT_exact_match
	,ParallelType
	,AlignerTableType
	//,QueryParserType
	,FormatType
	,SearchReturnType
	,ARGS...
	>
	:public Searcher_impl< 
		Aligner_types::BWT_Aligner
		,ParallelType
		,AlignerTableType
		//,QueryParserType
		,FormatType
		,SearchReturnType
		,ARGS...
	 >
{
private:
	const INTTYPE c_functions_interval;
	mutable std::pair<INTTYPE, INTTYPE> start_end_pos_;
public:
	typedef Searcher_impl<
		Aligner_types::BWT_Aligner
		,ParallelType
		,AlignerTableType
		//,QueryParserType
		,FormatType
		,SearchReturnType
		,ARGS...
	> SearchImplType;
	 
	std::array<char, 4> all_char;
	std::tuple<std::atomic<INTTYPE>, std::atomic<INTTYPE>, std::atomic<INTTYPE>, std::atomic<INTTYPE>, std::atomic<INTTYPE>> reads_map_count_;
							 				 
	Searcher(AlignerTableType& table)
		: SearchImplType (table, c_functions_interval)
		, c_functions_interval(12)
		, all_char{ {'A', 'C', 'G', 'T'} }
	{}
	
	typedef Sam< 
		std::tuple <
			std::string, //QNAME
			int, //SAM_FLAG, //FLAG
			std::string, //RNAME
			uint64_t, //POS
			int, //MAPQ
			std::string, //CIGAR
			std::string, //RNEXT
			uint64_t, //PNEXT
			int64_t, //TLEN
			std::string, //SEQ
			std::string, //QUAL
			int //NH_TAG
		>
	> SAM_FORMAT_TYPE;
	
	/// @brief 定義 Input type 單筆
	typedef std::map<int, FormatType > IN_DATA_TYPE_1;						 				 
	
	/// @brief 定義 Input type 多筆
	typedef std::map<int, std::vector<FormatType> > IN_DATA_TYPE_2;
	
	/// @brief 定義 Input type as file_list
	typedef std::vector<std::string> IN_DATA_TYPE_3;
	
	/// @brief 定義 Output type
	typedef std::vector< SAM_FORMAT_TYPE > OUT_DATA_TYPE_1;
	
	/// @brief 定義 Output type
	typedef std::stringstream OUT_DATA_TYPE_2;
	
	//for fast multiple cpus
	SearchReturnType
	search(IN_DATA_TYPE_3 &in_data, std::ostream &out, int limitNumber = 0, int strand = 2)
	{
		std::vector <uint64_t> fastq_size_vec(0);
		
		typedef std::tuple <std::string, std::string, std::string, std::string > TUPLETYPE;
		typedef FileReader < ParallelTypes::NORMAL, Fastq, TUPLETYPE, SOURCE_TYPE::IFSTREAM_TYPE> QueryParserType;
		
		QueryParserType fastq_parser(in_data, fastq_size_vec);
		
		Job_distributer <ParallelTypes::M_T, QueryParserType> jd;
		jd.distribute_jobs(fastq_parser, out, 0,
			[this, limitNumber, strand](FormatType &format_data, OUT_DATA_TYPE_2 &out_buffer)
			{
				OUT_DATA_TYPE_1 out_sam_vector;
				
				INTTYPE NH_tag = this->start_sbwt_match(format_data, out_sam_vector, limitNumber, strand);
				std::get<4>(reads_map_count_) += out_sam_vector.size();
				
				for(auto& sam : out_sam_vector)
				{
					out_buffer << sam;
				}
					
				std::get<0>(reads_map_count_) ++;
				if(NH_tag == 1)
				{
					std::get<1>(reads_map_count_) ++;
					std::get<2>(reads_map_count_) ++;
				}
				else if(NH_tag > 1)
				{
					std::get<1>(reads_map_count_) ++;
					std::get<3>(reads_map_count_) ++;
				}
			}
		);
	}
	
	struct SBWT_VAR
	{
		int maybe_result_number;
		int strand;
		std::vector< SAM_FORMAT_TYPE > &sam_result;
		std::string &query;
		std::string &_query;
		std::vector< std::tuple<std::string, INTTYPE, bool> > &result;
		FormatType &fq;
	};
	
	inline int start_sbwt_match( FormatType &fq, std::vector< SAM_FORMAT_TYPE >& sam_result, int limitNumber, int strand)
	{
		std::pair<INTTYPE,INTTYPE> start_end_pos_;
		std::string _query = fq.getSeq();
		/* reverse complement query string */
		std::string query {_query.crbegin(), _query.crend()};
		for (auto & c : query) 
		{
			switch (c) {
				case 'A': c = 'T'; break;
				case 'T': c = 'A'; break;
				case 'C': c = 'G'; break;
				case 'G': c = 'C'; break;
				default : throw "illegal char";
			}
		} /* end of RC */
		
		//init_exact_match( query.back() );
		this->init_exact_match( start_end_pos_, query.substr(query.size()-c_functions_interval,c_functions_interval) );
		
		for (int i=query.length()-1-c_functions_interval; i>=0 && start_end_pos_.first < start_end_pos_.second; i--)
		//for (int i=query.length()-1-1; i>=0 && start_end_pos_.first < start_end_pos_.second; i--)
		{
			this->exec_exact_match(start_end_pos_, query[i]);	
		}
		
		int maybe_result_number = start_end_pos_.second - start_end_pos_.first;
		
		//once there is any hit, this loop will run and print every hits
		if (start_end_pos_.second - start_end_pos_.first > limitNumber || start_end_pos_.second - start_end_pos_.first <= 0)
			return 0;


		std::vector< std::tuple<std::string, INTTYPE, bool> > result;
		result.reserve(maybe_result_number);
		
		SBWT_VAR VARS {maybe_result_number, strand, sam_result, query, _query, result, fq};

		for(INTTYPE i = start_end_pos_.first; i < start_end_pos_.second; ++i)
		{
			if(this->abwt_table_.fbwt[i] != 0)
				continue;
			auto pp = std::lower_bound( 
				this->abwt_table_.location_table.begin(),
				this->abwt_table_.location_table.end(), 
				std::pair<INTTYPE, INTTYPE>{ i, 0}, 
					[]( const std::pair<INTTYPE, INTTYPE>& a, const std::pair<INTTYPE, INTTYPE>& b)
					{
						return a.first < b.first;
					} 
			);
			if(pp->first == i)
			{
				push_position(VARS, pp->second);
				if(start_end_pos_.second - start_end_pos_.first == 1)
					break;
			}
		}
		//std::cout << fq.getName () << "\t" << _query << "\t" << maybe_result_number << "\t" << "CC" << std::endl;
		
		for(int cn(0); cn < all_char.size(); ++cn)
		{		
			this->find_possible(VARS, 1, all_char[cn], start_end_pos_.first, start_end_pos_.second);
		}
		
		
		//convert to sam
		if(result.size()==0)
			return 0;

		return position2sam(VARS);
		
	}
		
	inline void find_possible(SBWT_VAR &VARS, int len, char c, INTTYPE pos_f, INTTYPE pos_e)
	{
		if(len >= this->abwt_table_.interval)
			return;
		if (VARS.result.size() == VARS.maybe_result_number)
			return;
		
		INTTYPE n_pos_f = this->abwt_table_.c_function[ c ] + this->abwt_table_.get_occ_using_jbwt( pos_f, c );
		INTTYPE n_pos_e = this->abwt_table_.c_function[ c ] + this->abwt_table_.get_occ_using_jbwt( pos_e, c );
		
		//INTTYPE n_pos_f = this->abwt_table_.c_function[ c ] + this->abwt_table_.get_occ( pos_f, c );
		//INTTYPE n_pos_e = this->abwt_table_.c_function[ c ] + this->abwt_table_.get_occ( pos_e, c );
		
		if(n_pos_f >= n_pos_e) // no result
		{
			// no result
			return;
		}

		for(INTTYPE i = n_pos_f; i < n_pos_e; ++i)
		{
			if(this->abwt_table_.fbwt[i] != 0)
				continue;
			auto pp = std::lower_bound( 
				this->abwt_table_.location_table.begin(), 
				this->abwt_table_.location_table.end(), 
				std::pair<INTTYPE, INTTYPE>{ i, 0}, 
					[]( const std::pair<INTTYPE, INTTYPE>& a, const std::pair<INTTYPE, INTTYPE>& b)
					{
						return a.first < b.first;
					} 
			);
			if(pp->first == i)
			{
				push_position(VARS, len + pp->second);
				if(n_pos_e - n_pos_f == 1 )
					return;
			}
		}
		for(int cn(0); cn < all_char.size(); ++cn)
		{
			this->find_possible(VARS, len+1, all_char[cn], n_pos_f, n_pos_e);
		}
	}
	inline void push_position(SBWT_VAR &VARS, INTTYPE position)
	{
		bool isRC;
		if (position >= this->abwt_table_._realSize)
		{
			isRC = true; //0
			position = this->abwt_table_._realSize*2 - position - VARS._query.size();
		}
		else
		{
			isRC = false; //16
		}
		if( (VARS.strand==0 && !isRC) || (VARS.strand==1 && isRC) )
			return;//continue;
		
		auto lowerIter = this->abwt_table_.chr_start_pos.upper_bound (position);
		std::advance (lowerIter, -1);
		auto chr = lowerIter->second;

		auto lowerIter3 = this->abwt_table_.chr_start_pos.upper_bound (position + VARS._query.size() - 1);
		std::advance (lowerIter3, -1);
		auto chr3 = lowerIter3->second;
		if (chr != chr3) return;//continue;

		auto NLowerIter = this->abwt_table_.chr_umbiguous_starting_length.upper_bound (position);
		std::advance (NLowerIter, -1);

		auto NLowerIter3 = this->abwt_table_.chr_umbiguous_starting_length.upper_bound (position + VARS._query.size() - 1);
		std::advance (NLowerIter3, -1);
		if (NLowerIter != NLowerIter3) return;//continue;
		
		position = position - lowerIter->first + NLowerIter->second;
		
		VARS.result.push_back(std::make_tuple(chr, position, isRC) );

	}
	inline int position2sam(SBWT_VAR &VARS)
	{
		
		auto NH_tag = VARS.result.size();
		
		for (auto &result_position : VARS.result)
		{
			auto &chr = std::get<0>(result_position);
			auto &position = std::get<1>(result_position);
			auto &isRC = std::get<2>(result_position);
			
			if (!isRC)
			{	
				VARS.sam_result.emplace_back
				(
					std::move(
						std::make_tuple(
							std::move(VARS.fq.getName ()),
							SAM_FLAG::REVERSE_COMPLEMENTED, //16
							std::move (chr),
							position+1,
							255,
							std::to_string (VARS._query.size ()) + 'M',
							"*",
							0,
							0,
							VARS.query,
							std::move(VARS.fq.getRevQuality()),
							NH_tag
						)
					)
				);
			}
			else
			{
				VARS.sam_result.emplace_back
				(
					std::move(
						std::make_tuple(
							std::move(VARS.fq.getName ()),
							SAM_FLAG::MAPPED,
							std::move (chr),
							position+1,
							255,
							std::to_string (VARS._query.size ()) + 'M',
							"*",
							0,
							0,
							VARS._query,
							std::move(VARS.fq.getQuality ()),
							NH_tag
						)
					)
				);
			}//if
		
		}//for
		return NH_tag;
	}
	
};



#endif
