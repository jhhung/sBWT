#ifndef JOB_DISTRIBUTER_HPP_
#define JOB_DISTRIBUTER_HPP_

#include <vector>
#include <fstream>
#include <iostream>

#include "../constant_def.hpp"
#include "../thread_pool_update.hpp"

template<ParallelTypes ParallelType, class QueryParserType>
class Job_distributer
{};

template<class QueryParserType>
class Job_distributer<ParallelTypes::NORMAL, QueryParserType>
{
public:
	Job_distributer();
};

class VectorParser
{};

template<class QueryParserType>
class Job_distributer<ParallelTypes::M_T, QueryParserType>
{
private:
	std::mutex reader_mutex, writer_mutex;
	INTTYPE in_group_reads_number;
public:
	Job_distributer()
		:in_group_reads_number(1000)
	{};
	
	
	void distribute_jobs(
			QueryParserType &file_parser
		, std::ostream &out
		, std::size_t nthreads
		, std::function<void(typename QueryParserType::format_type &, std::stringstream &)> functor
		)
	{
//		GlobalPool.ChangePoolSize(nthreads);

		bool eof_flag(false);
		int file_idx(0);
		while( !eof_flag )
		{
			std::vector< typename QueryParserType::format_type > *group = new std::vector< typename QueryParserType::format_type >();
			group->reserve( in_group_reads_number );
			
			for(INTTYPE i(0); i<in_group_reads_number; ++i)
			{
				group->emplace_back( std::move (file_parser.get_next_entry (file_idx) ) );
				if ( group->back().eof_flag )
				{
					//remove last nothing fastq
					group->pop_back();
					++file_idx;
					if(file_parser.file_num_ == file_idx)
					{
						eof_flag = true;
						break;
					}
				}
			}
			if(group->size() == 0)
				break;
			GlobalPool.JobPost(
				[this, group, &out, &functor]()
				{
					decltype(group) GROUP(group);
					std::stringstream out_buffer;
					for ( auto &fp_v : *GROUP)
					{
						functor(fp_v, out_buffer);
					}
					{
						std::lock_guard<std::mutex> lock(this->writer_mutex);
						out << out_buffer.rdbuf();
						out_buffer.str("");
					}
					delete GROUP;
				}, std::vector<size_t>(0) 
			);
		}
		GlobalPool.FlushPool();
		
	}
};
#endif
