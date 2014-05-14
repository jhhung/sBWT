/// @file file_reader.hpp
/// @brief define a FileReader_impl interface, for achieving get_next_entry selectively via NORMAL, multi-thread, and MPI policies 
/// @author C-Salt Corp.
#ifndef FILE_READER_HPP_
#define FILE_READER_HPP_
#include <vector>
#include <map>
#include <iostream>
#include <boost/ref.hpp>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/thread.hpp>
#include <functional>
#include "constant_def.hpp"
#include "file_reader_impl.hpp"
#include "tuple_utility.hpp"
#include "wrapper_tuple_utility.hpp"
#include "thread_pool_update.hpp"

/// @brief generic form of FileReader < ParallelTypes T, template <typename> class FORMAT, typename TUPLETYPE > is specialized into: \n 1) a normal version, i.e. a serial operation version; 2) a multi-thread version; and 3) a MPI version by assiging the corresponding enum value of ParallelTypes::NORMAL, ParallelTypes::M_T, and ParallelTypes::M_P_I. \n The parallel versions, i.e. the multi-thread and the MPI versions are respectively implemented with our thread_pool.hpp and  mpi_pool.hpp.
/// @tparam none-type parameter of int / enum of ParallelTypes, indicating what policies that the FileReader is going to employ
/// @tparam template template parameter of FORMAT, designed to be any one of current data formats
/// @tparam TUPLETYPE indicating data format of the got piece of data
template < ParallelTypes ParaType, template < typename> class FORMAT, typename TUPLETYPE, SOURCE_TYPE STYPE, typename... ARGS >
class FileReader 
	: public FileReader_impl < FORMAT, TUPLETYPE, STYPE, ARGS... >//< T, FORMAT, TUPLETYPE >
{};

/// @brief specialized form of the FileReader, with the ParallelType specialized to be NORMAL, indicating non-parallel scheme
template < template < typename > class FORMAT, typename TUPLETYPE, SOURCE_TYPE STYPE, typename... ARGS >
class FileReader < ParallelTypes::NORMAL, FORMAT, TUPLETYPE, STYPE, ARGS... >
	: public FileReader_impl < FORMAT, TUPLETYPE, STYPE, ARGS...  >
{
private:
/// @memberof FileReader
/// @brief Getting entries stored in the file_handles object of data_source object
/// @param std::vector<FORMAT<TUPLETYPE>> for holding the got Num entries, from the corresponding file_handles
/// @param size_t indicating the number of entries to be got
	void Read_impl (size_t index, std::vector< FORMAT <TUPLETYPE> > & a, size_t Num)
	{
		std::vector< FORMAT<TUPLETYPE> > temp;
		for (size_t i=0; i!=Num; ++i)
		{
			auto r = FileReader_impl<FORMAT, TUPLETYPE, STYPE, ARGS... >::get_next_entry (index);
			if (r.eof_flag)
				break;
			temp.push_back (r);
		}
			a.swap (temp);
	}

	void Read_impl (size_t index, std::vector< FORMAT <TUPLETYPE> > & a, size_t Num, size_t overloading_for_N_skip)
	{
		std::vector< FORMAT<TUPLETYPE> > temp;
		for (size_t i=0; i!=Num; ++i)
		{
			auto r = FileReader_impl<FORMAT, TUPLETYPE, STYPE, ARGS... >::get_next_entry (index);
			if (r.eof_flag)
				break;
			if (r.getSeq().find ('N') != std::string::npos || r.getSeq().find ('n') != std::string::npos )
				continue;
			else
				temp.push_back (r);	
		}
			a.swap (temp);
	}

protected:
	std::map < int, std::vector< FORMAT<TUPLETYPE> > >* result2;
//	int Job_count;

public:
/// @brief default constructor
	FileReader ( std::vector <std::string>& uni_file_in, 
				 std::map <int, std::vector< FORMAT<TUPLETYPE> > >* content,
				 std::vector <uint64_t> sizein = std::vector<uint64_t>(0) )  
		: FileReader_impl < FORMAT, TUPLETYPE, STYPE, ARGS...  > (uni_file_in, sizein) 
		, result2 ( content )
	{}
	
	FileReader ( std::vector <std::string>& uni_file_in,
				 std::vector <uint64_t> sizein = std::vector<uint64_t>(0) )  
		: FileReader_impl < FORMAT, TUPLETYPE, STYPE, ARGS...  > (uni_file_in, sizein) 
	{}

/// @brief getting entries stored in the files corresponding to the url/file_path
/// @param Num, defaulted to be 1, indicating how many entries the Read function is going to get in a single run
/// @return a pointer pointing to a map <int, vector< FORMAT<TUPLETYPE> > >, so as to return access result of the Read function
	std::map< int, std::vector< FORMAT <TUPLETYPE> > >* Read ( size_t Num = 1 )
	{	
		size_t Job_count=0;
		for ( auto itr = 0; itr != this->source_num_; ++itr )
		{
			Read_impl (itr, (*result2)[Job_count], Num);
			++Job_count;
		}
		return result2;
	}


	bool Read (std::map< int, std::vector< FORMAT <TUPLETYPE> > >* result, size_t Num = 1)
	{	
		size_t Job_count=0;
		for ( auto itr = 0; itr != this->source_num_; ++itr )
		{
			Read_impl (itr, (*result)[Job_count], Num);
			++Job_count;
		}
		bool flag = false;
		std::for_each ( this -> file_handle.begin(), this -> file_handle.end(), [&] ( decltype (*(this->file_handle.begin())) Q)//std::stringstream* Q )
        	{   if ( ! Q->good() )  flag = true;   });
		return flag;		
	}

	bool Read_NSkip (std::map< int, std::vector< FORMAT <TUPLETYPE> > >* result, size_t Num = 1)
	{	
		size_t Job_count=0;
		for ( auto itr = 0; itr != this->source_num_; ++itr )
		{
			Read_impl (itr, (*result)[Job_count], Num, 5566);
			++Job_count;
		}
		bool flag = false;
		std::for_each ( this -> file_handle.begin(), this -> file_handle.end(), [&] ( decltype (*(this->file_handle.begin())) Q)//std::stringstream* Q )
        	{   if ( ! Q->good() )  flag = true;   });
		return flag;		
	}


/// @brief getting entries stored in the files corresponding to the url/file_path
/// @tparam FUNCTOR, indicating the type of the to be executed functor 
/// @param Q, a functor for achieving certain operation on the accessed entries
/// @param Num, defaulted to be 1, indicating how many entries the Read function is going to get in a single run
/// @return a pointer pointing to a map <int, vector< FORMAT<TUPLETYPE> > >, so as to return access result of the Read function
	template <typename FUNCTOR>
	std::map< int, std::vector <FORMAT <TUPLETYPE> > >*  Read_combo ( FUNCTOR Q, size_t Num = 1 )
	{
		size_t Job_count=0;
		for ( auto itr = 0; itr != this->source_num_; ++itr )
		{
			Read_impl (itr, (*result2)[Job_count], Num);
			++Job_count;
		}
		Q();
		return result2;
	}
};

#endif


