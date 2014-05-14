/**
 * @file fasta_reader_impl.hpp
 * @brief Specialized version of class FileReader_impl for Fasta, for providing operation detail for getting a piece of FastA data from an ifstream object
 * This class provides a function to read the fasta format file.
 *
 * @author JHH corp.
 */
#ifndef FASTA_READER_IMPL_HPP_
#define FASTA_READER_IMPL_HPP_
#include <tuple>
#include "fasta.hpp"
#include "../file_reader_impl.hpp"
#include "boost/algorithm/string/trim.hpp"
#include "boost/algorithm/string/erase.hpp"
/**
 * @brief Read Fasta format data, including member functions: io_end (), get_next_entry (), and get_all_entry ()
 * @tparam Fasta template template parameter Fasta indicating current FileReader_impl is a specilization form of the base FileReader_impl template class for Fasta format
 * @tparam TUPLETYPE indicating Fasta data format, defaulted as tuple < string, string >
 * 
 */
template<class TUPLETYPE, SOURCE_TYPE STYPE> //Specialization for FASTA
class FileReader_impl < Fasta, TUPLETYPE, STYPE >
	: public DataSource < STYPE >
{
public:
typedef Fasta<TUPLETYPE> format_type;
typedef TUPLETYPE tuple_type;
typedef SOURCE_TYPE source_type; 
};

template< class TUPLETYPE>//ifstream Specialization for FASTA
class FileReader_impl <Fasta, TUPLETYPE, SOURCE_TYPE::IFSTREAM_TYPE>
    : public DataSource < SOURCE_TYPE::IFSTREAM_TYPE >
{
public:
	size_t file_num_;
///@typedef io_iterator in format of formatIoIterator < Fasta <TUPLETYPE>, FileReader_impl <Fasta, TUPLETYPE> > 

	typedef FormatIoIterator<
							Fasta<TUPLETYPE>, 
							FileReader_impl< Fasta, TUPLETYPE, SOURCE_TYPE::IFSTREAM_TYPE > 
							> io_iterator;
// @memberof FileReader_impl<Fasta, TUPLETYPE>
/// @brief return a pointer pointing to a Fasta <TUPLETYPE> with its eof_flag setted as true
/// @return shared_ptr < Fasta <TUPLETYPE> > (new Fasta <TUPLETYPE> (true) )
	io_iterator io_end ()
	{
		return io_iterator (std::shared_ptr < Fasta <TUPLETYPE> > (new Fasta <TUPLETYPE> (true) ) );
	}

    FileReader_impl ()
        : DataSource < SOURCE_TYPE::IFSTREAM_TYPE > ()
        , file_num_ (0)
    {}

    FileReader_impl (std::vector<std::string>& file_path, std::vector <uint64_t> sizein = std::vector<uint64_t>(0) )
        :  DataSource < SOURCE_TYPE::IFSTREAM_TYPE > (file_path, sizein)
        , file_num_ ( file_path.size() )
    {}

/// @memberof FileReader_impl<Fasta, TUPLETYPE>
/// @brief return an object of Fasta <TUPLETYPE> by means of reading an std::ifstream& object; \n
/// @param file_handle an std::ifstream& object 
/// @return Fasta <TUPLETYPE> 
	Fasta<TUPLETYPE> get_next_entry (size_t index)//(std::ifstream& file_handle)
	{
		std::string data, line, head;
		int is_head=1; //a flag indicating whether is the head line.
		Fasta<TUPLETYPE> result;
//		file_handle.peek();
		this->file_handle [index] -> peek();
		if ( this->file_handle[index]->eof() || this->file_handle[index]->fail() || this->file_handle[index]->bad() )
		{
        // The if (file_handle.fail() || file_handle.eof() || file_handle.bed() ) function detects ifstream fail/eof/bad flags, and 
        // return Fasta <TUPLETYPE> with enabled eof_flag when any one of the aformentioned flags indicates true. 
        // Whenever the returned Fasta <TUPLETYPE> has an enabled eof_flag, it is indicated that at least one of fail/eof/bad/ situations occurs
			result.eof_flag = true;
			return result;
		}
		while ( std::getline ( *(this->file_handle[index]), line ) )
		{
			boost::trim (line);
			if( line.front()  == '>' ) //the head line
			{
				if( is_head != 0 ) //confirm it's the first tuple
				{
//					head = line.substr(1); //get the head without '>'
std::get<0>(result.data) = line.substr(1);
					is_head = 0;
				}
				else
				{
					this->file_handle[index]->seekg(static_cast<long long int>(this->file_handle[index]->tellg()) - static_cast<long long int>(line.length()) - 1); 
					//move the pointer to the right position to the end of this tuple.
					break;
				}
			
			}
			else
			{	
				//data += line; //sequence data
std::get<1>(result.data)+=line;
			}
		}
//		boost::erase_all(data, " ");
//		boost::erase_all(data, "\t");
//		result.data = std::make_tuple ( std::move(head), std::move(data) );

		return result;
	}

/// @memberof FileReader_impl<Fasta, TUPLETYPE>
/// @brief receive a std::string of file_name, and accordingly read all Fasta <TUPLETYPE> entries thereof and update the call by reference parameter of data_cluster to record the got Fasta <TUPLETYPE> data entries
/// @param file_name in format of std::string 
/// @param data_cluster passed into get_all_entry function by reference, so as to keep the got Fasta <TUPLETYPE> data entries
//	static void get_all_entry (std::string& file_name, std::vector< Fasta<TUPLETYPE> >& data_cluster )
//	{
//		std::ifstream file_handle (file_name);
//		while (true)
//		{
//			auto x = get_next_entry (file_handle);
//			if ( x.eof_flag)
//				break;
//			data_cluster.push_back (x);
//		}
//	}
};

#endif
