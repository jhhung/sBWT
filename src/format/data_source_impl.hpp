/**																																								  
 *  @file data_source_impl.hpp
 *  @brief provide interface class to select from file_path and curl data sources 
 *  @author C-Salt Corp.
 */
#ifndef DATA_SOURCE_IMPL_HPP_
#define DATA_SOURCE_IMPL_HPP_  
//#include "curl_recv_wrapper.hpp"

/**
 * @class DataSource
 * @brief provide generic form of data source
 * @tparam SOURCE_TYPE an enum served as an interface to indicate the data source is a curl_gz_format, curl_plain_format, or a ifstream_plain_format
 */
template < SOURCE_TYPE, typename... >
class DataSource 
{};

/**
 * @class DataSource < SOURCE_TYPE::IFSTREAM_TYPE >
 * @brief Specialized version of class DataSource exclusively for the situation that the data sources are in the format of ifstream
 */
template <>
class DataSource < SOURCE_TYPE::IFSTREAM_TYPE >
{
public:
	size_t source_num_;
	std::vector <uint64_t> file_size_;
	std::vector < std::ifstream* > file_handle;

/**
 * @brief constructor
 */
	DataSource ()
	{}

	~DataSource (void)
	{
		for ( auto& j : file_handle )
			delete j;
	}
/**
 * @brief constructor
 */
	DataSource ( std::vector <std::string> file_path, std::vector <uint64_t> sizein = std::vector<uint64_t>(0) )
		: source_num_ ( file_path.size () )
		, file_size_ (sizein)
	{
		for ( auto& i : file_path )
			file_handle.push_back ( new std::ifstream( i ) );																								
	}

/**
 * @brief main interface to get line from the index-th file_handle and return a std::string
 */
	std::string source_get_line (size_t index)
	{
		std::string temp;
		std::getline ( (*file_handle[index]), temp );
		return temp;
	}

/**
 * @brief main interface to get line from each of the file_handle and return a std::vector<std::string>
 */
	std::vector<std::string> source_get_line (void)
	{
		std::vector < std::string> return_vec;
		std::string temp;
		for ( auto i : file_handle )
		{
			std::getline ( (*i), temp );
			return_vec.push_back (temp);
		}
		return return_vec;
	}
};

#endif
