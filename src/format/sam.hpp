/// @file sam.hpp
/// @brief Provide definition for Sam Format
#ifndef SAM_HPP_
#define SAM_HPP_
#include <tuple>
#include <string>
#include <cstdint>
#include <vector>
#include <bitset>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/mpl/string.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/tuple/tuple.hpp> 
#include <boost/preprocessor/repetition/enum_params.hpp> 
#include <boost/mpl/copy.hpp> 
#include <boost/algorithm/string/detail/trim.hpp>

#include "../constant_def.hpp"
#include "../tuple_utility.hpp"
#include "is_tuple_type.hpp"

/**
 * @struct Sam
 * @brief This is a Sam struct, including member elements:
 *  static type (enum format_type);
 *  data (TUPLETYPE);
 *  eof_flag (bool)
 * This struct can store the Sam format data, provide constructor for move and copy. We also overload the operator<< for print the fasta's data.\n
 */ 
/// @tparam TUPLETYPE defaulted as tuple < string, int, string, uint32_t, int, string, string, uint32_t, int, string >,indicate type of member element: data, \n
template < 
	class TUPLETYPE = std::tuple 
	<
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
		int // NH
	>
>
struct Sam
{
	static_assert( IsTupleType<TUPLETYPE>::value == true, "ARGUMENT IS NOT A TUPLE");
	TUPLETYPE data;//_f, data_b;

	static const format_types type = format_types::SAM;
	bool eof_flag;
	typedef TUPLETYPE TupleType;
	//default constructor
	Sam ()	//empty construct - can't pass unit_test
		: data ( TUPLETYPE() ), eof_flag (false) //_f ( TUPLETYPE() ), data_b ( TUPLETYPE() ), eof_flag (false)
	{}
	//Copy constructor
	Sam (const Sam<>& Sm)	
		: data (Sm.data), eof_flag (Sm.eof_flag)
	{}
	//Assignment operator
	Sam<>& operator=( const Sam<>& other)
	{
//		static_assert( is_tuple_type<decltype(other.data)>::value == true, "assignment's argument is not a tuple");
		eof_flag = other.eof_flag;
		data = other.data;
		return *this;
	}
	Sam (TUPLETYPE &data_in)
		: data (data_in), eof_flag (false) 
	{}
	Sam (TUPLETYPE &&data_in) 
		: data (std::move (data_in)), eof_flag (false) 
	{}
	//Move constructor
	Sam (Sam<> && other) 
		: data (std::move(other.data)), eof_flag (std::move(other.eof_flag))
	{}
	//Move assignment operator
	Sam<>& operator= (Sam<> && other)
	{
		data = std::move(other.data);
		eof_flag = std::move (other.eof_flag);
		return *this;
	}
	Sam (bool EofFlag) 
		: data ( TUPLETYPE() ),  eof_flag (EofFlag) 
	{}
	//String (copy) constructor
	Sam(std::string &line)
		:data( TUPLETYPE() ), eof_flag (false)
	{
		line2sam(line, data);
	}
	//String (move) constructor
	Sam(std::string &&line)
		:data( TUPLETYPE() ), eof_flag (false)
	{
		line2sam(line, data);
	}
	
	void line2sam(std::string & line, TUPLETYPE & tuple)
	{
		std::vector<std::string> tmp_splits;
		boost::split( tmp_splits, line, boost::is_any_of( "\t" ));

		tuple = 
		std::move(
			std::make_tuple
			(
				boost::lexical_cast< std::string > (tmp_splits[0]), //QNAME
	            boost::lexical_cast< int > (tmp_splits[1]), //SAM_FLAG, //FLAG
	            boost::lexical_cast< std::string > (tmp_splits[2]), //RNAME
	            boost::lexical_cast< uint64_t > (tmp_splits[3]), //POS
	            boost::lexical_cast< int > (tmp_splits[4]), //MAPQ
	            boost::lexical_cast< std::string > (tmp_splits[5]), //CIGAR
	            boost::lexical_cast< std::string > (tmp_splits[6]), //RNEXT
	            boost::lexical_cast< uint64_t > (tmp_splits[7]), //PNEXT
	            boost::lexical_cast< int64_t > (tmp_splits[8]), //TLEN
	            boost::lexical_cast< std::string > (tmp_splits[9]), //SEQ
	            boost::lexical_cast< std::string > (tmp_splits[10]), //QUAL
	            boost::lexical_cast< int > (tmp_splits[11]) //NH
			)
		);
	}
	
	friend std::ostream& operator<< (std::ostream& out, Sam& s)
	{
				
		std::string result;
		TupleUtility< TupleType , 12 >//std::tuple_size< TupleType >::value >
			::PrintTuple(s.data, result, 1);	//calling the overload version for tab delimited data structure
		result+='\n';
		out << result;
		return out;
	}

	std::string str()
	{
		std::string result;
		TupleUtility< TupleType , 12 >//std::tuple_size< TupleType >::value >
			::PrintTuple(data, result, 1);	//calling the overload version for tab delimited data structure
		return result;
	}

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		TupleUtility < TUPLETYPE, std::tuple_size<TUPLETYPE>::value >
			:: SerializeTuple (data, ar, version);
		ar & eof_flag;
	}
};
#endif
