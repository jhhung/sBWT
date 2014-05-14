/**
 * @file fasta.hpp
 * @brief Provide definition for Fasta Format. 
 *
 * @author JHH corp.
 */
#ifndef FASTA_HPP_
#define FASTA_HPP_
#include <string>
#include <tuple>
#include "../constant_def.hpp"
#include "../tuple_utility.hpp"
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include "is_tuple_type.hpp"
/**
 * @struct Fasta
 * @brief This is a Fasta struct, including member elements:
 *	static type (enum format_type);
 *	data (TUPLETYPE);
 *	eof_flag (bool)
 * This struct can store the fasta format data, provide constructor for move and copy. We also overload the operator<< for print the fasta's data.\n
 */

/**
 * @tparam TUPLETYPE defaulted as tuple< string, string  >, indicate type of member element: data \n
 * defaulted as tuple of: \n
 * 1.string: indicating a single-line description\n 
 * 2.string: indicating lines of sequence data
 */	

template <class TUPLETYPE = std::tuple<std::string, std::string> >
struct Fasta
{
    static_assert( IsTupleType<TUPLETYPE>::value == true, "ARGUMENT IS NOT A TUPLE");
	TUPLETYPE data;
	///@brief An enumerated static identifier, indicating Fasta<TUPLETYPE> format struct
	static const format_types type = format_types::FASTA;
	///@brief A flag indicates whether the Fasta data got by fasta_reader's get_next_entry function reaches the end-of-file.
	bool eof_flag;  //indicating whether file_handle reaches EOF
	///@typedef TupleType as a alias for TUPLETYPE
	typedef TUPLETYPE TupleType;
	
	///@brief Default constructor
	Fasta () 
		: data(), eof_flag( false )   //default constructor
	{}
	///@brief Copy constructor
	Fasta (const Fasta& in)
		: data( in.data ), eof_flag( in.eof_flag ) //Copy Constructor
	{}
	///@brief Assignment operator
	///@param empty
	Fasta& operator=(const Fasta& in)
	{
		data = in.data;
		eof_flag = in.eof_flag;
		//std::cerr << "Assignment operator" << std::endl;
		return *this;
	}
	/**
	 * @brief construct a Fasta <TUPLETYPE> with data_in_move in TUPLETYPE format
	 * @param data_in_move
	 */
	//Move Constructor for specific TUPLETYPE
	Fasta (TUPLETYPE && data_in_move) 
		: data(std::move (data_in_move) ), eof_flag (false)
	{}
	///@brief construct a Fasta <TUPLETYPE> with data_in in TUPLETYPE format
	///@param data_in in TUPLETYPE format
	Fasta (TUPLETYPE & data_in) 
		: data (data_in), eof_flag (false) 
	{}
	/**@brief Move construct for Fasta
	 * @param data 
	 */
	Fasta ( Fasta && other ) 
		: data( std::move(other.data) ), eof_flag (std::move(other.eof_flag)) //Move constructor
	{}
	///@brief Move assignment operator
	///@param other
	Fasta& operator=(Fasta && other) //Move assignment operator
	{
		data = other.data;
		eof_flag = other.eof_flag;
		//std::cerr << "Move Assignment operator" << std::endl;
		return *this;
	}
	/**@brief construct a Fasta <TUPLETYPE> with an EOFFLAG in bool format
	 * @param EofFlag in bool format
	 */
	//end of file flag for fasta_reader io_end()
	Fasta (bool EofFlag) 
		: eof_flag (EofFlag)
	{}
	/**@brief overload operator<< 
	 * @param out ostream object
	 * @param s Fasta format data
	 * @return struct Fasta <TUPLETYPE>
	 */
	friend std::ostream& operator<< (std::ostream& out, const Fasta& in)
	{
		std::string result;
		TupleUtility< std::tuple<std::string, std::string>, std::tuple_size< TupleType >::value >::PrintTuple(in.data, result);
		out << ">" << result;
		return out;
	}
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		TupleUtility < TUPLETYPE, std::tuple_size<TUPLETYPE>::value > 
			:: SerializeTuple (data, ar, version);
		ar & eof_flag;
	}
	
	inline std::string& getName()
	{
		return std::get<0>(data);
	}
	inline std::string& getSeq()
	{
		return std::get<1>(data);
	}
	inline std::string& getName2()
	{
		return std::get<0>(data);
	}
	inline std::string getQuality()
	{
		return std::string( std::get<1>(data).size(), 'I' );
	}
};
#endif
