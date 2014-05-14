#ifndef GENOME_PRE_HANDLER_HPP_
#define GENOME_PRE_HANDLER_HPP_

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <string>

#include "../file_reader.hpp"

struct Segment 
{
public:
	INTTYPE _offset {0};
	INTTYPE _len {0};
	std::string _seq {};
	Segment& operator = (const Segment&);
public:
	explicit Segment (std::istream& is) {
		bool finish = false;
		char c;

		while (is.peek() !='>' && is.good () && !finish) {
			c = is.get ();
			switch (c) {
			case 'A': case 'a': case 'C': case 'c':
			case 'G': case 'g': case 'T': case 't': case 'U': case 'u':
				is.putback(c); finish = true; break;
			case 'N': case 'n':
				++ _offset; break;
			case '\r': case '\n': case ' ':
				break;
			default:
				std::cerr << "unexpected char: " << c << std::endl; break;
			}
		}

		finish = false;

		while (is.peek () != '>' && is.good() && !finish) {
			c = is.get ();
			switch(c) {
			case 'A': case 'a': case 'C': case 'c':
			case 'G': case 'g': case 'T': case 't': case 'U': case 'u':
				++ _len;
				_seq += std::toupper (c);
				break;
			case 'N': case 'n':
				is.putback (c);
				finish = true;
				break;
			case '\r': case '\n': case ' ':
				break;
			default:
				std::cerr << "unexpected char: " << c << std::endl; break;
			}
		}
		//std::cerr << "_len:\t" << _len << "\t_offset:\t" << _offset << std::endl;
	}
	
	explicit Segment (std::string& is, uint64_t &idx) {
		bool finish = false;
		char c;

		while (!finish && idx<is.size()) {
			c = is[idx];
			switch (c)
			{
				case 'A': case 'a': case 'C': case 'c':
				case 'G': case 'g': case 'T': case 't': case 'U': case 'u':
					finish = true;
					break;
				case 'N': case 'n':
					++ _offset;
					++idx;
					break;
				case '\r': case '\n': case ' ':
					++idx;
					break;
				default:
					std::cerr << "unexpected char: " << c << std::endl; 
					++idx;
					break;
			}
		}

		finish = false;

		while (!finish && idx<is.size()) {
			c = is[idx];
			switch(c)
			{
				case 'A': case 'a': case 'C': case 'c':
				case 'G': case 'g': case 'T': case 't': case 'U': case 'u':
					++ _len;
					_seq += std::toupper (c);
					++idx;
					break;
				case 'N': case 'n':
					finish = true;
					break;
				case '\r': case '\n': case ' ':
					++idx;
					break;
				default:
					++idx;
					std::cerr << "unexpected char: " << c << std::endl;
					break;
			}
		}
		//std::cerr << "_len:\t" << _len << "\t_offset:\t" << _offset << std::endl;
	}

	Segment (const Segment&) = delete ;

	Segment (Segment&& other):
		_offset (other._offset),
		_len (other._len),
		_seq {std::move (other._seq)}
	{}

	Segment& operator = (Segment&& other) {
		if (this != &other) {
			_offset = other._offset ;
			_len = other._len;
			_seq.swap (other._seq);
		}
		return *this;
	}
}; /* end of class Segment definition */



/* definition of Genome_pre_handler class */
template < int GenomeStrandType, template < typename T > class FORMAT, typename TUPLETYPE, SOURCE_TYPE STYPE >
class Genome_pre_handler 
	//:public FileReader <ParallelTypes::NORMAL, Fasta , std::tuple <std::string, std::string>, SOURCE_TYPE::IFSTREAM_TYPE >
	:public FileReader_impl <FORMAT , TUPLETYPE, STYPE >
{
private:
	std::string _name {};
	std::vector<Segment> _sequences {}; /// each segment is N...NACGT...ACGT until next N
	INTTYPE _length {0};
	INTTYPE _lengthNoN {0};
	Genome_pre_handler& operator = (const Genome_pre_handler&);
	
public:
	typedef Fasta<> format_type;
	class badGenome_pre_handler {};
	Genome_pre_handler () = default ;
	//
	explicit Genome_pre_handler (std::vector<std::string> is)
		:FileReader_impl <FORMAT , TUPLETYPE, STYPE > (is)
	{}
	
	explicit Genome_pre_handler (std::vector<std::string> is, std::function<void(void)> cb )
		:FileReader_impl <FORMAT , TUPLETYPE, STYPE> (is)
	{}
	
	Genome_pre_handler (const Genome_pre_handler&) = delete;
	Genome_pre_handler (Genome_pre_handler&& other):
		_name (std::move (other._name)),
		_sequences (std::move (other._sequences))
	{}
	Genome_pre_handler& operator = (Genome_pre_handler&& other) {
		if (this != &other) {
			_name.swap (other._name);
			_sequences.swap (other._sequences);
		}
		return *this;
	}

	void run_handler (std::string prefixName, std::string &seq)
	{	
		/* running accumulator recording the length of each chr */
		INTTYPE tempLen {0}, accumulatedLength {0};
		/* for concatenated seq */
		std::map <INTTYPE, INTTYPE> NPosLen { };
	
		/* file to store which regions has which chr*/
		std::ofstream chrStartPos {prefixName + "chrStart"};
		/* file to store the length of each chr */
		std::ofstream chrLen {prefixName + "chrLen"};
		/* read in each fasta and make two string */
		
		
		parser_combo(
			[&seq, &chrStartPos, &chrLen, &accumulatedLength, &tempLen, &NPosLen, this]()
			{
				/* store start position of each chr */
				chrStartPos << this->getName () << '\t' << accumulatedLength << '\n';
				/* get chr length */
				tempLen = this->getLengthNoN ();
				/* store chr length */
				chrLen << this->getName () << '\t' << tempLen << '\n';
				/* update accumulated length */
				accumulatedLength += tempLen;
				/* update NPosLen */
				this->updateNpos (NPosLen);
				seq += this->getSeqNoN ();	
			}
		);
		chrStartPos.close ();
		chrLen.close ();
		
		
		/* resize to enough space for the reverse complemetary sequence and a $ sign */
		seq.resize (seq.size () * 2 + 1); // TODO: resize does mallocating the extra space and also initialization, the later is not necessary
		auto iter = seq.begin ();
		std::advance (iter, (seq.size ()-1)/2);
		auto iter2 = iter;
		--iter2;
		do {
			switch (*iter2) {
			case 'A':
				*iter = 'T'; break;
			case 'T':
				*iter = 'A'; break;
			case 'G':
				*iter = 'C'; break;
			case 'C':
				*iter = 'G'; break;
			}
			++iter;
		} while (iter2-- != seq.begin ());
		*iter = '$';
		/* writing NPosLen to file */
		
		{
			
			std::ofstream fp(prefixName + "NposLen.z", std::ios::binary);
			boost::archive::binary_oarchive archive_fp( fp );
			archive_fp & NPosLen;
			fp.close();
			//std::cerr << "NPosLen size : " << NPosLen.size() << std::endl;
			//for (auto i : NPosLen)
			//	std::cerr << i.first << ":" << i.second << std::endl;

		}
	}

	void parser_combo(std::function<void(void)> callback)
	{
		std::cerr << "file_number  " << this->file_num_ <<'\n';
		for (auto file_idx=0; file_idx != this->file_num_; ++file_idx)
		{
			//std::cerr << "current index " << file_idx<<'\n';
			while (true)
			{
				_sequences.clear();
				_length = 0;
				_lengthNoN = 0;
				
				Fasta < std::tuple <std::string, std::string > > object = this->get_next_entry (file_idx);
				bool flag = true;
				if ( object.eof_flag )
					break;
				_name = std::get<0> (object.data);
				
				std::string &seq = std::get<1> (object.data);
				
				uint64_t seq_idx(0);
				while(true && seq_idx < seq.size())
				{
					_sequences.emplace_back(seq,seq_idx);
					auto& _seg = _sequences.back ();
					_length += _seg._len + _seg._offset ;
					_lengthNoN += _seg._len;
				}
				
				callback();
			}
		}
	}
	
	
	void updateNpos (std::map <INTTYPE, INTTYPE>& NposLen) const {
		auto seg = _sequences.begin ();
		auto stopIter = _sequences.end(); advance (stopIter, -1);
		std::pair <std::map <INTTYPE, INTTYPE>::iterator, bool> lastPos2  {};
	// adding first segment
		/* first chromosome */
		if (NposLen.empty ()) {
			lastPos2 = NposLen.insert ( std::make_pair (0, seg->_offset) ); /// seg->_offset is consumed in current cycle
//			std::cerr << "Just inserted:\t" << (lastPos2.first)->first <<'\t' << (lastPos2.first)->second << '\n';
		}
		/* non first chromosome */
		else {
			auto lastPos = NposLen.rbegin ();  /// get the last entry, update its N
			lastPos->second = seg->_offset; /// the ->first of last entry is gonna used by this new Genome_pre_handler, but we won't get the information of the N until now, hereby we update it
			lastPos2 = NposLen.insert ( std::make_pair (lastPos->first, seg->_offset) ); /// insert a new one, update its _offset. (but _len won't get updated until next segment...)
		}
	// adding middle segment, every time we meet a new segment, we immediately use its _offset. But its _len won't get used until next circle
		while (seg != stopIter) {
//			std::cerr << "(lastPos2.first)->first:\t" << (lastPos2.first)->first << '\n';
//			std::cerr << "seg->_len:\t" << seg->_len << '\n';
			auto tmp = (lastPos2.first)->first + seg->_len;
			lastPos2 = NposLen.insert (std::make_pair ( tmp , (lastPos2.first)->second + (++seg)->_offset  ));
//			std::cerr << "Just inserted:\t" << (lastPos2.first)->first <<'\t' << (lastPos2.first)->second << '\n';
		}
		if (seg->_len != 0) { /// if last segment has ATGC, we have to update their _len
			NposLen.insert (std::make_pair ( (lastPos2.first)->first + seg->_len,  (lastPos2.first)->second )); // this will be updated in next chr
		}

//		for (const auto x : NposLen) {
//			std::cerr << x.first << '\t' << x.second << std::endl;
//		}
	}

	friend std::ostream& operator << (std::ostream& os, const Genome_pre_handler& Genome_pre_handler) {
		os << '>' << Genome_pre_handler._name << '\n';
		for (const Segment& seg : Genome_pre_handler._sequences) {
			for (int i = 0; i < seg._offset; ++i) {
				os << 'N';
			}
			os << seg._seq;
		}
		os << '\n';
		return os;
	}
	std::string getName () const {
		return _name;
	}
	std::string getSeq () const {
		std::string _sequence {};
		for (const Segment& seg : _sequences) {
			for (int i = 0; i < seg._offset; ++i) {
				_sequence += 'N';
			}
			_sequence += seg._seq;
		}
		return _sequence;
	}
	std::string getSeqNoN () const {
		std::string _sequence;
		for (const auto& seg : _sequences) {
			_sequence += seg._seq;
		}
		return _sequence;
	}
	std::string getReverseSeq () const {
		std::string tmp {this->getSeqNoN()};
		return std::string {tmp.crbegin(), tmp.crend()};
	}

	std::string getReverseComplementSeq () const {
		std::string _sequence {this->getReverseSeq()};
		for (auto iter = _sequence.begin(); iter!= _sequence.end(); ++ iter) {
			switch (*iter) {
			case 'A':
				*iter = 'T'; break;
			case 'T':
				*iter = 'A'; break;
			case 'G':
				*iter = 'C'; break;
			case 'C':
				*iter = 'G'; break;
			case 'N':
				break;
			default:
				throw badGenome_pre_handler ();
			}
		}
		return _sequence;
	}
	INTTYPE getLength () const {
		return _length;
	}
	INTTYPE getLengthNoN () const {
		return _lengthNoN;
	}
}; /* end of class Genome_pre_handler definition */

#endif
