#ifndef ALIGNER_HPP_
#define ALIGNER_HPP_

#include "aligner_trait.hpp"
/*
template<class Alinger_trait_>
class Aligner
	:public Alinger_trait_::Aligner_type
{};
*/


/**
 * @class Aligner 
 * @brief 組裝 Aligner，繼承 Indexer, Searcher，並以 Alinger_table為參數建構他們
 * @tparam Alinger_trait_ class Aligner_trait，需要有三大組成 : Alinger_trait_::IndexerType, Alinger_trait_::SearcherType, Alinger_trait_::AlignerTableType
 */
template<class Alinger_trait_>
class Aligner
	:public Alinger_trait_::IndexerType
	,public Alinger_trait_::SearcherType
{
public:

	/**
	 * @brief 用 Alinger_trait_::AlignerTableType 建構空的 aligner_table_
	 */
	static typename Alinger_trait_::AlignerTableType aligner_table_;
	
	/**
	 * @fn void Aligner::Aligner (void)
	 * @brief Aligner 建構子，用 aligner_table_ 建構繼承來的 Indexer 與 Searcher，因為他們需要共用 aligner_table
	 * @return void
	 */
	Aligner()
		:Alinger_trait_::IndexerType(aligner_table_)
		,Alinger_trait_::SearcherType(aligner_table_)
	{};
};

template<class Alinger_trait_>
typename Alinger_trait_::AlignerTableType Aligner<Alinger_trait_>::aligner_table_;

#endif