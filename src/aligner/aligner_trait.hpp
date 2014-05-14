/**
 *  @file aligner_trait.hpp
 *  @brief 放置 Aligner 會用的 trait，也就是"設定"集，設定 Aligner
 *  @author C-Salt Corp.
 */
 
#ifndef ALIGNER_TRAIT_HPP_
#define ALIGNER_TRAIT_HPP_

#include "../constant_def.hpp"
#include "genome_pre_handler.hpp"
#include "aligner_table.hpp"
#include "searcher.hpp"
#include "indexer.hpp"

/**
 * @class Aligner_trait<> 
 * @brief Aligner_trait 範型版本，此class主要為定義 Aligner參數，主要有三大要件： AlignerTableType, IndexerType, SearcherType
 *  
 * @tparam Aligner_type (enum) Aligner_type 決定 Aligner type，目前有 BWT_Aligner, BLAST_Aligner，預設為 BWT_Aligner。注意！此預設參數會直接讓編譯器選擇 Aligner_trait<Aligner_types::BWT_Aligner> 特化版本
 * @tparam Parallel_type (enum) 平行的方法，預設為 M_T
 * @tparam QueryParserType (class) Searcher 讀檔的 File_parser type，預設為 Fastq
 */
//For Tailor
template< Aligner_types Aligner_type = Aligner_types::BWT_Aligner
				 ,ParallelTypes Parallel_type = ParallelTypes::M_T
				 ,class QueryParserType = FileReader_impl
					< Fastq 
					 ,std::tuple <std::string, std::string, std::string, std::string>
					 ,SOURCE_TYPE::IFSTREAM_TYPE
					>
				>
class Aligner_trait
{};

/**
 * @class Aligner_trait<Aligner_types::SBWT_Aligner>
 * @brief Aligner_trait Aligner_types::SBWT_Aligner 特化版本，定義 SBWT Aligner參數
 *  
 * @tparam Aligner_type::SBWT_Aligner (enum) BWT 特化版本
 * @tparam Parallel_type (enum) 平行的方法
 * @tparam QueryParserType (class) Searcher 讀檔的 File_parser type
 *  
 */
template<ParallelTypes Parallel_type, class QueryParserType >
class Aligner_trait<Aligner_types::SBWT_Aligner, Parallel_type, QueryParserType>
{
public:
	
	/**
	 * @brief 定義必要元件 Aligner table type，這邊使用 BWT table
	 */
	typedef Aligner_table<Aligner_types::BWT_Aligner> AlignerTableType;
	
	/**
	 * @brief 定義必要元件 Indexer build index 完後的 return type
	 */	
	typedef void BuildIndexReturnType;
	
	/**
	 * @brief 定義必要元件 Searcher search 完後的 return type
	 */
	typedef void SearchReturnType;
	
	
	/**
	 * @brief 定義 Indexer 裡面會用到的 SeqType， ABSequence對sequence 做一些處理，像是壓縮，這邊使用的為普通 std::string（不壓縮）
	 */
	typedef ABSequence< std::string > SeqType;
	
	/**
	 * @brief 定義 Indexer 裡面會用到的 ContainerType，也就是存放 BWT index (數字) 的資料型態
	 */
	typedef std::vector< INTTYPE > ContainerType;
	
	/**
	 * @brief 定義 Indexer 裡面會用到的 GenomePreHandlerType，這邊指定 genome 為 fasta 格式。Genome_pre_handler主要為記錄 >name 與 sequence NNNN...，並且回傳只有 ATCG的字串
	 */
	typedef Genome_pre_handler
					< Genome_strand_types::Dual_strand
					 ,Fasta 
					 ,std::tuple <std::string, std::string>
					 ,SOURCE_TYPE::IFSTREAM_TYPE
					> GenomePreHandlerType;
	
	/**
	 * @brief 定義 Indexer 裡面會用到的 DcsType 裡面的 DcsSortType，也就是建構 DCS 表時，所用的 sort type
	 */
	typedef Multikey_quicksort
					<	SeqType,
						ContainerType,
						Record_rank_enable,
						Compare_default,
						Sort_small_n_disable,
						Bucket_sort
					>	DcsSortType;

	/**
	 * @brief 定義 Indexer 裡面會用到的 DcsType，也就是建構 DCS 的方法，這邊不建構 DCS
	 */
	typedef PSEUDO_DCS_FOR_SBWT < SeqType, DcsSortType > DcsType;
	
	/**
	 * @brief 定義 Indexer 裡面會用到的 Split_sort 裡面的 MainSortType, 此為使用 DCS table，與其他綁定參數
	 */
	typedef Multikey_quicksort
					<	SeqType,
						ContainerType,
						Record_rank_disable,
						Compare_default,
						Sort_small_n_enable,
						Bucket_sort
					>	MainSortType;
	
	/**
	 * @brief 定義 Indexer 裡面會用到的 Split_sort, 第四個參數為分群大小，攸關記憶體，預設為 30,000,000
	 */
	typedef Split_sort
					< SeqType
					 ,ContainerType
					 ,MainSortType
					 ,30000000
					> SplitSortType;
	
	/**
	 * @brief 定義必要元件 Indexer 裡面會用到的 Split_sort, 第四個參數為分群大小，攸關記憶體，預設為 30,000,000
	 */
	typedef Indexer
					< Aligner_types::BWT_Aligner
					 ,Parallel_type
					 ,AlignerTableType
					 ,GenomePreHandlerType
					 ,BuildIndexReturnType
					 ,SeqType
					 ,DcsType
					 ,SplitSortType
					> IndexerType;
	
	
	/**
	 * @brief 定義 Searcher 裡面會用到的 QueryParserType, 也就是searcher在讀取reads時所需的 file parser。現在已經提出為 tparam
	 */
	
	/*
typedef FileReader_impl
				< Fastq 
				 ,std::tuple <std::string, std::string, std::string, std::string>
				 ,SOURCE_TYPE::IFSTREAM_TYPE
				> QueryParserType;
	*/
	
	typedef Fastq< std::tuple<std::string,std::string,std::string,std::string> > FormatType;
	
	/**
	 * @brief 定義必要元件 Searcher，這邊是 Tailer
	 */
	typedef Searcher
					< Aligner_types::BWT_Aligner
					 ,Searcher_types::SBWT_exact_match
					 ,Parallel_type
					 ,AlignerTableType
					 //,QueryParserType
					 ,FormatType
					 ,SearchReturnType
					> SearcherType;
	
	
	
};


template<ParallelTypes Parallel_type, class QueryParserType>
class Aligner_trait<Aligner_types::BLAST_Aligner, Parallel_type, QueryParserType>
{};

#endif
