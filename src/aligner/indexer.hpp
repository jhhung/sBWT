/**
 *  @file indexer.h
 *  @brief Indexer of Aligner，建表 
 *  @author JHH Corp.
 */
#ifndef INDEXER_HPP_
#define INDEXER_HPP_

#include <string>
#include <vector>

#include "../compression/abit.hpp"
#include "../compression/jbit.hpp"
#include "../constant_def.hpp"
#include "../mkq_sort.hpp"
#include "../difference_cover.hpp"
#include "../split_sort.hpp"
#include "aligner_table.hpp"

/**
 * @class Indexer
 * @brief 範型版本 Indexer，主要處理 Genome 建 Aligner Searcher 所需要的表
 *  
 * @tparam AlignerType (enum) Aligner 版本，可以特化不同 Aligner 所需要的 Indexer, e.g.,BWT or BLAST
 * @tparam ParallelType (enum) 平行方法, e.g., M_T, NORMAL
 * @tparam AlignerTableType (class) Indexer 所需要的 Aligner table type
 * @tparam GenomePreHandlerType (class) 前處理 Genome sequence type，主要記錄 NNN..., >Name ...，並給予純 ACGT sequence 
 * @tparam BuildIndexReturnType (typename) 建完index後回傳型別
 * @tparam ARGS (typename)
 *  
 */
template<
	Aligner_types AlignerType
	,ParallelTypes ParallelType
	,class AlignerTableType
	,class GenomePreHandlerType
	,typename BuildIndexReturnType
	,typename... ARGS
>
class Indexer
{};

/**
 * @class Indexer
 * @brief 特化版本 BLAST Indexer，建構 Genome sequence BLAST
 *  
 * @tparam Aligner_types::BLAST_Aligner (enum) Aligner 版本，這邊為BLAST特化版板
 * @tparam ParallelType (enum) 平行方法, e.g., M_T, NORMAL
 * @tparam AlignerTableType (class) Indexer 所需要的 Aligner table type
 * @tparam GenomePreHandlerType (class) 前處理 Genome sequence type，主要記錄 NNN..., >Name ...，並給予純 ACGT sequence 
 * @tparam BuildIndexReturnType (typename) 建完index後回傳型別
 * @tparam ARGS (typename)
 *  
 */
template< 
	ParallelTypes ParallelType
	,class AlignerTableType
	,class GenomePreHandlerType
	,typename BuildIndexReturnType
	,typename... ARGS
>
class Indexer<
	Aligner_types::BLAST_Aligner
	,ParallelType
	,AlignerTableType
	,GenomePreHandlerType
	,BuildIndexReturnType
	,ARGS... 
>
{};


/**
 * @class Indexer
 * @brief 特化版本 BWT Indexer，建構 Genome sequence BWT
 *  
 * @tparam Aligner_types::BWT_Aligner (enum) Aligner 版本，這邊為BWT演算法特化版板
 * @tparam ParallelType (enum) 平行方法, e.g., M_T, NORMAL
 * @tparam AlignerTableType (class) Indexer 所需要的 Aligner table type
 * @tparam GenomePreHandlerType (class) 前處理 Genome sequence type，主要記錄 NNN..., >Name ...，並給予純 ACGT sequence 
 * @tparam BuildIndexReturnType (typename) 建完index後回傳型別
 * @tparam SeqType (class) 這邊要傳入 Genome sequence 型別，壓縮與否
 * @tparam DcsType (class) 建構 DCS表的方法以及實作，需要能提供一個 compare的 function
 * @tparam SplitSortType (class) 主要 sort genome suffix 的方法，並且決定是否一次 sort完或分小部分 sort(省記憶體)
 * @tparam ARGS (typename)
 *  
 */
template<
	ParallelTypes ParallelType
	,class AlignerTableType
	,class GenomePreHandlerType
	,class BuildIndexReturnType
	,class SeqType
	,class DcsType
	,class SplitSortType
>
class Indexer<
	Aligner_types::BWT_Aligner
	,ParallelType
	,AlignerTableType
	,GenomePreHandlerType
	,BuildIndexReturnType
	,SeqType
	,DcsType
	,SplitSortType
>
{
public:

	/**
	 * @brief Indexer 所需要的 Aligner table type，此為 reference，為了與 searcher 共用
	 */
	AlignerTableType &abwtt_;
	
	/**
	 * @brief Genome 前處理完後，純ACGT序列，存放容器，也是將會是SeqType的容器
	 */
	std::string sequence;
	
	/// @brief 計時用參數
	clock_t start, stop;
	
	
	Indexer(AlignerTableType& table)
		:abwtt_(table)
	{};
	
	/**
	 * @fn void Indexer::build (void)
	 * @brief 留給外面使用的功能，build，包含所有的建表動作
	 * @param[in] filelist Genome檔案清單
	 * @param[in] prefix_name 要建表的檔名前贅字元
	 * @param[in] pre_sort_length Genome suffix 只sort前N個字，決定建DSC表的長度大小，也決定Genome suffix只需sort前N個字
	 * @param[in] interval_size 建表時的間隔，攸關表的大小，記憶體大小，搜尋速度
	 * @return BuildIndexReturnType
	 */
	BuildIndexReturnType
	build(std::vector<std::string> filelist , std::string prefix_name="", INTTYPE pre_sort_length=256, INTTYPE interval_size=64)
	{
		if (prefix_name.back () != '.') {
			prefix_name += '.';
		}
		/// @brief 建構GenomePreHandlerType，一種 file parser，給予一個檔案清單
		/// @brief 前處理 Genome sequence type，主要記錄 NNN..., >Name ...，並給予純 ACGT sequence 
		GenomePreHandlerType handler(filelist);
		handler.run_handler(prefix_name, sequence);
		std::cout << "seq len " << sequence.size() << std::endl;
		/// @brief 建構SeqType壓縮方法，使用 reference方式，原始sequence不能被消滅
		SeqType genome ( sequence );
		
		return build_index(genome, pre_sort_length, interval_size, prefix_name);
	}

private:

	/**
	 * @fn void Indexer::build_index (void)
	 * @brief 建構 Aligner table，主要分為三個部分，1: 建構 DCS表, 2:Sort suffix with dcs, 3:Build index with sorted suffix array
	 * @param[in] sequence 乾淨的Genome序列
	 * @param[in] pre_sort_length Genome suffix 只sort前N個字，決定建DSC表的長度大小，也決定Genome suffix只需sort前N個字
	 * @param[in] interval_size 建表時的間隔，攸關表的大小，記憶體大小，搜尋速度
	 * @param[in] prefix_name 要建表的檔名前贅字元
	 * @return void
	 */
	BuildIndexReturnType
	build_index(SeqType &sequence, INTTYPE pre_sort_length=256, INTTYPE interval_size=64, std::string prefix_name="")
	{
		
		/// @brief part 1
		start = clock();
		
		/// @brief 建構 DCS 表，主要必須提供 dcs::compare function
		DcsType dcs(sequence, pre_sort_length);
		
		stop = clock();
		std::clog << "===================================" << std::endl;
		std::clog << "Create dcs table, time:" << double(stop - start) / CLOCKS_PER_SEC << std::endl;
		std::clog << "===================================" << std::endl;
		
		/// @brief part 2
		start = clock();
		
		/// @brief sort suffix using split sort
		SplitSortType split_sort (sequence, pre_sort_length);
		
		/// @brief 當初測試tuple type留下的功能，提供指定 tuple type哪個元素的方法
		split_sort.getTableV = [](INTTYPE& a) {return a;};
		
		/// @brief Copy dcs::compare 為 split_sort.dcs_compare，split sort 會在 sort到 pre_sort_length 後，使用此function 來比較大小
		split_sort.dcs_compare = std::bind( &DcsType::compare, dcs, std::placeholders::_1, std::placeholders::_2);	
		
		/// @brief Start split sort with tparam size, default 30,000,000
		split_sort.split_by_tparam_size();
		
		/// @brief copy sorted result (bwt index, filelist) to seq_table(bwt)
		std::vector<std::string> filenames = split_sort.archive_name;
		
		/// @brief release dcs memory
		dcs.release();
		split_sort.release();
		
		stop = clock();
		std::clog << "===================================" << std::endl;
		std::clog << "Create BWT, time:" << double(stop - start) / CLOCKS_PER_SEC << std::endl;
		std::clog << "===================================" << std::endl;
		
		/// @brief part 3
		
		start = clock();
		
		/// @brief 設定建表的間隔
		abwtt_.set_interval(interval_size);
		
		/// ＠brief 開始建表
		abwtt_.createAllTable(sequence, filenames);
		
		/// @brief 寫入前贅字為 prefix_name 的檔案
		abwtt_.saveTable( prefix_name + "t_table.bwt");
		
		stop = clock();
		std::clog << "===================================" << std::endl;
		std::clog << "Create abwt Table, time:" << double(stop - start) / CLOCKS_PER_SEC << std::endl;
		std::clog << "===================================" << std::endl;
		
		/// @brief 清除split sort 留下的暫存檔案
		split_sort.clean_up ();
		
		return;
	}

};

#endif
