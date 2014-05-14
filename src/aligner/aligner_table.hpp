/**
 *  @file aligner_table.h
 *  @brief Aligner Searcher與Indexer，所必須的 table，包含建表以及使用方法
 *  @author C-Salt Corp.
 */
#ifndef ALIGNER_TABLE_HPP_
#define ALIGNER_TABLE_HPP_
#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>
#include <utility>
#include <memory>

#include "../constant_def.hpp"
#include "../compression/abit.hpp"
#include "../compression/jbit.hpp"

#include "boost/serialization/utility.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"
#include "boost/iostreams/filtering_stream.hpp"
#include "boost/iostreams/device/file.hpp"
#include "boost/iostreams/filter/zlib.hpp"
#include "boost/serialization/map.hpp"
#include "boost/serialization/vector.hpp"

/**
 * @class Aligner_table< int AlignerType >
 * @brief Aligner_table 範型版本，可以特化不同 Aligner 所需要的 table, e.g.,BWT or BLAST
 * @tparam AlignerType (emun)，用來特化用的參數
 */

template<int AlignerType>
class Aligner_table
{};


/**
 * @class Aligner_table< Aligner_types::BWT_Aligner >
 * @brief Aligner Table 的 BWT 特化版本，專門處理建表與基本搜尋會使用的方法
 * @tparam AlignerType (emun)，BWT_Aligner 版本
 */

template<>
class Aligner_table<Aligner_types::BWT_Aligner>
{
private:
	/**
	 * @brief char_size 為內部使用的變數，256，決定一些表的大小，此參數不可變動
	 */
	INTTYPE char_size;
	
	/**
	 * @brief powerV 為內部使用變數，為2的n次方，為了加速計算
	 */
	INTTYPE powerV;
	
	/**
	 * @brief first_location 為內部使用變數，因為$在壓縮bwt後資訊會消失，也就損失suffix index 在起始的位置，因此要特別記錄
	 */
	INTTYPE first_location;
	
	/**
	 * @brief occ_char 為內部使用變數，{'A', 'C', 'G', 'T'}，方便使用
	 */
	std::vector<int> occ_char;
	
	/**
	 * @brief mtable 為內部使用變數，256大小，存放字元對應欲壓縮的數字，A=0, C=1, G=2, T=3
	 */
	std::vector<INTTYPE> mtable;
	
	/**
	 * @brief jbwt_idx_char 為內部使用變數，為了還原壓縮 jbwt，以及方便計算所記的表
	 */
	std::string jbwt_idx_char;
	
	/**
	 * @brief jbwt 為內部使用變數，指向 JBit class
	 */
	std::shared_ptr<JBit> jbwt;

	
public:
	
	typedef ABSequence<std::string> SEQTYPE;
	
	static int is_table_loaded;
	
	/**
	 * @brief interval 為建表時的間隔大小，間隔越小記憶體使用越大，搜尋越慢
	 */
	INTTYPE interval;
	
	/**
	 * @brief bwt 舊版用來儲存 bwt 字串的變數，現在已經捨棄不用。改用 jbwt 壓縮過的
	 */
	std::string bwt;
	
	/**
	 * @brief BWT 演算法中，要復原字串時所記的 C 表，記錄排序後，ACGT在suffix第一個字的起始位置
	 */
	std::vector<INTTYPE> c_function;
	
	/**
	 * @brief BWT 演算法中，要復原字串時所記的 C 表，不一樣的是，這邊不只記錄一個字，而是記錄一串字(12字，窮舉)的起始位置
	 */
	std::vector<INTTYPE> c_functions;
	
	/**
	 * @brief BWT 演算法中，要復原字串時所記的 OCC 表，記錄某個 bwt index 之上，累積有幾個 A, C, G or T
	 */
	std::vector< std::vector<INTTYPE> > occ_function;
	
	/**
	 * @brief BWT 演算法中，要復原字串時所記的 OCC 表，因為最後採用壓縮 jbwt，所以累積字元 A, C, G or T，需要重新記錄，為三維陣列
	 * 1d: size=256 但只存放 ACGT的ascii，2d
	 * 2d: size=5 要查詢的壓縮字串中的第n個字，前累積有C字幾個, 01234, 4比較特殊，就是下一個壓縮字的 0 所記的數
	 * 3d: size=256 壓縮字元，256個字
	 */
	std::vector< std::vector< std::vector<INTTYPE> > > occ_jbwt;
	
	/**
	 * @brief BWT 演算法中，記錄最後要將 bwt index 轉成真正在 Genome 上對應的 location。此表用來加速，避免 trace back 太多字。
	 */
	std::vector< std::pair<INTTYPE, INTTYPE> > location_table;
	
	/**
	 * @brief 超快速找 location表，此表記錄某個 bwt index trace back 幾次就一定在location表中查的到對應的位置，可以避免無謂的 binary search
	 */
	std::vector<uint8_t> fbwt;
	
	/**
	 * @brief 記錄壓縮後的 bwt，一個 uint8_t 代表四個 A, C, G or T 
	 */
	std::vector<uint8_t> jbwt_seq;
	
	/**
	 * @brief FIXME: this is just temperarily storing the real size of the genome
	 */
	INTTYPE _realSize = 0;
	
	/**
	 * @brief starting site -> chromosome
	 */
	std::map <INTTYPE, std::string> chr_start_pos {};
	
	/**
	 * @brief chr -> chr size
	 */
	std::map <std::string, INTTYPE> chr_length {};
	
	/**
	 * @brief unambiguous segment sequence starting position of each chr
	 */
	std::map <INTTYPE, INTTYPE > chr_umbiguous_starting_length {}; 
	
	
	Aligner_table(INTTYPE iv)
		: interval(iv)
		, char_size(256) 
		, c_function(char_size, 0)//, occ_function(bwt_size/interval+1, std::vector<INTTYPE>(char_size, 0))//, location_table(bwt_size/interval+1,0)
		, c_functions()
		, mtable(256)
		, occ_jbwt(256, std::vector< std::vector<INTTYPE> >(5,std::vector<INTTYPE>(256,0) ) )
		, occ_char({'A','C','G','T'})
	{
		set_interval(interval);
		mtable['A']=0;
		mtable['C']=1;
		mtable['G']=2;
		mtable['T']=3;
	}
	Aligner_table()
		: interval(0)
		, char_size(256) 
		, c_function(char_size, 0)//, occ_function(bwt_size/interval+1, std::vector<INTTYPE>(char_size, 0))//, location_table(bwt_size/interval+1,0)
		, c_functions()
		, mtable(256)
		, occ_jbwt(256, std::vector< std::vector<INTTYPE> >(5,std::vector<INTTYPE>(256,0) ) )
		, occ_char({'A','C','G','T'})
	{
		mtable['A']=0;
		mtable['C']=1;
		mtable['G']=2;
		mtable['T']=3;
	}
	
	void readChrStartPos (const std::string& chrStartPosFile);

	void readChrLen (const std::string& chrLenFile);
	void readNPosLen (const std::string& fileName);
	void set_interval(INTTYPE iv);
	
	/**
	 * @fn void Aligner_table::using_jbwt (void)
	 * @brief 建構 occ_jbwt 與 jbwt_idx_char 表
	 * @return void
	 */
	void using_jbwt();

	/**
	 * @fn void Aligner_table::saveTable (void)
	 * @brief 將所有Aligner會用到的表儲存archive到檔案
	 * @param[in] filename 儲存的檔案名稱
	 * @return void
	 */
	void saveTable(std::string filename);
	
	/**
	 * @fn void Aligner_table::readTable (void)
	 * @brief 將所有Aligner會用到的表從檔案讀取進記憶體
	 * @param[in] filename 讀取的檔案名稱
	 * @return void
	 */
	void readTable(std::string filename);
	
	/**
	 * @fn void Aligner_table::createAllTable (void)
	 * @brief 會根據排序過後的bwt index 與原始genome，建構出所有 Searcher 需要用到的表
	 * @tparam SEQTYPE genome sequence 型別
	 * @param[in] seq 原始 genome sequence
	 * @param[in] filenames split sort後，儲存起來的bwt index檔案集合
	 * @return void
	 */
	void createAllTable(SEQTYPE &seq, std::vector<std::string>& filenames);


	/**
	 * @fn void Aligner_table::str_idx_compare (void)
	 * @brief 加速用比較 str idx 大小，內部使用
	 * @tparam SEQTYPE genome sequence 型別
	 * @param[in] a sequence a 位置為起始
	 * @param[in] b sequence b 位置為起始
	 * @param[in] len 比較幾個字，超過為平手
	 * @return bool
	 */
	bool str_idx_compare(SEQTYPE &seq, INTTYPE a, INTTYPE b, INTTYPE len);
	
	/**
	 * @fn void Aligner_table::get_c (INTTYPE)
	 * @brief 查C表，回傳某字元起始的 bwt index 
	 * @param[in] i bwt index
	 * @return INTTYPE bwt index 
	 */
	INTTYPE get_c(INTTYPE i) const;
	
	/**
	 * @fn void Aligner_table::get_c (char)
	 * @brief 查C表，回傳某字元起始的 bwt index 
	 * @param[in] c A,C,G or T 字元
	 * @return INTTYPE ，查C表，回傳某字元起始的 bwt index 
	 */
	INTTYPE get_c(char c) const;
	
	/**
	 * @fn void Aligner_table::get_jbwt_char (INTTYPE)
	 * @brief 解壓縮。從 jbwt 壓縮字串，取得某 bwt index 位置的字元
	 * @param[in] i bwt index
	 * @return char 字元 A, C, G or T
	 */
	char get_jbwt_char(INTTYPE i) const;
	
	/**
	 * @fn void Aligner_table::get_occ_using_jbwt (INTTYPE, char, int)
	 * @brief 壓縮版本，取得某bwt index位置之上有幾個 A, C, G or T。先查表，取得最近距離數值後，再重新累加計算
	 * @param[in] i bwt index
	 * @param[in] c 看某個字，沒給值則為 bwt index上的字
	 * @return INTTYPE 累積數量
	 */
	INTTYPE get_occ_using_jbwt(INTTYPE i, char c);
	
	/**
	 * @fn void Aligner_table::get_occ (INTTYPE, char, int)
	 * @brief 未壓縮版本，取得某bwt index位置之上有幾個 A, C, G or T。先查表，取得最近距離數值後，再重新累加計算
	 * @param[in] i bwt index
	 * @param[in] c 看某個字，沒給值則為 bwt index上的字
	 * @return INTTYPE 累積數量
	 */
	INTTYPE get_occ(INTTYPE i, char c) const;
	
	/**
	 * @fn void Aligner_table::back_tracking_using_jbwt (INTTYPE)
	 * @brief 壓縮版本，解壓縮BWT，往前 trace 一個字
	 * @param[in] i bwt index
	 * @return INTTYPE 前面一個字的新 bwt index
	 */
	INTTYPE back_tracking_using_jbwt(INTTYPE i);
	
	/**
	 * @fn void Aligner_table::back_tracking (INTTYPE)
	 * @brief 未壓縮版本，解壓縮BWT，往前 trace 一個字
	 * @param[in] i bwt index
	 * @return INTTYPE 前面一個字的新 bwt index
	 */
	INTTYPE back_tracking(INTTYPE i) const;
};

#endif
