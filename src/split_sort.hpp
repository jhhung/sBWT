/**
 *  @file split_sort.h
 *  @brief 負責分群、分批sort suffix array
 *  @author JHH Corp.
 */
#ifndef SPLIT_SORT_HPP_
#define SPLIT_SORT_HPP_
#include <fstream>
#include <map>
#include <list>
#include <random>
#include <algorithm>
#include <functional>
#include <boost/ref.hpp>
#include <utility>
#include <cstdlib>
#include <queue>
#include <tuple>
#include <ctime>
#include <memory>

#include "boost/serialization/vector.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/filesystem.hpp"

#include "mkq_sort.hpp"
#include "compression/abit.hpp"




/**
 * @class Split_sort 
 * @brief 負責分群、分批sort suffix array並且存入檔案，但是不負責 Sort 的方法，方法是由外面傳進來的
 *  
 * @tparam SeqType (class) 要sort的原始字串, sequence
 * @tparam ContainerType (class) 要 sort的 suffix index 所存放的容器型別
 * @tparam SortType (class) Sort的方法
 * @tparam GROUPSIZE (int) 分群分批 sort的每群大約大小，與記憶體、速度有關，default 30,000,000
 *  
 */
template <class SeqType, class ContainerType, class SortType, int GROUPSIZE = 30000000>
class Split_sort
{
private:
	
	/// @brief 原本的sequence
	SeqType &seq;
	
	/// @要存放的 suffix index 容器，vec = self_vec
	ContainerType &vec, self_vec;
	
	/**
	 * @brief 平均每群大約的大小，此數值會決定大需的分群數量 split_num
	 */
	INTTYPE average_size;
	
	/**
	 * @brief 記錄分群的數量
	 */
	INTTYPE split_num;
	
	/**
	 * @brief 將會 sampling 的數量
	 */
	INTTYPE random_num;
	
	/**
	 * @brief 存放隨機 sampling 的 suffix idx 容器
	 */
	std::vector<INTTYPE> random_table;
	
	/**
	 * @brief 存放最後真正用來分類的 sampling suffix idx，數量是 split_num -1
	 */
	std::vector<INTTYPE> split_table;
	
	/**
	 * @brief 存放 suffix index 的容器，通常為 vector< vector< uint_32t > >，group< suffix_index< int > >
	 */
	std::vector< ContainerType > SeqTables;
	
	/**
	 * @brief sort 的字串深度，通常為 dcs size in bwt
	 */
	INTTYPE limit_depth;

public:

	
	/**
	 * @brief 給特殊 container 取值用的，也就是非預設 vector< uint_32t >，像是 vector< tuple<uint32_t, uint32_t> > 時可以用到
	 */
	std::function<INTTYPE(typename ContainerType::value_type&)> getTableV;
	
	/**
	 * @brief 在 sort 字串深度達到限制後，最後使用的比較大小，通常 dcs compare table 可以直接分出勝負。sbwt 則不需要分出勝負，就全回傳 true or false 即可
	 */
	std::function<bool(INTTYPE, INTTYPE)> dcs_compare;
	
	/**
	 * @brief 記錄寫出檔案名稱 list 
	 */
	std::vector<std::string> archive_name;
	
	/**
	 * @brief 記錄每個 group 寫入幾次 archive，讀檔的時候會用到
	 */
	std::vector<INTTYPE> object_in_each_archive;
	
	/**
	 * @brief 記錄每個 group 存放幾個 suffix index
	 */
	std::vector<INTTYPE> size_in_each_archive;
	std::vector< std::shared_ptr<std::ofstream> > ofile_list;
	std::vector< std::shared_ptr<std::ifstream> > ifile_list;
	std::vector< std::shared_ptr< boost::archive::binary_oarchive> > oarchive_list;
	std::vector< std::shared_ptr< boost::archive::binary_iarchive> > iarchive_list;
	
	

	Split_sort(SeqType& sequence, ContainerType& table, INTTYPE ld)
		:seq(sequence), vec(table), limit_depth(ld)
	{}
	
	Split_sort(SeqType& sequence, INTTYPE ld)
		:seq(sequence), vec(self_vec), limit_depth(ld)
	{}
	
	template <typename T>
	void free( T & t ) {
			T tmp;
			t.swap( tmp );
	}
	void release()
	{
		free(vec);
		free(self_vec);
		free(random_table);
		free(split_table);
		free(SeqTables);
		free(archive_name);
		free(object_in_each_archive);
		free(size_in_each_archive);
	}
	void split_by_tparam_size()
	{
		split_by_size_init(GROUPSIZE);
	}
	
	void split_by_size_init(const INTTYPE s)
	{
		clock_t start, stop;
		average_size = s;
		split_num = seq.size() / average_size;
		if(split_num == 0) split_num = 1;
		random_num = split_num * 4;
		SeqTables.resize(split_num, std::vector<typename ContainerType::value_type> () );
		
		make_random_table();
		sort_random_table();
		make_split_table();
		
		start = clock();
		std::clog << "classify start:" <<std::endl;
		classify_seq_tables();
		stop = clock();
		std::clog << "classify end: "<< double(stop - start)/CLOCKS_PER_SEC << std::endl;
		
		//INTTYPE i(0);
		//for(auto &K : SeqTables)
		//{
		//	std::clog << "group: "<< i << ":" << K.size() << std::endl;
		//	++i;
		//}
		
		start = clock();
		std::clog << "sort final start:" <<std::endl;
		
		mkq_sort();
		
		stop = clock();
		std::clog << "sort final end: "<< double(stop - start)/CLOCKS_PER_SEC << std::endl;

		//merge();
	}
	
	inline void make_random_table()
	{
		for(INTTYPE i(0); i < random_num; ++i)
		{
			INTTYPE val = std::rand() % seq.size();
			random_table.push_back(val);
		}
		std::sort(random_table.begin(),random_table.end());
		std::unique(random_table.begin(), random_table.end());
	}
	
	inline void sort_random_table()
	{
		
		std::sort(random_table.begin(), random_table.end(),
			[&]( typename ContainerType::value_type A, typename ContainerType::value_type B){
				INTTYPE a(getTableV(A)), b(getTableV(B));
				for(INTTYPE i(0); i <= limit_depth; ++i)
				{	
					if( seq[(a+i)] > seq[(b+i)] )
						return false;
					else if( seq[(a+i)] < seq[(b+i)] )
						return true;
				}
				return dcs_compare(a,b);
			}
		);
	}
	
	
	
	inline void make_split_table()
	{
		
		for(INTTYPE i(1); i < (split_num); ++i)
		{
			INTTYPE idx = random_table.size() / (split_num) * i;
			split_table.push_back(random_table[idx]);
		
		}
		for (INTTYPE i(0); i < (split_num); ++i)
		{
				archive_name.push_back( std::string("split_")+boost::lexical_cast<std::string>(i) );
				ofile_list.push_back( std::shared_ptr<std::ofstream> ( new std::ofstream(archive_name[i], std::ios::binary) ) );
				oarchive_list.push_back( std::shared_ptr<boost::archive::binary_oarchive> ( new boost::archive::binary_oarchive( *ofile_list[i] ) ) );
				object_in_each_archive.push_back(0);
				size_in_each_archive.push_back(0);
		}
		std::cerr << "table size : " <<	SeqTables.size() << std::endl;
	}

	void clean_up () const {
		for(INTTYPE i(0); i < (split_num); ++i) {
			boost::filesystem::remove_all (std::string("split_") + std::to_string(i));
		}
	}

	inline void classify_seq_tables()
	{	
		clock_t start, stop;
		std::vector<INTTYPE>::iterator split_table_it(split_table.begin());
		INTTYPE idx(0);
		
		for(INTTYPE K(0); K<seq.size(); ++K)
		{
			if( K % 67108863 == 0)
				std::cout << "progress: " << K << "/" << seq.size() << std::endl;
			split_table_it = std::lower_bound(split_table.begin(), split_table.end(), K ,
				[&](INTTYPE a, INTTYPE b){
					for(INTTYPE i(0); i <= limit_depth; ++i)
					{	
						if( seq[(a+i)] > seq[(b+i)] )
							return false;
						else if( seq[(a+i)] < seq[(b+i)] )
							return true;
					}
					return dcs_compare(a,b);
				}
			);
			idx = split_table_it - split_table.begin() ;
			SeqTables[idx].push_back(K);
			
			if (SeqTables[idx].size()>1000)
			{
					size_in_each_archive[idx] += SeqTables[idx].size();
					*oarchive_list[idx] & SeqTables[idx];
					object_in_each_archive[idx]++;
					SeqTables[idx].clear();
			};
			
		}
		
		idx=0;
		for(auto &KK : SeqTables)
		{
			size_in_each_archive[idx] += KK.size();
			*oarchive_list[idx] & KK;
			ofile_list[idx]->close();
			std::cerr<<"done outputing: "<<archive_name[idx]<< "size: "<<size_in_each_archive[idx]<<std::endl;
			object_in_each_archive[idx]++;
			KK.clear();
			idx++;
		}
		
		
	}
	
	inline void mkq_sort()
	{
		clock_t start, stop;
		INTTYPE i(0);
		for(ContainerType &K : SeqTables)
		{
			start = clock();
			std::ifstream in( archive_name[i], std::ios::binary );
			std::cerr << "reading: " << archive_name[i] << " " << object_in_each_archive[i] <<std::endl;
			boost::archive::binary_iarchive temp_archive(in);
			for (INTTYPE j=0; j<object_in_each_archive[i]; ++j)
			{
				ContainerType temp_vector;
				temp_archive & temp_vector;
				K.insert( K.end(), temp_vector.begin(), temp_vector.end()	); 
			}
			in.close();
			
			std::clog << "Sort start:" << i << " size : " << K.size() <<std::endl;
			
			SortType sorter(seq, K, limit_depth, dcs_compare);
																							
			sorter.getTableV = [](INTTYPE& a){
					return a;
			};
			sorter.sort(0,K.size(),0);
			
			stop = clock();

			std::cerr << "Sort end: "<< double(stop - start)/CLOCKS_PER_SEC << std::endl;
			
			std::ofstream out( archive_name[i], std::ios::binary );
			boost::archive::binary_oarchive temp_archive_o(out);
			temp_archive_o & K;
			out.close();
			
			
			ContainerType temp;
			K.swap(temp);
			
			std::clog << "Sort end: "<< double(stop - start)/CLOCKS_PER_SEC << std::endl;
			++i;
		}
	}
};


#endif
