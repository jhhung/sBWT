#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/functional.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <algorithm>
#include "cudaec.cpp"
#include <sys/time.h>
#include <unistd.h>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/lexical_cast.hpp>
#include <memory>

#include "my_genome_pre_handler.hpp"

typedef unsigned short Type;
typedef long long INTTYPE;
#define CARD_MEMORY_LIMIT 6442450944
// #define CARD_MEMORY_LIMIT 1048576
#define SUFFIXLEN 256 //###
//#define SUFFIXLEN 16 //###
#define TBLSTRIDE 32 //###
//#define TBLSTRIDE 2 //###
#define TYPEBITS 16
#define TYPEBYTES sizeof(Type)
#define ENCODEBITS 2
#define ENCODELEN (TYPEBITS/ENCODEBITS)
#define ENCODERATIO (8/ENCODEBITS)
#define MAXREADLEN 128
#define DIMB 16
#define LEN_HANDLE_CHAR 32
#define LEN_CHR 64

#define OVERNUM 1 

#define MAX_DIM_THD 1024
#define MAX_DIM_BLK 65536

//#define MYFILE 1

//#define IS_CREATE_TBL 1
//#define IS_FIND_READS 1

enum DNA_enum{A = 0, C, G, T, DOLLAR};
enum STRAND_opt_enum{FIND_N_STRAND = 0, FIND_P_STRAND, FIND_BOTH_STRAND};
enum SAM_FLAG{MAPPED = 0, REVERSE_COMPLEMENTED = 16};

struct Timer {
	char *topic;
	timeval start, end;
	
	Timer(char *itopic = NULL): topic(itopic) {
		fprintf(stderr, "[start: %s]\n", topic);
		gettimeofday(&start, NULL);
	}

	~Timer() {
		gettimeofday(&end, NULL);
		const float timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
		fprintf(stderr, "[timer: %s] %f\n", topic, (timeuse / 1000000));
	}
};


template <class Value, class Compare>
__device__ __host__
int upper_bound(Value data[], int it_first, int it_last, const Value &value, Compare comp)
{
	int it;
	int count, step;
	count = it_last - it_first;

	while (count > 0)
	{
		it = it_first;
		step = count / 2;
		it = it + step;
		
		if (!comp(value, data[it]))
		{
			it_first = ++it;
			count -= (step + 1);
		}
		else
		{
			count = step;
		}
	}
	
	return it_first;
}


struct INTTYPEComp
{
	__host__ __device__
	bool operator() (const INTTYPE &a, const INTTYPE &b)
	{
		return (a < b);
	}
};


struct PositionChecker
{
	int _realSize;
	INTTYPE *chr_start_pos_key;
	INTTYPE *NPosLen_key;
	INTTYPE *NPosLen_val;
	int len_chr_start_pos;
	int len_NPosLen;

	__host__ __device__
	bool is_vaild(int len_read, int &lowerIter, int &position, bool &isRC, int strand)
	{
		if (position >= _realSize)
		{
			isRC = true;
			position = _realSize * 2 - position - len_read;
		}
		else
		{
			isRC = false;
		}

		if ((strand == 0 && !isRC) || (strand == 1 && isRC))
			return false;

		lowerIter = upper_bound(chr_start_pos_key, 0, len_chr_start_pos, (INTTYPE)position, INTTYPEComp());
		if (lowerIter > 0) lowerIter = lowerIter - 1;
		
		int lowerIter3 = upper_bound(chr_start_pos_key, 0, len_chr_start_pos, (INTTYPE)(position + len_read - 1), INTTYPEComp());
		if (lowerIter3 > 0) lowerIter3 = lowerIter3 - 1;

		if (lowerIter != lowerIter3)
			return false;

		int NLowerIter = upper_bound(NPosLen_key, 0, len_NPosLen, (INTTYPE)position, INTTYPEComp());
		if (NLowerIter > 0) NLowerIter = NLowerIter - 1;

		int NLowerIter3 = upper_bound(NPosLen_key, 0, len_NPosLen, (INTTYPE)(position + len_read - 1), INTTYPEComp());
		if (NLowerIter3 > 0) NLowerIter3 = NLowerIter3 - 1;

		if (NLowerIter != NLowerIter3)
			return false;

		position = position - chr_start_pos_key[lowerIter] + NPosLen_val[NLowerIter];

		return true;
	}
};


template <class T>
void archive_save(std::string filename, T &data) {
	typedef boost::archive::binary_oarchive oarchive;

	std::ofstream outf(filename.c_str());
	if (!outf.is_open()) {
		std::cerr << "oArchive file open failed\n";
	}

	oarchive oa(outf);
	oa << data;
	
	outf.close();
}

template <class T>
void archive_load(std::string filename, T &data) {
	typedef boost::archive::binary_iarchive iarchive;

	std::ifstream inf(filename.c_str(), std::ios::binary);
	if (!inf.is_open()) {
		std::cerr << "iArchive file open failed\n";
	}

	iarchive ia(inf);
	ia >> data;

	inf.close();
}


inline __host__ __device__ int get_token_id(char *token) {
	switch (*token) {
		case 'A':
			return A;
		case 'C':
			return C;
		case 'G':
			return G;
		case 'T':
			return T;
		case '$': //因為read都沒有包含$所以這邊就不管回傳什麼了
			return DOLLAR;
		default: return -1;
	}
}

inline __host__ __device__ Type  my_atoi(const char *str) {
	Type ret;
	unsigned char tmp = 0;

	if (str[0] == 'C') {
		tmp = 1;
	} else if (str[0] == 'G') {
		tmp = 2;
	} else if (str[0] == 'T') {
		tmp = 3;
	} else {
		tmp = 0;
	}

	ret = tmp;

	for (int i = 1; i < ENCODELEN; i++) {
		ret = ret << ENCODEBITS;

		if (str[i] == 'C') {
			tmp = 1;
		} else if (str[i] == 'G') {
			tmp = 2;
		} else if (str[i] == 'T') {
			tmp = 3;
        } else if ( str[i] == '\0' ) { // this condition change the spec
            tmp = 0;
            break;
		} else {
			tmp = 0;
		}

		ret = ret | tmp;
	}

	return ret;
}

inline __host__ __device__ char ch_encoding(const char ch) {
	switch(ch) {
		case 'A':
			return '0';
		case 'C':
			return '1';
		case 'G':
			return '2';
		case 'T':
			return '3';
		default:
			printf("So wrong %c!!!\n", ch);
			return '0';
	}
}

__global__ void seq_encoding(char *dna_seq, const int max_task) {
	int globalID = blockIdx.x * blockDim.x + threadIdx.x;
	if (globalID >= max_task) return;

	char *seq = dna_seq + (globalID * LEN_HANDLE_CHAR);
	for (int i = 0; i < LEN_HANDLE_CHAR; i++) {
		seq[i] = ch_encoding(seq[i]);
	}
}

__global__ void bseq_encoding(char *dna_seq, Type *bseq_ptr, const int bseq_size, const size_t segment_offset) {
	const int startID = (blockIdx.y*gridDim.x+blockIdx.x)*(blockDim.x*blockDim.y);
	int globalID = segment_offset + startID + (threadIdx.y*blockDim.x+threadIdx.x);
	if (globalID >= bseq_size) return;

	bseq_ptr[globalID] = my_atoi(dna_seq + (globalID * ENCODELEN));
	printf("gid %d is %u\n", globalID, bseq_ptr[globalID]);
	printf("gid %d offset %d\n", globalID, globalID * ENCODELEN);
}
__global__ void bseq_encoding_check_seq (char *dna_seq, Type *bseq_ptr, const int64_t seq_size, const int64_t segment_offset) {
	const int64_t startID = (blockIdx.y*gridDim.x+blockIdx.x)*(blockDim.x*blockDim.y);
	int64_t globalID = startID + (threadIdx.y*blockDim.x+threadIdx.x);
	if (globalID >= seq_size) return;
    char* const tmp = dna_seq + (globalID * ENCODELEN);
	bseq_ptr[globalID + segment_offset] = my_atoi ( tmp );
	// printf("gid %d is %u\n", globalID, bseq_ptr[globalID + segment_offset]);
	// printf("gid %d offset %d\n", globalID, globalID * ENCODELEN);
}

inline __host__ __device__ int get_occ_table(int position, int token, int occ_table_reduce[][4], char *sbwt_string) {
	int last_record = position / TBLSTRIDE;
	//char tail_char;
	int token_counter;

	token_counter = occ_table_reduce[last_record][token]; 

	char key_token;

	if (token == 0) key_token = 'A';
	else if (token == 1) key_token = 'C';
	else if (token == 2) key_token = 'G';
	else if (token == 3) key_token = 'T';
	else printf("Wrong key token!!!\n");

	last_record = last_record * TBLSTRIDE;
	for (int i = last_record + 1; i <= position; i++) {
		//tail_char = sbwt_string[i - 1];
		//if (tail_char == key_token) token_counter++;
		if (sbwt_string[i - 1] == key_token) token_counter++;
	}

	return token_counter;
}


inline __host__ __device__ int match_reduce(int token, int position, int c_table[], int occ_table_reduce[][4], char *sbwt_string) {
	return c_table[token] + get_occ_table(position, token, occ_table_reduce, sbwt_string);
}

inline __host__ __device__ int binary_search(int search, int data[], int len_data) {
	int low = 0;
	int high = len_data - 1;
	int mid;

	while (low <= high) {
		mid = (low + high) / 2;

		if (data[mid] == search) {
			return mid;
		} else if (search < data[mid]) {
			high = mid - 1;
		} else if (search > data[mid]) {
			low = mid + 1;
		}	
	}

	return -1;
}


inline void __host__ __device__ find_possible_nr3(char *read, int upper, int bottom,  int c_table[], int occ_table_reduce[][4], int len_location_table, int location_table_key[], int location_table_val[], char *sbwt_string, int maxsize, int &resultcnt, const int &resultbase, int resultz[], bool result_rc[], int result_it[], char fbwt_loc_mark[], PositionChecker *poschk, int strand_opt, const int len_read) {
	int his_upper[TBLSTRIDE];
	int his_bottom[TBLSTRIDE];
	int his_token[TBLSTRIDE];

	const int r = TBLSTRIDE - 1;

	his_upper[0] = upper;
	his_bottom[0] = bottom;
	his_token[1] = 0;

	int ret;
	int lowerIter = 0;
	int position;
	bool isRC;
	for (int level = 1; level <= r; ) {
		while (his_token[level] <= T) {
			if (resultcnt == maxsize) return;
			his_upper[level] = match_reduce(his_token[level], his_upper[level-1], c_table, occ_table_reduce, sbwt_string);
			his_bottom[level] = match_reduce(his_token[level], his_bottom[level-1], c_table, occ_table_reduce, sbwt_string);

			bool go_down = true;

			if (his_upper[level] >= his_bottom[level]) go_down = false;
			else {
				for (int x = his_upper[level]; x < his_bottom[level]; x++) {
					if (fbwt_loc_mark[x] == 0) continue;
					if ((ret = binary_search(x, location_table_key, len_location_table)) != -1) {
						position = location_table_val[ret] + level;
						if (poschk->is_vaild(len_read, lowerIter, position, isRC, strand_opt))
						{
							const int idx_result = resultbase + resultcnt;
							//resultz[idx_result] = location_table_val[ret];
							resultz[idx_result] = position;
							result_rc[idx_result] = isRC;
							result_it[idx_result] = lowerIter;
							resultcnt++;
						}
						if (his_bottom[level] - his_upper[level] == 1) go_down = false;
					}
				}
			}

			if (go_down) { // down
				if (level == r) { // end-boundary continue
					his_token[level]++;
					continue;
				}
				his_token[++level] = 0;
				break;
			} else { // same level continue
				his_token[level]++;
			}
		}
		if (resultcnt >= maxsize) return;
		// i >= np
		// next level start
		if (his_token[level] == 0) continue;
		// up
		if (level == 1) return;
		his_token[--level]++;
	}
}

__global__ void find_read(
		char *reads,
		int c_table[],
		int occ_table_reduce[][4],
		int len_occ_table_reduce,
		char *sbwt_string,
		int numSuffix,
		int len_location_table,
		int location_table_key[],
		int location_table_val[],
		int result[],
		int resultz[],
		bool result_rc[],
		int result_it[],
		int max_threads,
		char fbwt_loc_mark[],
		int read_length,
		PositionChecker *poschk,
		int strand_opt
) {
	int globalID = blockIdx.x * blockDim.x + threadIdx.x;
	if(globalID >= max_threads) return;

	char *read = (reads + (globalID * read_length));
	//printf("findqq: %s\n", read);

	const int len_read = read_length - 1;

	int token = get_token_id(&read[len_read - 1]);

	int upper = c_table[token];
	int bottom = c_table[token + 1];

	for (int i = len_read - 2; i >= 0; i--) {
		token = get_token_id(&read[i]);

		upper = match_reduce(token, upper, c_table, occ_table_reduce, sbwt_string);
		bottom = match_reduce(token, bottom, c_table, occ_table_reduce, sbwt_string);
		if (upper >= bottom) {
			result[globalID] = 0;
			return;
		}
	}

	if ((bottom - upper) > OVERNUM) {
		result[globalID] = 0;
		return;
	}

	result[globalID] = bottom - upper;

	int ret;
	int resultcnt = 0;
	int resultbase = globalID * OVERNUM;
	int lowerIter = 0;
	int position;
	bool isRC;
	for (int i = upper; i < bottom; i++) {
		if (fbwt_loc_mark[i] == 0) continue;
		if ((ret = binary_search(i, location_table_key, len_location_table)) != -1) {
			position = location_table_val[ret];
			if (poschk->is_vaild(len_read, lowerIter, position, isRC, strand_opt))
			{
				const int idx_result = resultbase + resultcnt;
				//resultz[idx_result] = location_table_val[ret];
				resultz[idx_result] = position;
				result_rc[idx_result] = isRC;
				result_it[idx_result] = lowerIter;
		  		resultcnt++;
		  	}
			if (bottom - upper == 1) break;
		}
	}

	find_possible_nr3(read, upper, bottom, c_table, (int (*)[4])occ_table_reduce, len_location_table, location_table_key, location_table_val, sbwt_string, (bottom - upper), resultcnt, resultbase, resultz, result_rc, result_it, fbwt_loc_mark, poschk, strand_opt, len_read);

	result[globalID] = resultcnt;
}


inline __host__ __device__ Type  bseq_segment(const int ref_pos, Type *bseq, const int bseq_size) {
	const int bseq_pos = ref_pos / ENCODELEN;
	const int bseq_offset = ref_pos % ENCODELEN;
	Type s1, s2;

	s1 = bseq[bseq_pos];
	s2 = bseq[bseq_pos + 1];

	const int bit_offset = (bseq_offset * ENCODEBITS); 
	s1 = s1 << bit_offset;
	s2 = s2 >> (TYPEBITS - bit_offset);

	return s1 | s2;
}


__global__ void cpy_in_cuda_bseq(Type *key, Type *bseq, int *idx, int round, int max, const int bseq_size) {
	const int startID = (blockIdx.y*gridDim.x+blockIdx.x)*(blockDim.x*blockDim.y);
	int globalID = startID+(threadIdx.y*blockDim.x+threadIdx.x);
	if(globalID >= max) return;

	const int ref_pos = idx[globalID] + round * ENCODELEN;
	key[globalID] = bseq_segment(ref_pos, bseq, bseq_size);
}


inline void simple_cuda_config(dim3 &dim_grid, dim3 &dim_block, int total, int max_thd = MAX_DIM_THD) {
	dim_block.y = 1;
	dim_grid.y = 1;

	dim_block.x = total;
	dim_grid.x = 1;
	if (dim_block.x > max_thd) {
		dim_grid.x = (dim_block.x + max_thd - 1) / max_thd;
		dim_block.x = max_thd;
	}
}

inline void advance_cuda_config(dim3 &dim_grid, dim3 &dim_block, int total, int max_thd = 16, int max_blk = 256) {
	dim_block.x = max_thd;
	dim_block.y = max_thd;
	
	dim_grid.x = max_blk;
	
	const int count = max_blk * max_thd * max_thd;
	dim_grid.y = (total + count - 1) / count;
}


__global__ void sbwt_chain(int suffix_sorted[],char reference[],  char sbwt_string[], int num_suffix) {
	const int startID = (blockIdx.y*gridDim.x+blockIdx.x)*(blockDim.x*blockDim.y);
	int globalID = startID+(threadIdx.y*blockDim.x+threadIdx.x);

	if (globalID >= num_suffix) return;
	
	if (suffix_sorted[globalID] == 0)
		sbwt_string[globalID] = '$';
	else
		sbwt_string[globalID] = *(reference + suffix_sorted[globalID] - 1);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// api funciton
//
///////////////////////////////////////////////////////////////////////////////////////////////////////
void sbwt_string_create(char *sbwt_string, char *d_reference_char, thrust::host_vector<int> &suffix_array)
{
	dim3 dim3_grid, dim3_block;

	char *d_sbwt_string;
	int num_suffix = suffix_array.size();
	thrust::device_vector<int> d_vals(suffix_array);
	int *d_vals_ptr = thrust::raw_pointer_cast(&d_vals[0]);

	CudaSafeCall( cudaMalloc((void **)&d_sbwt_string, sizeof(char) * num_suffix) );

	advance_cuda_config(dim3_grid, dim3_block, num_suffix);

	sbwt_chain<<<dim3_grid, dim3_block>>>(d_vals_ptr, d_reference_char, d_sbwt_string, num_suffix);
	CudaSafeCall( cudaMemcpy(sbwt_string, d_sbwt_string, sizeof(char) * num_suffix, cudaMemcpyDeviceToHost) );
	cudaFree(d_sbwt_string);
}
// d_reference_char, bseq_ptr和bseq_size都是回傳用的
void compress_reference_check(const char *reference_char, Type **bseq_ptr, int &bseq_size)
{
    const char* ref_end ( reference_char + strlen(reference_char) );
    std::vector<Type> bseq_result;
    std::vector<Type> bseq_answer;
    bseq_result.reserve(bseq_size);
    bseq_answer.reserve(bseq_size);
	CudaSafeCall( cudaMemcpy( bseq_result.begin().base(), *bseq_ptr, bseq_size * sizeof(Type), cudaMemcpyDeviceToHost ) );
	CudaCheckError();
    int k = 0;
    for(const char* itr = reference_char; itr < ref_end; itr += ENCODELEN )
    {
        bseq_answer.push_back(my_atoi(itr));
        k ++;
    }
    // for(size_t i = 0; i < (strlen(reference_char) + ENCODELEN - 1 )/ ENCODELEN; i++)
    // {
    //     bseq_answer[i] = my_atoi(reference_char + (i * ENCODELEN));
    // }
    for(size_t i = 0; i < bseq_size; i ++ )
    {
        assert(bseq_result[i] == bseq_answer[i]);
    }
}
void compress_reference(const char *reference_char, char **d_reference_char, Type **bseq_ptr, int &bseq_size)
{
	const int64_t len_dna_padded = strlen(reference_char);
    std::cout << "len_dna_padded : " << len_dna_padded << std::endl;
    int64_t ref_total_need_mem = sizeof(char) * len_dna_padded + 1;
    std::cout << "ref_total_need_mem : " << ref_total_need_mem << std::endl;
    const int64_t mem_limit = getCudaFreeMemSize() * 0.6 ; // 60% for reference uncompressed data and round to 4's mul
    std::cout << "mem_limit : " << mem_limit << std::endl;
    const int64_t mem_limit_in_char =  ( mem_limit / sizeof(char) ) / ENCODELEN * ENCODELEN; 
    std::cout << "mem_limit_in_char : " << mem_limit_in_char << std::endl;
    int64_t aggrigate_char_num = 0;
	bseq_size = ( len_dna_padded + ENCODELEN - 1) / ENCODELEN;
    std::cout << "bseq_size : " << bseq_size << std::endl;
    std::cout << "ENCODELEN : " << ENCODELEN << std::endl;
	static thrust::device_vector<Type> d_bseq(bseq_size + 1);
	d_bseq[bseq_size] = 0;
	*bseq_ptr = thrust::raw_pointer_cast(&d_bseq[0]);
	// CudaSafeCall( cudaMalloc((void **)d_reference_char, sizeof(char) * len_dna_padded) );
	CudaSafeCall( cudaMalloc( (void **)d_reference_char, std::min ( 
        (int64_t)(mem_limit_in_char * sizeof(char))
        , ref_total_need_mem 
    ) ) );
	CudaCheckError();
    size_t aggrigate_ref_offset = 0;
    int compress_segment_count(0);
    for( int64_t remain(len_dna_padded) ; remain > 0; remain -= mem_limit_in_char )
    {
	    // CudaSafeCall( cudaMemcpy(*d_reference_char, reference_char, sizeof(char) * len_dna_padded, cudaMemcpyHostToDevice) );
        std::cout << "aggrigate_ref_offset : " << aggrigate_ref_offset << std::endl;
        std::cout << "copy size : " << std::min ( 
            (int64_t)(mem_limit_in_char * sizeof(char))
            , remain + 1 )  // + 1 for \0
            << std::endl;
        // std::cout << __FILE__ << __LINE__ << reference_char + aggrigate_ref_offset << std::endl;
	    CudaSafeCall( cudaMemcpy( 
            *d_reference_char
            , reference_char + aggrigate_ref_offset
            , std::min ( 
                (int64_t)(mem_limit_in_char * sizeof(char))
                , remain + 1 ) // + 1 for \0
            , cudaMemcpyHostToDevice) );

	    // bseq is d_reference_char_binary

	    dim3 dim3_grid, dim3_block;

	    // advance_cuda_config(dim3_grid, dim3_block, mem_limit_in_char);
        // size_t bseq_seg_size ( ( mem_limit_in_char + ENCODELEN - 1) / ENCODELEN );
        size_t bseq_seg_size ( 
            std::min(
                mem_limit_in_char / ENCODELEN
                , ( remain + ENCODELEN - 1 )/ ENCODELEN));
        size_t seq_seg_size ( std::min(mem_limit_in_char, remain) );
        std::cout << "seq_seg_size : " << seq_seg_size << std::endl;
        std::cout << "bseq_seg_size : " << bseq_seg_size << std::endl;
	    advance_cuda_config(dim3_grid, dim3_block, bseq_seg_size);
        std::cout << "aggrigate_char_num : " << aggrigate_char_num << std::endl;
	    bseq_encoding_check_seq<<<dim3_grid, dim3_block>>> ( 
            *d_reference_char
            , *bseq_ptr
            // , seq_seg_size
            , bseq_seg_size
            , aggrigate_char_num
        );
        // release d_reference_char
	    CudaCheckError();

        aggrigate_char_num += bseq_seg_size;
        aggrigate_ref_offset += mem_limit_in_char;
        compress_segment_count ++;
    }
    std::cout << "aggrigate_char_num : " << aggrigate_char_num << std::endl;
    std::cout << "compress_segment_count : " << compress_segment_count << '\n';
    CudaSafeCall(cudaFree( *d_reference_char ));
}


// suffix_array 要先放入要排的suffix
thrust::host_vector<int> suffix_sort(int sort_length, thrust::host_vector<int> &suffix_array, Type *bseq_ptr, int bseq_size)
{
	dim3 dim3_grid, dim3_block;

	const int num_suffix = suffix_array.size();

	thrust::host_vector<int> &h_vals = suffix_array;
	{
		thrust::device_vector<Type> d_keys(num_suffix);
		thrust::device_vector<int> d_vals;

		if (SUFFIXLEN % ENCODELEN != 0)
			std::cerr << "Warning... mismatch suffixlen\n";

		d_vals = h_vals;

		Type *d_keys_ptr = thrust::raw_pointer_cast(&d_keys[0]);
		int *d_vals_ptr = thrust::raw_pointer_cast(&d_vals[0]);

		advance_cuda_config(dim3_grid, dim3_block, num_suffix);
		//std::cerr << "num of be sorted suffix: " << num_suffix << "\n";

		// 由LSB排到MSB
		const int PART = (sort_length + ENCODELEN - 1)/ ENCODELEN;
		for (int round = PART - 1; round >= 0; round--) {

			// 複製壓縮序列到排序的key陣列
			cpy_in_cuda_bseq<<<dim3_grid, dim3_block>>>(d_keys_ptr, bseq_ptr, d_vals_ptr, round, num_suffix, bseq_size);
			thrust::sort_by_key(d_keys.begin(), d_keys.end(), d_vals.begin(), thrust::less<Type>());
		}

		h_vals = d_vals;
		showCudaUsage();
		CudaCheckError();

	}

	return h_vals;
}

////////////////////////////////////////////////////////////////////////


std::vector<int> make_random_table(int num_suffix, int random_num)
{
	srand(time(NULL));

	std::set<int> random_table;

	for (int i = 0; i < random_num; i++)
	{
		random_table.insert(rand() % num_suffix);
	}

	std::vector<int> random_table_uniq(random_table.begin(), random_table.end());

	return random_table_uniq;
}


void make_split_table(std::vector<int> &split_table, int split_num, const thrust::host_vector<int> &random_table, std::vector<std::string> &archive_name)
{
	// Sampling
	for (int i = 1; i < split_num; i++)
	{
		int idx = random_table.size() / split_num * i;
		split_table.push_back(random_table[idx]);	
	}
	
	for (int i = 0; i < split_num; i++)
	{
		archive_name.push_back( std::string("split_")+boost::lexical_cast<std::string>(i) );
	}
}


__device__
int bseq_cmp(int pos_a, int pos_b, Type *bseq, int bseq_size)
{
	// SUFFIXLEN必須為(TYPEBYTES * ENCODERATIO)的倍數
	const int len_cmp = (SUFFIXLEN + (TYPEBYTES * ENCODERATIO) - 1) / (TYPEBYTES * ENCODERATIO);

	Type a, b;
	for (int i = 0; i < len_cmp; i++)
	{
		a = bseq_segment(pos_a, bseq, bseq_size);
		b = bseq_segment(pos_b, bseq, bseq_size);
		if (a == b) continue;
		else
		{
			if (a > b) return 1;
			else return -1;
		}
	}
	return 0;
}


__device__
int bseq_lower_bound(int search, int data[], int len_data, Type *bseq, int bseq_size) // data is d_sampler_ptr
{
	int first = 0, last = len_data;
	int count = last - first, step;
	int it;

	while (count > 0)
	{
		it = first;
		step = count / 2;
		it += step;
		if (bseq_cmp(search, data[it], bseq, bseq_size) > 0)
		{
			first = ++it;
			count -= step + 1;
		}
		else
		{
			count = step;
		}
	}
	return first;
}


__global__
void classify_seq_tables_cuda(
      int *d_suffixs_ptr
    , int size_suffix
    , int *d_result_ptr
    , int *d_sampler_ptr
    , int size_sampler
    , Type *bseq_ptr
    , int bseq_size
    , const int start_of_this_stage)
{
	const int startID = (blockIdx.y*gridDim.x+blockIdx.x)*(blockDim.x*blockDim.y);
	const int globalID = startID+(threadIdx.y*blockDim.x+threadIdx.x);
	const int accmulated_globalID = globalID + start_of_this_stage;

	if(accmulated_globalID >= size_suffix) return;

	d_result_ptr[globalID] = bseq_lower_bound(d_suffixs_ptr[globalID],  d_sampler_ptr, size_sampler, bseq_ptr, bseq_size);
}


void classify_seq_tables(
      std::vector<int> &split_table
    , Type *bseq_ptr, int bseq_size
    , thrust::host_vector<int> &suffix_array
    // , std::vector< thrust::host_vector<int> > &SeqTables
    , std::vector< std::pair< std::string, uint64_t > >& SeqTables
    , std::vector<std::string> archive_name)
{
    std::vector< std::ofstream* > SeqTablesf;
    SeqTablesf.reserve(SeqTables.size());
    for ( int i = 0; i < SeqTables.size(); i ++ )
    {
        std::pair<std::string, uint64_t>& p = SeqTables[i];
        std::cout << "open file : " << p.first << std::endl;
        SeqTablesf.push_back( new std::ofstream(p.first.c_str()) );
        p.second = 0;
    }
	dim3 dim3_grid, dim3_block;

	const int num_suffix = suffix_array.size();

	thrust::host_vector<int> &h_suffixs = suffix_array;
	thrust::host_vector<int> h_sampler(split_table);


	thrust::device_vector<int> d_sampler(h_sampler);
	int *d_sampler_ptr = thrust::raw_pointer_cast(&d_sampler[0]);
	//showCudaUsage();

	// can not send data over cuda total memory size
	const int64_t MAX_NUM_TO_CUDA = 200000000;

	int64_t start, end = 0;

	for (int i = 0; end < h_suffixs.size(); i++)
	{
		{
			Timer tm("Classify each");
			
			start = i * MAX_NUM_TO_CUDA;
			end = (i + 1) * MAX_NUM_TO_CUDA;
			if (end > h_suffixs.size())
				end = h_suffixs.size();

			advance_cuda_config(dim3_grid, dim3_block, MAX_NUM_TO_CUDA);

			thrust::device_vector<int> d_suffixs((h_suffixs.begin() + start), (h_suffixs.begin() + end));
			//showCudaUsage();
			thrust::device_vector<int> d_result(end - start);
			//showCudaUsage();

			int *d_suffixs_ptr = thrust::raw_pointer_cast(&d_suffixs[0]);
			int *d_result_ptr = thrust::raw_pointer_cast(&d_result[0]);

			{
				Timer in("Classify inner");
				classify_seq_tables_cuda<<<dim3_grid, dim3_block>>>(
                      d_suffixs_ptr
                    , h_suffixs.size()
                    , d_result_ptr
                    , d_sampler_ptr
                    , h_sampler.size()
                    , bseq_ptr
                    , bseq_size
                    , start
                );
			}

			thrust::host_vector<int> h_result = d_result;

			for (int i = 0; i < h_result.size(); i++)
			{
				// SeqTables[h_result[i]].push_back(h_suffixs[i + start]);
				*SeqTablesf[h_result[i]] << h_suffixs[i + start] << '\n';
                // std::cout << h_suffixs[i + start] << std::endl;
                SeqTables[h_result[i]].second ++ ;
			}
		}
	}
    for ( int i = 0; i < SeqTablesf.size(); i ++ )
    {
        SeqTablesf[i]->flush();
        SeqTablesf[i]->close();
        delete SeqTablesf[i];
    }
    
}

void mkq_sort(
      std::vector<std::string> &archive_name
    , thrust::host_vector<int> &suffix_array
    // , std::vector< thrust::host_vector<int> > &SeqTables
    , std::vector< std::pair< std::string, uint64_t > > &SeqTables
    , Type *bseq_ptr
    , int bseq_size
)
{
	for (int i = 0; i < SeqTables.size(); i++)
	{
		//std::cerr << "uuu\n";
		// thrust::host_vector<int> &sub_suffix_array = SeqTables[i];
		thrust::host_vector<int> sub_suffix_array;
        sub_suffix_array.reserve(SeqTables[i].second );
        std::string seq_table_input_line;
        std::ifstream f(SeqTables[i].first.c_str());
        while(std::getline(f, seq_table_input_line))
        {
            sub_suffix_array.push_back(atoi(seq_table_input_line.c_str()));
        }
        f.close(); std::remove(SeqTables[i].first.c_str());
		if (sub_suffix_array.size() != 0) 
		{
			sub_suffix_array = suffix_sort(SUFFIXLEN, sub_suffix_array, bseq_ptr, bseq_size);
		}

		showCudaUsage();
		{
			Timer tt("Split archive");
			// std::vector<int> each_group(SeqTables[i].begin(), SeqTables[i].end());
			std::vector<int> each_group(sub_suffix_array.begin(), sub_suffix_array.end());
			archive_save(archive_name[i], each_group);

			//thrust::host_vector<int> empty;
			//empty.swap(SeqTables[i]);
		}
	}
		
}


void split_sort(
      std::vector<std::string> &archive_name
    , Type *bseq_ptr
    , int64_t bseq_size
    , int64_t num_suffix
    , thrust::host_vector<int> &suffix_array
    , int64_t average_size = 100000000
    , int len_compare = SUFFIXLEN
)
{
	int64_t split_num = num_suffix / average_size;
	if (split_num == 0) split_num = 1;
	std::cerr << "split_num: " << split_num << "\n";
	int64_t random_num = split_num * 4;

	// split_table is used for sampling
	std::vector<int> split_table;
	// std::vector< thrust::host_vector<int> > SeqTables(split_num);
	std::vector< std::pair<std::string, uint64_t> > SeqTables(split_num);
    for( int i = 0; i < split_num; i ++ )
    {
        SeqTables[i].first = "seq_table_" + boost::lexical_cast<std::string>(i);
    }

	std::vector<int> random_table = make_random_table(num_suffix, random_num);
	//std::vector<int> random_sort_w_suffix = suffix_sort(len_compare, random_table, bseq_ptr, bseq_size);
	thrust::host_vector<int> h_random_table(random_table.begin(), random_table.end());
	thrust::host_vector<int> random_sort_w_suffix = suffix_sort(len_compare, h_random_table, bseq_ptr, bseq_size);

	{
		Timer tm("Make split table");
		make_split_table(split_table, split_num, random_sort_w_suffix, archive_name);
	}

	{
		Timer tm("Classify");
		classify_seq_tables(split_table, bseq_ptr, bseq_size, suffix_array, SeqTables, archive_name);
	}

	{
		thrust::host_vector<int> &tmp = suffix_array;
		thrust::host_vector<int> empty;
		empty.swap(suffix_array);
	}


	thrust::host_vector<int> new_suffix_array;
	{
		Timer tm("MKQ sort");
		mkq_sort(archive_name, new_suffix_array, SeqTables, bseq_ptr, bseq_size);
	}

	archive_save("archive_name.archive", archive_name);
}


// input: chrStartPosFile
// output: chr_start_pos
void readChrStartPos (const std::string& chrStartPosFile, std::map<INTTYPE, std::string> &chr_start_pos) {
	std::ifstream in(chrStartPosFile.c_str());
	std::string line, chr;
	INTTYPE startPos = 0;
	while (getline (in, line)) {
		std::stringstream ss(line);
		ss >> chr >> startPos;
		chr_start_pos.insert(std::make_pair (startPos, chr));
	}
}


// input: chrLenFile
// output: _realSize, chr_length
void readChrLen (const std::string& chrLenFile, int &_realSize, std::map<std::string, int> &chr_length) {
	std::ifstream in(chrLenFile.c_str());
	std::string line, chr;
	int length = 0;
	while (getline (in, line)) {
		std::stringstream ss(line);
		ss >> chr >> length;
		// FIXME: replace _realSize...
		_realSize += length;
		if (chr_length.find(chr) == chr_length.end())
			chr_length.insert(std::make_pair (chr, length));
		else {
			std::cerr << "Error: duplicated chromosome name" << std::endl;
			exit (1);
		}
	}
}


// input: fileName
// output: chr_umbiguous_starting_length
void readNPosLen (const std::string& fileName, std::map<INTTYPE, INTTYPE> &chr_umbiguous_starting_length) {
	std::ifstream fp(fileName.c_str(), std::ios::binary);
	boost::archive::binary_iarchive archive_fp(fp);
	archive_fp & chr_umbiguous_starting_length;
	fp.close();
}


template<class Iter>
void PrintCollection(Iter first, Iter last,
		const char* separator="\n",
		const char* arrow="->",
		const char* optcstr="") 
{
	typedef Iter iter_type;
	std::cerr << optcstr;
	for (iter_type begin = first, it = begin, end = last;
			it != end; ++it) {
		if (it != begin) {
			std::cerr << separator;
		}
		std::cerr << it->first << arrow << it->second;
	}
	std::cerr << std::endl;
}


template <class FileType, class StrType>
void fastq_reader(FileType &reads_fs, StrType &read_name, StrType &read_body, StrType &read_opt, StrType &read_quality)
{
	reads_fs >> read_name;
	reads_fs >> read_body;
	reads_fs >> read_opt;
	reads_fs >> read_quality;
}


std::string rc_seq(char *pseq, int len_seq)
{
	std::string seq(pseq, pseq + len_seq);
	std::reverse(seq.begin(), seq.end());
	for (int i = 0; i < seq.length(); i++)
	{
		if (seq[i] == 'A') seq[i] = 'T';
		else if (seq[i] == 'C') seq[i] = 'G';
		else if (seq[i] == 'G') seq[i] = 'C';
		else if (seq[i] == 'T') seq[i] = 'A';
	}

	return seq;
}


inline std::string r_seq(std::string seq)
{
	std::reverse(seq.begin(), seq.end());

	return seq;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// test funciton
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

__global__
void test_bseq(Type *bseq_ptr, int bseq_size)
{
	for (int i = 0; i < 17; i++)
	{
		printf("test bseq %d %u\n", i, bseq_segment(i, bseq_ptr, bseq_size));
	}
}

__global__
void test_poschk(PositionChecker *poschk)
{
	printf("==[ test poschk ]==\n");
	printf("realSize = %d\n", poschk->_realSize);
	printf("len_chr_start_pos %d\n", poschk->len_chr_start_pos);
	printf("len_NPosLen %d\n", poschk->len_NPosLen);
	for (int i = 0; i < poschk->len_chr_start_pos; i++)
	{
		printf("%lld -> (not needed in cuda anymore)\n", poschk->chr_start_pos_key[i]);
	}

	for (int i = 0; i < poschk->len_NPosLen; i++)
	{
		printf("%lld -> %lld\n", poschk->NPosLen_key[i], poschk->NPosLen_val[i]);
	}
}

void test_split_sort(thrust::host_vector<int> &suffix_array_qq, Type *bseq_ptr, int bseq_size)
{
	const int num_suffix = suffix_array_qq.size();

	std::cerr << "==[ test split sort start]==\n";

	std::vector<std::string> archive_name;
	split_sort(archive_name, bseq_ptr, bseq_size, num_suffix, suffix_array_qq, 40000000);

	std::vector<int> suffix_array;

	for (int i = 0; i < archive_name.size(); i++)
	{
		std::vector<int> sub_suffix_array;
		archive_load(archive_name[i], sub_suffix_array);
		suffix_array.insert(suffix_array.end(), sub_suffix_array.begin(), sub_suffix_array.end());
	}

	showCudaUsage();
	std::cerr << "wwww\n";


	thrust::host_vector<int> suffix_array_golden(num_suffix);
	int from_back = num_suffix - 1;
	for (int i = 0; i < num_suffix; i++) {
		suffix_array_golden[i] = from_back--;
	}
	
	showCudaUsage();
	std::cerr << "wwww\n";

	suffix_array_golden = suffix_sort(SUFFIXLEN, suffix_array_golden, bseq_ptr, bseq_size);

	std::cerr << "suffix_array size " << suffix_array.size() << "\n";
	std::cerr << "suffix_array_golden size " << suffix_array_golden.size() << "\n";

	for (int i = 0; i < suffix_array.size(); i++)
	{
		if (suffix_array[i] != suffix_array_golden[i])
		{
		  std::cerr << i << " " << suffix_array[i] << " " << suffix_array_golden[i] << "\n";
		}
	}
	
	std::cerr << "==[ test split sort finish]==\n";
}


void test_load_archive_my_genome_pre_handler
(
	std::string charStartPosFile, 
	std::string chrLenFile,
	std::string startingLenFile
)
{
	std::map<INTTYPE, std::string> chr_start_pos;
	std::map<std::string, int> chr_length;
	std::map<INTTYPE, INTTYPE> chr_umbiguous_starting_length;
	int _realSize = 0;

	readChrStartPos(charStartPosFile, chr_start_pos);
	readChrLen(chrLenFile, _realSize, chr_length);
	readNPosLen(startingLenFile, chr_umbiguous_starting_length);

	std::cerr << "===[ ChrStartPos ]===\n";
	PrintCollection(chr_start_pos.begin(), chr_start_pos.end());
	std::cerr << "===[ Realsize] ===\n";
	std::cerr << _realSize << "\n";
	std::cerr << "===[ ChrLen ]===\n";
	PrintCollection(chr_length.begin(), chr_length.end());
	std::cerr << "===[ NPosLen ]===\n";
	PrintCollection(chr_umbiguous_starting_length.begin(), chr_umbiguous_starting_length.end());
}



////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// main funciton
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

int STRAND_OPT = FIND_BOTH_STRAND;

void build(int argc, char *argv[])
{
#ifdef IS_CREATE_TBL
	std::string non_reference;
	{
		Timer tm("Pre genome handler");
		my_genome_pre_handler(argv[1], non_reference);
		//std::cerr << non_reference << "\n";
		//test_load_archive_my_genome_pre_handler("my_gph_chrStart", "my_gph_chrLen", "my_gph_NposLen.z");
	}

	dim3 dim3_grid;
	dim3 dim3_block;

	std::string dna;
	dna.swap(non_reference);

	dna.resize(dna.size() - 1);

	char *dna_seq = (char *)dna.c_str();
	const int64_t len_dna_padded = strlen(dna_seq) + SUFFIXLEN;

	const size_t numSuffix = len_dna_padded - (SUFFIXLEN - 1);
	std::cerr << "numSuffix " << numSuffix << "\n";



	char *d_reference_char = NULL;
	Type *bseq_ptr = NULL;
	int bseq_size = 0;

	{
	Timer tm("Table total");
	dna.append(SUFFIXLEN, 'A');

	
	{
		Timer tm("Compress");
		compress_reference(dna.c_str(), &d_reference_char, &bseq_ptr, bseq_size);
        // compress_reference_check(dna.c_str(), &bseq_ptr, bseq_size);
		// test_bseq<<<1, 1>>>(bseq_ptr, bseq_size);
        // exit(1);
	}

	{
		// delete dna
		std::string &tmp = dna;
		std::string tmp2;
		tmp2.swap(tmp);
	}

	thrust::host_vector<int> suffix_array(numSuffix);
	int from_back = numSuffix - 1;
	for (int i = 0; i < numSuffix; i++) {
		suffix_array[i] = from_back--;
	}

//////////////////////////// [ Split ] ///////////////////////////////////////

	//test_split_sort(suffix_array, bseq_ptr, bseq_size);
	//return 0;

	std::vector<std::string> archive_name;
	{
		Timer tm("Split sort total");
		split_sort(archive_name, bseq_ptr, bseq_size, numSuffix, suffix_array, 200000000);
	}


//////////////////////////// [ Split ] ///////////////////////////////////////

	{
		const int len_c_table = 5;
		int *c_table = new int[len_c_table];
		for (int i = 0; i < len_c_table; i++) {
			c_table[i] = 0;
		}
		int len_occ_table_reduce = (numSuffix + TBLSTRIDE -1) / TBLSTRIDE  + 1; // 1 is for ceil
		int *occ_table_reduce = new int[len_occ_table_reduce * 4]; // 4 for 4 char
		int len_location_table = (numSuffix + TBLSTRIDE - 1) / TBLSTRIDE + 1; // ceil and last
		int *location_table_key = new int[len_location_table];
		int *location_table_val = new int[len_location_table];
		int location_cnt = 0;
		int occ_cnt = 0;
		std::vector<char> fbwt_loc_mark_tmp(numSuffix, 0);

		std::string sbwt_string;
		int pitch_accumlate = 0;;
		for (int i = 0; i < archive_name.size(); i++)
		{
			std::vector<int> each_split;
			archive_load(archive_name[i], each_split);
			thrust::host_vector<int> sub_suffix(each_split.begin(), each_split.end());
			char *sub_sbwt = new char[each_split.size() + 1];
			sub_sbwt[each_split.size()] = 0;

			sbwt_string_create(sub_sbwt, d_reference_char, sub_suffix);
			sbwt_string.append(sub_sbwt);

			free(sub_sbwt);

			int pitch_base = pitch_accumlate;
			pitch_accumlate += each_split.size();

			for (int i = 0; i < each_split.size(); i++) {
				if (sub_suffix[i] % TBLSTRIDE == 0) {
					location_table_key[location_cnt] = i + pitch_base;  // key is SBWT index
					location_table_val[location_cnt] = sub_suffix[i]; // val is reference index
					location_cnt++;
					fbwt_loc_mark_tmp[i + pitch_base] = 1;
				}
			}

			for (int i = 0; i < each_split.size(); i++) {
				if (sub_suffix[i] % TBLSTRIDE == 0) {
					fbwt_loc_mark_tmp[i + pitch_base] = 1;
				}
			}


		}

		int xa = 0, xc = 0, xg = 0, xt = 0;
		for (int i = 0; i < numSuffix; i++) {

			if (i % TBLSTRIDE == 0) {
				const int occ_base = occ_cnt * 4;
				occ_table_reduce[occ_base + A] = xa;
				occ_table_reduce[occ_base + C] = xc;
				occ_table_reduce[occ_base + G] = xg;
				occ_table_reduce[occ_base + T] = xt;
				occ_cnt++;
			}

			if (sbwt_string[i] == 'A')
				xa++;
			else if (sbwt_string[i] == 'C')
				xc++;
			else if (sbwt_string[i] == 'G')
				xg++;
			else if (sbwt_string[i] == 'T')
				xt++;

		}

		len_location_table = location_cnt;

		c_table[A] = 1;
		c_table[C] = c_table[A] + xa;
		c_table[G] = c_table[C] + xc;
		c_table[T] = c_table[G] + xg;
		c_table[4] = numSuffix;
	
		///////////////////// Archive save ////////////////////////

		archive_save("sbwt.archive", sbwt_string);

		std::vector<int> occ_table_reduce_v(occ_table_reduce, occ_table_reduce + len_occ_table_reduce * 4);
		archive_save("occ_table.archive", occ_table_reduce_v);

		std::vector<int> c_table_tmp(c_table, c_table + len_c_table);
		archive_save("c_table.archive", c_table_tmp);

		std::vector<int> loc_tbl_key_tmp(location_table_key, location_table_key + len_location_table);
		archive_save("loc_tbl_k.archive", loc_tbl_key_tmp);

		std::vector<int> loc_tbl_val_tmp(location_table_val, location_table_val + len_location_table);
		archive_save("loc_tbl_v.archive", loc_tbl_val_tmp);

		archive_save("fbwt_loc_mark.archive", fbwt_loc_mark_tmp);

		//free(sbwt_string);
		free(occ_table_reduce);
		free(c_table);
		free(location_table_key);
		free(location_table_val);
	}

	}
#endif
}

void map(int argc, char *argv[])
{
#ifdef IS_FIND_READS
	{
		dim3 dim3_grid;
		dim3 dim3_block;

		Timer tmr("Search total");
	///////////////////// Pre genome handler load ////////////////////////

		std::map<INTTYPE, std::string> chr_start_pos;
		std::map<std::string, int> chr_length;
		std::map<INTTYPE, INTTYPE> NPosLen;
		int _realSize = 0;

		readChrStartPos("my_gph_chrStart", chr_start_pos);
		readChrLen("my_gph_chrLen", _realSize, chr_length);
		readNPosLen("my_gph_NposLen.z", NPosLen);

		PositionChecker poschk;
		poschk._realSize = _realSize;
		poschk.len_chr_start_pos = chr_start_pos.size();
		poschk.len_NPosLen = NPosLen.size();
		thrust::device_vector<INTTYPE> dv_chr_start_pos_key(poschk.len_chr_start_pos);
		thrust::device_vector<INTTYPE> dv_NPosLen_key(poschk.len_NPosLen);
		thrust::device_vector<INTTYPE> dv_NPosLen_val(poschk.len_NPosLen);

		std::map<INTTYPE, std::string>::iterator it_chr_start_pos = chr_start_pos.begin();
		for (int i = 0; i < poschk.len_chr_start_pos; i++)
		{
			dv_chr_start_pos_key[i] = it_chr_start_pos->first;
			it_chr_start_pos++;
		}

		std::map<INTTYPE, INTTYPE>::iterator it_NPosLen = NPosLen.begin();
		for (int i = 0; i < poschk.len_NPosLen; i++)
		{
			dv_NPosLen_key[i] = it_NPosLen->first;
			dv_NPosLen_val[i] = it_NPosLen->second;
			it_NPosLen++;
		}

		poschk.chr_start_pos_key = thrust::raw_pointer_cast(dv_chr_start_pos_key.data());
		poschk.NPosLen_key = thrust::raw_pointer_cast(dv_NPosLen_key.data());
		poschk.NPosLen_val = thrust::raw_pointer_cast(dv_NPosLen_val.data());

		PositionChecker *d_poschk;
		CudaSafeCall( cudaMalloc((void **)&d_poschk, sizeof(PositionChecker)) );
		CudaSafeCall( cudaMemcpy(d_poschk, &poschk, sizeof(PositionChecker), cudaMemcpyHostToDevice) );

		//test_poschk<<<1, 1>>>(d_poschk);
		//return 0;

	///////////////////// Archive loc ////////////////////////
		std::string sbwt_arc;
		archive_load("sbwt.archive", sbwt_arc);
		char *p_sbwt_arc = (char *)sbwt_arc.c_str(); // dangerous from const to non-const
		const int sz_sbwt = sbwt_arc.size();
		const int num_suffix = sz_sbwt + 1;

		std::vector<int> occ_table_arc;
		archive_load("occ_table.archive", occ_table_arc);
		int *p_occ_table_arc = &occ_table_arc[0];
		const int sz_occ_table = occ_table_arc.size();

		std::vector<int> c_table_arc;
		archive_load("c_table.archive", c_table_arc);
		int *p_c_table_arc = &c_table_arc[0];
		const int sz_c_table = c_table_arc.size();

		std::vector<int> loc_tbl_key_arc;
		archive_load("loc_tbl_k.archive", loc_tbl_key_arc);
		int *p_loc_tbl_k_arc = &loc_tbl_key_arc[0];
		const int sz_loc_tbl_k = loc_tbl_key_arc.size();

		std::vector<int> loc_tbl_val_arc;
		archive_load("loc_tbl_v.archive", loc_tbl_val_arc);
		int *p_loc_tbl_v_arc = &loc_tbl_val_arc[0];
		const int sz_loc_tbl_v = loc_tbl_val_arc.size();

		std::vector<char> fbwt_loc_mark_arc;
		archive_load("fbwt_loc_mark.archive", fbwt_loc_mark_arc);
		char *p_fbwt_loc_mark_arc = &fbwt_loc_mark_arc[0];
		const int sz_fbwt_loc_mark = fbwt_loc_mark_arc.size();

	/////////////////////         Searching      ////////////////////////

		std::fstream reads_fs;
		reads_fs.open(argv[2]);
		if (!reads_fs.is_open()) {
			std::cerr << "Reads file open failed." << std::endl;
			return;
		}

		std::ofstream out_result(argv[6]);

		int total_num_reads = atoi(argv[5]);

		if (total_num_reads < 1) {
			std::cerr << "Load reads number failed." << std::endl;
			return;
		}


		//////// : transfer to CUDA here
		int *d_result;
		int *d_resultz;
		bool *d_result_rc;
		int *d_result_it;
		char *d_reads;
		int *d_c_table;
		int *d_occ_table_reduce;
		char *d_sbwt_string;
		int *d_location_table_key;
		int *d_location_table_val;
		char *d_fbwt_loc_mark;


		CudaSafeCall( cudaMalloc((void **)&d_c_table, sizeof(int) * sz_c_table) );
		CudaSafeCall( cudaMemcpy(d_c_table, p_c_table_arc, sizeof(int) * sz_c_table, cudaMemcpyHostToDevice) );
		CudaSafeCall( cudaMalloc((void **)&d_occ_table_reduce, sizeof(int) * sz_occ_table) );
		CudaSafeCall( cudaMemcpy(d_occ_table_reduce, p_occ_table_arc, sizeof(int) * sz_occ_table, cudaMemcpyHostToDevice) );
		CudaSafeCall( cudaMalloc((void **)&d_sbwt_string, sizeof(char) * (sz_sbwt + 1)) );
		CudaSafeCall( cudaMemcpy(d_sbwt_string, p_sbwt_arc, sizeof(char) * (sz_sbwt + 1), cudaMemcpyHostToDevice) );
		CudaSafeCall( cudaMalloc((void **)&d_location_table_key, sizeof(int) * sz_loc_tbl_k) );
		CudaSafeCall( cudaMemcpy(d_location_table_key, p_loc_tbl_k_arc, sizeof(int) * sz_loc_tbl_k, cudaMemcpyHostToDevice) );
		CudaSafeCall( cudaMalloc((void **)&d_location_table_val, sizeof(int) * sz_loc_tbl_v) );
		CudaSafeCall( cudaMemcpy(d_location_table_val, p_loc_tbl_v_arc, sizeof(int) * sz_loc_tbl_v, cudaMemcpyHostToDevice) );
		CudaSafeCall( cudaMalloc((void **)&d_fbwt_loc_mark, sizeof(char) * sz_fbwt_loc_mark) );
		CudaSafeCall( cudaMemcpy(d_fbwt_loc_mark, p_fbwt_loc_mark_arc, sizeof(char) * sz_fbwt_loc_mark, cudaMemcpyHostToDevice) );

		int dim_blk = atoi(argv[3]);
		int dim_thd = atoi(argv[4]);

		//showCudaUsage();

		int num_reads = (dim_blk * dim_thd);

		int round = total_num_reads / num_reads;
		if (round == 0)
		{
			num_reads = total_num_reads;
			round = 1;
		}

		std::string peek_read;

		// the second line is read in fastq format
		reads_fs >> peek_read;
		reads_fs >> peek_read;
		const int read_length = peek_read.length() + 1;

		reads_fs.seekg(0);
		reads_fs.seekp(0);

		std::string read_name, read_body, read_opt, read_quality;

		std::vector<std::string> vec_read_name(num_reads);
		std::vector<std::string> vec_read_quality(num_reads);

		for (int y = 0; y < round; y++) {
			//# showCudaUsage();

			char *reads_array = new char[num_reads * read_length];
			char *reads_array_trans;

			int *result;
			int *resultz;
			bool *result_rc;
			int *result_it;
			{
				//@@Timer tt("Allocation");
				for (int i = 0; i < num_reads; i++) {
					reads_array_trans = (reads_array + i * read_length);
					fastq_reader(reads_fs, read_name, read_body, read_opt, read_quality);
					
					vec_read_name[i] = read_name;
					vec_read_quality[i] = read_quality;
					strcpy(reads_array_trans, read_body.c_str());
				}

				CudaSafeCall( cudaMalloc((void **)&d_reads, sizeof(char) * num_reads * read_length) );
				CudaSafeCall( cudaMemcpy(d_reads, reads_array, sizeof(char) * num_reads * read_length, cudaMemcpyHostToDevice) );
				//showCudaUsage();
				CudaSafeCall( cudaMalloc((void **)&d_result, sizeof(int) * num_reads) ); // result
				//showCudaUsage();
				CudaSafeCall( cudaMalloc((void **)&d_resultz, sizeof(int) * num_reads * OVERNUM) );
				//showCudaUsage();
				CudaSafeCall( cudaMalloc((void **)&d_result_rc, sizeof(bool) * num_reads * OVERNUM) );
				//showCudaUsage();
				CudaSafeCall( cudaMalloc((void **)&d_result_it, sizeof(int) * num_reads * OVERNUM) );
				showCudaUsage();

				result = new int[num_reads];
				resultz = new int[num_reads * OVERNUM];
				result_rc = new bool[num_reads * OVERNUM];
				result_it = new int[num_reads * OVERNUM];
			}
			{
				//@@Timer tt("Find once");
				find_read<<<dim_blk, dim_thd>>>(d_reads, d_c_table, (int (*)[4])d_occ_table_reduce, (sz_occ_table / 4), d_sbwt_string, num_suffix, sz_loc_tbl_k, d_location_table_key, d_location_table_val, d_result, d_resultz , d_result_rc, d_result_it, num_reads, d_fbwt_loc_mark, read_length, d_poschk, STRAND_OPT);

				CudaCheckError();

				CudaSafeCall( cudaMemcpy(result, d_result, sizeof(int) * num_reads, cudaMemcpyDeviceToHost) ); 
				CudaSafeCall( cudaMemcpy(resultz, d_resultz, sizeof(int) * num_reads * OVERNUM, cudaMemcpyDeviceToHost) );
				CudaSafeCall( cudaMemcpy(result_rc, d_result_rc, sizeof(bool) * num_reads * OVERNUM, cudaMemcpyDeviceToHost) );
				CudaSafeCall( cudaMemcpy(result_it, d_result_it, sizeof(int) * num_reads * OVERNUM, cudaMemcpyDeviceToHost) );
			}
			{
				//@@Timer tt("IO");
				int base;
				for (int i = 0; i < num_reads; i++) {
					base = i * OVERNUM;
					for (int j = 0; j < result[i]; j++) {
						//std::cout << resultz[base + j] << ' ' << (reads_array + i * read_length) << '\n';
						//std::cout << result_rc[base + j] << '\n';
						//std::cout << result_it[base + j] << '\n';

						std::map<INTTYPE, std::string>::iterator it = chr_start_pos.begin();
						std::advance(it, result_it[base + j]);
						if (result_rc[base + j ] == false)
						{
							std::cout << vec_read_name[i].erase(0, 1) << "\t"
									  << MAPPED << "\t"
									  << it->second << "\t"
									  << resultz[base + j] + 1 << "\t"
									  << 255 << "\t"
									  << (read_length - 1) << "M" << "\t"
									  << "*\t"
									  << 0 << "\t" << 0 << "\t"
									  << (reads_array + i * read_length) << "\t"
									  << vec_read_quality[i] << "\t"
									  << result[i] << "\n";
						}
						else
						{
							std::cout << vec_read_name[i].erase(0, 1) << "\t"
									  << REVERSE_COMPLEMENTED << "\t"
									  << it->second << "\t"
									  << resultz[base + j] + 1 << "\t"
									  << 255 << "\t"
									  << (read_length - 1) << "M" << "\t"
									  << "*\t"
									  << 0 << "\t" << 0 << "\t"
									  << rc_seq(reads_array + i * read_length, read_length - 1) << "\t"
									  << r_seq(vec_read_quality[i]) << "\t"
									  << result[i] << "\n";
						}

					}
				}
			}
			{
				//@@Timer tt("Deallocate");
				CudaSafeCall( cudaFree(d_result) );
				CudaSafeCall( cudaFree(d_resultz) );
				CudaSafeCall( cudaFree(d_result_rc) );
				CudaSafeCall( cudaFree(d_result_it) );
				CudaSafeCall( cudaFree(d_reads) );
				free(reads_array);
				free(result);
				free(resultz);
				free(result_rc);
				free(result_it);
			}
		}
	}
#endif
}

int main(int argc, char *argv[])
{
	build(argc, argv);
	map(argc, argv);

	return 0;
}


