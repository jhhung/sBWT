#include <bitset>
#include <iostream>
#include <deque>
#include <string>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <ctime>
#include <cstring>

#include "boost/dynamic_bitset.hpp"
#include "boost/utility/binary.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"

#include "constant_def.hpp"
#include "aligner/aligner.hpp"
#include "file_reader.hpp"

typedef Aligner< 
	Aligner_trait<
		Aligner_types::SBWT_Aligner
		, ParallelTypes::M_T
		, FileReader_impl< 
			Fastq 
			,std::tuple <std::string, std::string, std::string, std::string>
			,SOURCE_TYPE::IFSTREAM_TYPE
		>
	>
> ALIGNER_TYPE;

bool checkIndexIntact (const std::string& prefixName)
{
	if ( boost::filesystem::exists (prefixName + "t_table.bwt") &&
		boost::filesystem::exists (prefixName + "NposLen.z") &&
		boost::filesystem::exists (prefixName + "chrStart") &&
		boost::filesystem::exists (prefixName + "chrLen") )
	{
		return true;
	}
	else
	{
		return false;
	}
}


void build(int argc, char** argv)
{
	std::string usage = R"(

*********************************************************************************
+----+
|sBWT|
+----+
	sBWT uses BWT to perform genomic mapping.
	
    sBWT is freely avaible on github: jhhung.github.com/sBWT

# To generate index files from a fasta file.

>  sbwt build

*********************************************************************************

)";
	
	std::string file_fasta;
	std::string index_prefix;
	int sort_length;
	int interval;
	bool overwrite;
	
	boost::program_options::options_description opts (usage);
	
	try 
	{
		opts.add_options ()
			("help,h", "display this help message and exit")
			("input,i", boost::program_options::value<std::string>(&file_fasta)->required(), "The input fasta file.")
			("prefix,p", boost::program_options::value<std::string>(&index_prefix)->required(), "Prefix of index file to generate.")
			("sort_len,l", boost::program_options::value<int>(&sort_length)->default_value(256), "Sort suffix length")
			("interval,s", boost::program_options::value<int>(&interval)->default_value(64), "Index table density for occ and location table")
			("force,f", boost::program_options::bool_switch(&overwrite)->default_value(false), "Overwrite the existing index files if they already exist.")
		;
		boost::program_options::variables_map vm;
		boost::program_options::store (boost::program_options::parse_command_line(argc, argv, opts), vm);
		boost::program_options::notify(vm);
	} catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << opts << std::endl;
        exit (1);
    } catch (...) {
        std::cerr << "Unknown error!" << std::endl;
        std::cerr << opts << std::endl;
        exit (1);
    }
    
	/* test whether input file exist */
	if (!boost::filesystem::exists (file_fasta)) 
	{
		std::cerr << "Error: Input fasta file " << file_fasta << " does not exist! Please double check. Existing..." << std::endl;
		exit (1);
	}

	if (index_prefix.back () != '.') {
		index_prefix += '.';
	}
	
	/* check whether index already exist */
	if (checkIndexIntact (index_prefix) && !overwrite) 
	{
		std::cerr << "Error: index files already exist. If you want to overwrite them, please run it again with option -f.\nExisting..." << std::endl;
		exit (2);
	}
	
	std::cerr << "Parameters fasta: " << file_fasta << "\n";
	std::cerr << "Parameters prefix: " << index_prefix << "\n";
	std::cerr << "Parameters sort_len: " << sort_length << "\n";
	std::cerr << "Parameters interval: " << interval << std::endl;
	
	/* executing build SBWT */
	std::vector<std::string> filelist ({ file_fasta });
	ALIGNER_TYPE aligner;
	aligner.build(filelist, index_prefix, sort_length, interval);
}
void map(int argc, char** argv)
{
	std::string usage = R"(

*********************************************************************************
+----+
|sBWT|
+----+
	sBWT uses BWT to perform genomic mapping.
	
    sBWT is freely avaible on github: jhhung.github.com/sBWT
  
# To map sequences in a fastq file to against an index.

>  sbwt map
   

*********************************************************************************


)";

	std::string file_fastq;
	std::string file_sam;
	std::string index_prefix;
	int nthread;
	int limitNumber;
	int strand;//0+, 1-, 2+-
	
	boost::program_options::options_description opts {usage};
	try
	{
		opts.add_options ()
            ("help,h", "display this help message and exit")
            ("input,i", boost::program_options::value<std::string>(&file_fastq)->required(), "Input fastq file")
            ("index,p", boost::program_options::value<std::string>(&index_prefix)->required(), "Prefix of the index")
            ("output,o", boost::program_options::value<std::string>(&file_sam)->default_value(std::string{"stdout"}), "Output SAM file, stdout by default ")
            ("thread,n", boost::program_options::value<int>(&nthread)->default_value(1), "Number of thread to use; if the number is larger than the core available, it will be adjusted automatically")
            ("multiple,m", boost::program_options::value<int>(&limitNumber)->default_value(1000), "suppress all alignments if > <int> exist")
            ("strand,s", boost::program_options::value<int>(&strand)->default_value(2), "align 0=forward strand; 1=reverse strand; 2=both strands")
        ;
		boost::program_options::variables_map vm;
		boost::program_options::store (boost::program_options::parse_command_line(argc, argv, opts), vm);
		boost::program_options::notify(vm);
		if (vm.count("help") || argc < 4)	{ std::cerr << opts << std::endl; exit (1); }
	}
	catch (std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		std::cerr << opts << std::endl;
		exit (1);
	} catch (...) {
		std::cerr << "Unknown error!" << std::endl;
		std::cerr << opts << std::endl;
		exit (1);
	}
	
	/** check index **/
	if (index_prefix.back () != '.') 
	{
		index_prefix += '.';
	}
	if (!checkIndexIntact (index_prefix)) 
	{
		std::cerr << "Error: index files appear to be damaged. Please rebuild them.\nExiting..." << std::endl;
		exit (2);
	}
	/** check input fastq **/
	if (!boost::filesystem::exists (file_fastq)) {
		std::cerr << "Error: Input fastq file " << file_fastq << " does not exist! Please double check.\nExiting..." << std::endl;
		exit (1);
	}
	/** check output **/
	std::ostream* out{nullptr};
	if (file_sam == "stdout" || file_sam == "-")
	{
		out = &std::cout;
	}
	else
	{
		out = new std::ofstream {file_sam};
		if (!*out) {
			std::cerr << "Error: cannot creat output file " << file_sam << ".\nPlease double check.\nExiting..." << std::endl;
			exit (1);
		}
	}
	/** check thread **/
	auto nCore = boost::thread::hardware_concurrency();
	if ( nCore != 0 && nthread > nCore)
	{
		std::cerr << "Warning: the number of threads set (" << nthread << ") is larger than the number of cores available (" << nCore << ") in this machine.\nSo reset -n=" << nCore << std::endl;
		nthread = nCore;
	}
	
	std::cerr << "Parameters fastq: " << file_fastq << "\n";
	std::cerr << "Parameters prefix: " << index_prefix << "\n";
	std::cerr << "Parameters output: " << file_sam << "\n";
	std::cerr << "Parameters thread: " << nthread << "\n";
	std::cerr << "Parameters multiple: " << limitNumber << "\n";
	std::cerr << "Parameters strand: " << strand << std::endl;
	
	
	/** execute mapping **/
	GlobalPool.ChangePoolSize(nthread);
	typedef std::tuple <std::string, std::string, std::string, std::string > TUPLETYPE;
	std::vector <std::string> fastq_file_vec({ file_fastq });
	std::vector <uint64_t> fastq_size_vec(0);
	
	ALIGNER_TYPE aligner;
	aligner.load_table(index_prefix);
	int total_sam_number = 0;
	
	//version 2, more fast multiple cpus
	aligner.search(fastq_file_vec, *out, limitNumber, strand);
	total_sam_number = std::get<4>(aligner.reads_map_count_);
	
	if (out != &std::cout)
	{
		static_cast<std::ofstream*>(out)-> close ();
		delete out;
	}
	//
	std::cerr << "Total reads: " << std::get<0>(aligner.reads_map_count_) << "\t" << "Single mapped reads: " << std::get<2>(aligner.reads_map_count_) << "\t";
	std::cerr << "Multiple mapped reads: " << std::get<3>(aligner.reads_map_count_);
	std::cerr << "\t Mappability: " << ((double)std::get<1>(aligner.reads_map_count_) / std::get<0>(aligner.reads_map_count_) )*100 << "%" << std::endl;
	std::cerr << "Total sam results number: " << total_sam_number << std::endl;
}


int main(int argc, char** argv)
{
	std::string usage = R"(

*********************************************************************************

+----+
|sBWT|
+----+
	sBWT uses BWT to perform genomic mapping.
	
    sBWT is freely avaible on github: jhhung.github.com/sBWT

Usage:

 1. building index of the genome
>	sbwt build --help
  
 2. mapping fastq to the index
>	sbwt map --help


*********************************************************************************

)";
	if (argc < 2)
	{
		std::cerr << usage << std::endl;
		exit (1);
	}
	
	if (strcmp (argv[1], "build") == 0)
	{
		build (argc-1, argv+1);
	}
	else if (strcmp (argv[1], "map") == 0)
	{
		map (argc-1, argv+1);
	}
	else
	{
		std::cerr << "Error: unrecognized option " << argv[1] << std::endl;
		std::cerr << usage << std::endl;
		exit (1);
	}
	return 0;
};
