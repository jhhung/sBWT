sBWT aligner
=========
sBWT aligner is a Burrows–Wheeler transform (BWT) based fast aligner specialized in parallelized indexing from Next Generation Sequencing data. 

sBWT aligner is released under GPLv2 with additional restriction so that is only applicable to individuals and non-profits and that any for-profit company must purchase a different license.    

##INSTALL
*Only 64 bits systems are able to compile and run sBWT aligner.
    
### Run the binary directly without installation 
Please try the precompiled binaries first, most of the linux systems should be able to run Tailor without any troubles.
```bash
bin/sbwt_linux 		# for linux
```

### Compile from the source code
#### Install the dependencies
- 1.1 Relative recent C++ compiler that support most features of C++11. We recommend [GCC](http://gcc.gnu.org/).
- 1.2 [Boost](http://www.boost.org/users/download/)
- 1.3 [CMake](http://www.cmake.org/)

#### Get the latest version of the software
```bash
git clone git@github.com:jhhung/sBWT.git
```

#### Enter the folder sBWT and:
- Set environmental variable `$BOOST_ROOT` to the directory of boost if CMake cannot find boost automatically;
- Set environmental variable `$CC` and `$CXX` to the gcc/g++ compiler you want to use.	
```bash
cmake .
```

#### Compile the software by typing:
```bash
make
```

#### troubleshooting
- If you got linker error, it is possible that the default library in the lib/ is not suitable to your platform. 
 There are one librarie for linux available, rename the file that fit the best to "libabwt_table.a",
 and retype 
```bash
make
```
	
##USAGE

#### Build genomic index (similar to bowtie-build)
```bash
sbwt build -i [genome.fa] -p [ouput prefix]
```

#### Mapping 
```bash
sbwt map -i [reads.fq] -p [index prefix] -o [output.sam]
```

GPU version
===========

##Install
In sbwtCuda direction
```bash
bin/com_rf	#for mapping
bin/com_tbl	#for indexing
```
### Compile from the source code
#### Install the dependencies
- 1.1 Relative recent C++ compiler that support most features of C++11. We recommend [GCC](http://gcc.gnu.org/).
- 1.2 [Boost](http://www.boost.org/users/download/)
- 1.3 CUDA's Compute Capability 2.0, Tesla K20 is recommended


#### Enter the folder sbwtCuda and:
- Set environmental variable `$CC` and `$CXX` to the gcc/g++ compiler you want to use.
```bash
make com_rf		#for mapping
make com_tbl	#for indexing
```

##USAGE
```bash
com_tbl [reference.fa] xxx 3125 128 4000000					#for indexing
#the bwt file will generate in local diretion where you run the program.

com_rf xxx [reads.fq] 3125 128 [reads number] > [output.sam] 	#for mapping
#the "xxx" and numeric parameter is unessesery please just copy it.
```


##Contact
- Jui-Hung Hung <juihunghung@gmail.com>
- Min-Te Chou <poi5305@gmail.com>
- Ting-Wei Hong <thestyle1202@gmail.com>
- Chia-Hua Chang <CHChang810716@gmail.com>


