#APP=complete.cu
APP=complete
CC=nvcc
CFLAGS=-O3
NVFLAGS=-arch=sm_20
BOOSTFLAGS=-DBOOST_NOINLINE='__attribute__ ((noinline))'

SEFLAGS=-I/usr/include/boost -L/usr/lib -lboost_serialization

#all: $(APP).cu
	#$(CC) -o $(APP) $(APP).cu $(CFLAGS) $(NVFLAGS)
#all: $(APP).cu
#	$(CC) -o $(APP) $(APP).cu $(CFLAGS) $(NVFLAGS)
#com: complete.cu
#	nvcc complete.cu -o complete -arch=sm_20 -O3 $(SEFLAGS) -DIS_CREATE_TBL -DIS_FIND_READS

#com_tbl_32: complete.cu
#	nvcc complete.cu -o com_tbl_32 -arch=sm_20 -O3 $(SEFLAGS) -DIS_CREATE_TBL

#com_tbl_64: complete.cu
#	nvcc complete.cu -o com_tbl_64 -arch=sm_20 -O3 $(SEFLAGS) -DIS_CREATE_TBL

#com_rf_32: complete.cu
#	nvcc complete.cu -o com_rf_32 -arch=sm_20 -O3 $(SEFLAGS) -DIS_FIND_READS

#com_rf_64: complete.cu
#	nvcc complete.cu -o com_rf_64 -arch=sm_20 -O3 $(SEFLAGS) -DIS_FIND_READS

#mix: complete.cu my_genome_pre_handler.hpp
#	nvcc $(BOOSTFLAGS) complete.cu -o complete -arch=sm_20 -O3 -I/usr/local/include/boost -L/usr/local/lib -lboost_serialization -DIS_CREATE_TBL -DIS_FIND_READS my_genome_pre_handler.o -lboost_thread -lboost_system -lboost_regex -lboost_iostreams -lboost_program_options -lboost_filesystem -lboost_chrono -lboost_program_options -lpthread -lrt -lz -lcurl -w

#mix_tbl_32: complete.cu my_genome_pre_handler.hpp
#	nvcc $(BOOSTFLAGS) complete.cu -o com_tbl_32 -arch=sm_20 -O3 -I/usr/local/include/boost -L/usr/local/lib -lboost_serialization -DIS_CREATE_TBL my_genome_pre_handler.o -lboost_thread -lboost_system -lboost_regex -lboost_iostreams -lboost_program_options -lboost_filesystem -lboost_chrono -lboost_program_options -lpthread -lrt -lz -lcurl -w

#mix_rf_32: complete.cu my_genome_pre_handler.hpp
#	nvcc $(BOOSTFLAGS) complete.cu -o com_rf_32 -arch=sm_20 -O3 -I/usr/local/include/boost -L/usr/local/lib -lboost_serialization -DIS_FIND_READS my_genome_pre_handler.o -lboost_thread -lboost_system -lboost_regex -lboost_iostreams -lboost_program_options -lboost_filesystem -lboost_chrono -lboost_program_options -lpthread -lrt -lz -lcurl -w

#qq_tbl_32: complete.cu my_genome_pre_handler.hpp
#	nvcc -DBOOST_NOINLINE='__attribute__ ((noinline))' complete.cu -o com_tbl_32 -arch=sm_20 -O3 -I/usr/local/include/boost -L/usr/local/lib /usr/local/lib/libboost_serialization.a -DIS_CREATE_TBL my_genome_pre_handler.o /usr/local/lib/libboost_thread.a /usr/local/lib/libboost_system.a -lboost_regex -lboost_iostreams -lboost_program_options -lboost_filesystem -lboost_chrono -lpthread -lrt -lz -lcurl -w

com_tbl:
	nvcc -DBOOST_NOINLINE='__attribute__ ((noinline))' src/complete.cu -o bin/com_tbl -arch=sm_20 -O3 -I/usr/local/include/boost -L/usr/local/lib /usr/local/lib/libboost_serialization.a -DIS_CREATE_TBL src/my_genome_pre_handler.o /usr/local/lib/libboost_thread.a /usr/local/lib/libboost_system.a -lboost_regex -lboost_iostreams -lboost_program_options -lboost_filesystem -lboost_chrono -lpthread -lrt -lz -lcurl -w

com_rf: 
	nvcc -DBOOST_NOINLINE='__attribute__ ((noinline))' src/complete.cu -o bin/com_rf -arch=sm_20 -O3 -I/usr/local/include/boost -L/usr/local/lib /usr/local/lib/libboost_serialization.a -DIS_FIND_READS src/my_genome_pre_handler.o /usr/local/lib/libboost_thread.a /usr/local/lib/libboost_system.a -lboost_regex -lboost_iostreams -lboost_program_options -lboost_filesystem -lboost_chrono -lpthread -lrt -lz -lcurl -w
clean:
	rm $(APP)
