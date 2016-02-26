#APP=complete.cu
APP=complete
CC=nvcc
CFLAGS=-O3
NVFLAGS=-arch=sm_20
BOOSTFLAGS=-DBOOST_NOINLINE='__attribute__ ((noinline))'

SEFLAGS=-I/usr/include/boost -L/usr/lib -lboost_serialization

all: com_tbl com_rf

com_tbl:
	nvcc $(BOOSTFLAGS) src/complete.cu -o bin/com_tbl $(NVFLAGS) -O3 -I/usr/local/include/boost -L/usr/local/lib /usr/local/lib/libboost_serialization.a -DIS_CREATE_TBL src/my_genome_pre_handler.o /usr/local/lib/libboost_thread.a /usr/local/lib/libboost_system.a -lboost_regex -lboost_iostreams -lboost_program_options -lboost_filesystem -lboost_chrono -lpthread -lrt -lz -lcurl -w

com_rf: 
	nvcc $(BOOSTFLAGS) src/complete.cu -o bin/com_rf $(NVFLAGS) -O3 -I/usr/local/include/boost -L/usr/local/lib /usr/local/lib/libboost_serialization.a -DIS_FIND_READS src/my_genome_pre_handler.o /usr/local/lib/libboost_thread.a /usr/local/lib/libboost_system.a -lboost_regex -lboost_iostreams -lboost_program_options -lboost_filesystem -lboost_chrono -lpthread -lrt -lz -lcurl -w
clean:
	rm $(APP)
