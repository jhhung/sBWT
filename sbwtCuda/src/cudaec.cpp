#define CUDA_ERROR_CHECK
 
#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )
 
inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
    if ( cudaSuccess != err )
    {
        fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }
#endif
 
    return;
}
 
inline void __cudaCheckError( const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
	cudaError err = cudaGetLastError();
	if ( cudaSuccess != err )
	{
		fprintf( stderr, "cudaCheckError() failed at %s:%i : %s\n",
				file, line, cudaGetErrorString( err ) );
		exit( -1 );
	}

	// More careful checking. However, this will affect performance.
	// Comment away if needed.
	err = cudaDeviceSynchronize();
	if( cudaSuccess != err )
	{
		fprintf( stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",
				file, line, cudaGetErrorString( err ) );
		exit( -1 );
	}
#endif

	return;
}

void showCudaUsage() {
	size_t free_byte;
	size_t total_byte;
	cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;
	if ( cudaSuccess != cuda_status ){
		fprintf(stderr, "Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
		exit(1);
	}
	
	double free_db = (double)free_byte ;
	double total_db = (double)total_byte ;
	double used_db = total_db - free_db ;
	//printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
	// used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);
	std::cerr << "GPU memory usage: used = "
		<< used_db/1024.0/1024.0
		<< ", free = "
		<< free_db/1024.0/1024.0
		<< " MB, total = "
		<< total_db/1024.0/1024.0
		<< " MB\n";
}
