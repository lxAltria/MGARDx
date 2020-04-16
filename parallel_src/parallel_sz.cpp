#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>
#include "sz.h"
#include "mpi.h"
#include "io_utils.hpp"

int main(int argc, char ** argv)
{

	size_t n1 = 0, n2 = 0, n3 = 0;
	MPI_Init(&argc, &argv);

	int size = 0;
  int rank = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	
	if(argc < 5){
		printf("Test case: mpirun -np 16 ./a.out n1 n2 n3 eb\n");
		exit(0);
	}

	n1 = atoi(argv[1]);
	n2 = atoi(argv[2]);
	n3 = atoi(argv[3]);
	double rel_eb = atof(argv[4]);
	if(rank == 0){
		printf("Data dimension {%d, %d, %d}, relative eb = %.4g\n", n1, n2, n3, rel_eb);
		fflush(stdout);
	}

	int num_vars = 12;
	char file[12][50] ={"PRES-98x1200x1200.dat", "QG-98x1200x1200.dat", "QR-98x1200x1200.dat", "QV-98x1200x1200.dat", "T-98x1200x1200.dat", "V-98x1200x1200.dat", "QC-98x1200x1200.dat", "QI-98x1200x1200.dat", "QS-98x1200x1200.dat", "RH-98x1200x1200.dat", "U-98x1200x1200.dat", "W-98x1200x1200.dat"};

	size_t compressed_size[12] = {0};
	double abs_eb[12] = {0};
	std::string folder("/gpfs/alpine/scratch/jieyang/csc143/xin/parallel/SCALE");
	size_t num_elements = n1 * n2 * n3;
	size_t est_compression_ratio = 20;
	size_t est_compressed_size = n1 * n2 * n3 * sizeof(float) * num_vars / est_compression_ratio;
	unsigned char * compressed_output = (unsigned char *) malloc(est_compressed_size);
	unsigned char * compressed_output_pos = compressed_output;
	double elapsed_time = 0;
	double compression_time = 0, writing_time = 0, reading_time = 0, decompression_time = 0;
  if(rank == 0){
		printf("Start parallel compressing, #process = %d\n", size);
		fflush(stdout);
	}
	for(int i=0; i<num_vars; i++){
		std::string filename = folder + "/" + file[i];
		float * data = (float *) malloc(num_elements * sizeof(float));
		// Read Input Data
		/*
		if(rank == 0){
			elapsed_time = -MPI_Wtime();
			posix_read(filename, data, num_elements);
			elapsed_time += MPI_Wtime();
			printf("file %s read time: %.2f\n", filename.c_str(), elapsed_time);
			MPI_Barrier(MPI_COMM_WORLD);
			elapsed_time = -MPI_Wtime();
			MPI_Bcast(data, num_elements, MPI_FLOAT, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			elapsed_time += MPI_Wtime();
			printf("broadcast time: %.2f\n", elapsed_time);
			fflush(stdout);
		}
		else{
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(data, num_elements, MPI_FLOAT, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
		}
		*/
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank == 0) elapsed_time = -MPI_Wtime();
		posix_read(filename, data, num_elements);
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank == 0){
			elapsed_time += MPI_Wtime();
			printf("file %s read time: %.2f\n", filename.c_str(), elapsed_time);
		}
		// Compute value range
		float max_val = data[0];
		float min_val = data[0];
		for(int i=0; i<num_elements; i++){
			if(max_val < data[i]) max_val = data[i];
			if(min_val > data[i]) min_val = data[i];
		}
		abs_eb[i] = rel_eb * (max_val - min_val);
		// Compress Input Data
		if(rank == 0) printf("Compressing %s with absolute eb %.4g\n", filename.c_str(), abs_eb[i]);
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank == 0) elapsed_time = -MPI_Wtime();
		unsigned char * compressed = SZ_compress_args(SZ_FLOAT, data, &compressed_size[i], ABS, abs_eb[i], 0, 0, 0, 0, n1, n2, n3);
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank == 0){
			elapsed_time += MPI_Wtime();
			compression_time += elapsed_time;
		}
		free(data);
		memcpy(compressed_output_pos, compressed, compressed_size[i]);
		compressed_output_pos += compressed_size[i];
		free(compressed);
	}
	if(rank == 0){
		printf("Compressing using %d processes, elapsed time = %.4g seconds\n", size, compression_time);
		fflush(stdout);
	}
	std::string zip_filename = folder + "/" + "sz_compressed_" + std::to_string(rank) + ".dat";
	size_t total_size = compressed_output_pos - compressed_output;
	// Write Compressed Data
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0) elapsed_time = -MPI_Wtime();
	posix_write(zip_filename, compressed_output, total_size);
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){
		elapsed_time += MPI_Wtime();
		writing_time += elapsed_time;
	}
	free(compressed_output);
  if(rank == 0){
		printf("Writing compressed files, elapsed time = %.4g seconds, writing performance %.4g GB/s\n", writing_time, size * total_size * 1.0/ (1<<30) / writing_time);
    fflush(stdout);
  }
	clear_cache();
	// reverse reading order to avoid reading cache
	zip_filename = folder + "/" + "sz_compressed_" + std::to_string(size - 1 - rank) + ".dat";
	// Read Compressed Data
	compressed_output = (unsigned char *) malloc(total_size);
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0) elapsed_time = -MPI_Wtime();
	posix_read(zip_filename, compressed_output, total_size);
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){
		elapsed_time += MPI_Wtime();
		reading_time += elapsed_time;
	}
  if(rank == 0){
		printf("Reading compressed files, elapsed time = %.4g seconds, reading performance %.4g GB/s\n", reading_time, size * total_size * 1.0/ (1<<30) / reading_time);
    fflush(stdout);
  }
	compressed_output_pos = compressed_output;
	for(int i=0; i<num_vars; i++){
		// Decompress Compressed Data
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank == 0) elapsed_time = -MPI_Wtime();
		float * dec_data =(float *) SZ_decompress(SZ_FLOAT, compressed_output_pos, compressed_size[i], 0, 0, n1, n2, n3);
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank == 0){
			elapsed_time += MPI_Wtime();
			decompression_time += elapsed_time;
		}
		compressed_output_pos += compressed_size[i];
		free(dec_data);
	}
	free(compressed_output);
  if(rank == 0){
    printf("Decompressing using %d processes, elapsed time = %.4g seconds\n", size, decompression_time);
    fflush(stdout);
  }
	posix_delete(folder + "/" + "sz_compressed_" + std::to_string(rank) + ".dat");
	if (rank == 0){
		printf("SZ Finish parallel compressing, total compression ratio = %.4g, compressed size = %.4g GB.\n", 1.0*n1*n2*n3*sizeof(float)*num_vars / total_size, total_size * 1.0/ (1<<30));
		printf("Separate ratios: ");
		for(int i=0; i<num_vars; i++){
			printf("%.4g ", 1.0*n1*n2*n3*sizeof(float) / compressed_size[i]);
		}
		printf("\n");
		//printf("Writing performance = %.4g GB/s, reading performance = %.4g GB/s\n", size * total_size * 1.0/ (1<<30) / writing_time, size * total_size * 1.0/ (1<<30) / reading_time);
		//printf("Reading compressed files, elapsed time = %.4g seconds\n", reading_time);
		//printf("Writing compressed files, elapsed time = %.4g seconds\n", writing_time);
		//printf("Compressing using %d processes, elapsed time = %.4g seconds\n", size, compression_time);
		//printf("Decompressing using %d processes, elapsed time = %.4g seconds\n\n", size, decompression_time);
	}
	MPI_Finalize();

	return 0;
}
