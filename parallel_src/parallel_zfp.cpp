#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>
#include "zfp.h"
#include "mpi.h"
#include "io_utils.hpp"

unsigned char * zfp_compress_3D(float * array, double tolerance, size_t n1, size_t n2, size_t n3, size_t& out_size){
	int status = 0;    /* return value: 0 = success */
	zfp_type type;     /* array scalar type */
	zfp_field* field;  /* array meta data */
	zfp_stream* zfp;   /* compressed stream */
	void* buffer;      /* storage for compressed stream */
	size_t bufsize;    /* byte size of compressed buffer */
	bitstream* stream; /* bit stream to write to or read from */
	size_t zfpsize;    /* byte size of compressed stream */

	/* allocate meta data for the 3D array a[nz][ny][nx] */
	type = zfp_type_float;
	field = zfp_field_3d(array, type, n3, n2, n1);

	/* allocate meta data for a compressed stream */
	zfp = zfp_stream_open(NULL);

	/* set compression mode and parameters via one of three functions */
	/*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
	/*  zfp_stream_set_precision(zfp, precision); */
	zfp_stream_set_accuracy(zfp, tolerance);

	/* allocate buffer for compressed data */
	bufsize = zfp_stream_maximum_size(zfp, field);
	buffer = malloc(bufsize);

	/* associate bit stream with allocated buffer */
	stream = stream_open(buffer, bufsize);
	zfp_stream_set_bit_stream(zfp, stream);
	zfp_stream_rewind(zfp);

	zfpsize = zfp_compress(zfp, field);
    if (!zfpsize) {
      fprintf(stderr, "compression failed\n");
      status = 1;
    }	

	zfp_field_free(field);
	zfp_stream_close(zfp);
	stream_close(stream);
	out_size = zfpsize;
	return (unsigned char *)buffer;
}

float * zfp_decompress_3D(unsigned char * comp_data, double tolerance, size_t buffer_size, size_t n1, size_t n2, size_t n3){
	int status = 0;    /* return value: 0 = success */
	zfp_type type;     /* array scalar type */
	zfp_field* field;  /* array meta data */
	zfp_stream* zfp;   /* compressed stream */
	void* buffer;      /* storage for compressed stream */
	size_t bufsize;    /* byte size of compressed buffer */
	bitstream* stream; /* bit stream to write to or read from */
	size_t zfpsize;    /* byte size of compressed stream */

	/* allocate meta data for the 3D array a[nz][ny][nx] */
	float * array = (float *) malloc(n1 * n2 * n3 * sizeof(float));
	type = zfp_type_float;
	field = zfp_field_3d(array, type, n3, n2, n1);

	/* allocate meta data for a compressed stream */
	zfp = zfp_stream_open(NULL);

	/* set compression mode and parameters via one of three functions */
	/*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
	/*  zfp_stream_set_precision(zfp, precision); */
	zfp_stream_set_accuracy(zfp, tolerance);

	/* allocate buffer for compressed data */
	bufsize = zfp_stream_maximum_size(zfp, field);
	// buffer = malloc(bufsize);
	buffer = (void *) comp_data;
	bufsize = buffer_size;

	/* associate bit stream with allocated buffer */
	stream = stream_open(buffer, bufsize);
	zfp_stream_set_bit_stream(zfp, stream);
	zfp_stream_rewind(zfp);

    if (!zfp_decompress(zfp, field)) {
      fprintf(stderr, "decompression failed\n");
      status = 1;
    }
	zfp_field_free(field);
	zfp_stream_close(zfp);
	stream_close(stream);
	return array;
}

// USAGE
// mpirun -np 16 parallel sz.config folder_num n3 n2 n1
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
		unsigned char * compressed = zfp_compress_3D(data, abs_eb[i], n1, n2, n3, compressed_size[i]);
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
	std::string zip_filename = folder + "/" + "zfp_compressed_" + std::to_string(rank) + ".dat";
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
	zip_filename = folder + "/" + "zfp_compressed_" + std::to_string(size - 1 - rank) + ".dat";
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
		float * dec_data = zfp_decompress_3D(compressed_output_pos, abs_eb[i], compressed_size[i], n1, n2, n3);
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
	posix_delete(folder + "/" + "zfp_compressed_" + std::to_string(rank) + ".dat");
	if (rank == 0){
		printf("ZFP Finish parallel compressing, total compression ratio = %.4g, compressed size = %.4g GB.\n", 1.0*n1*n2*n3*sizeof(float)*num_vars / total_size, total_size * 1.0/ (1<<30));
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
