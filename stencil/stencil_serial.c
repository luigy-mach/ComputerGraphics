#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 30
//#define N 100000
//#define BLOCK_SIZE 32
#define BLOCK_SIZE 8
//#define NUM_BLOCKS N+BLOCK_SIZE-1/BLOCK_SIZE 
#define RADIUS 3

//int temp[BLOCK_SIZE + 2 * RADIUS];

void fill_vec(int *a, int n){
	int i;
	for(i=0;i<n;i++){
		a[i] = i+1;
	}
}

void print_vec(int *a, int n){
	int i;
	for(i=0;i<n;i++){
		printf("%d ",a[i]);
	}
	printf("\n");
}

void fill_stencil(int* in,int* out)
{	
	int tamBlock = BLOCK_SIZE + 2 * RADIUS;

	int temp[tamBlock];
	int index = 0;
	int i;
	int j;
	int k;

	memcpy(temp,in,tamBlock+1);

	for(i=index+RADIUS;N+RADIUS;i++)
	{

		for(j=i-RADIUS;j<RADIUS;j++)
		{
			for(k=)
			out[j] +=  in[j];

		}

	}
}

int main()
{

	clock_t start_t, end_t, total_t;

	int* input  ;
	int* output ;
	int  size = N * sizeof(int);

	input 	= (int*) malloc(size);
	output	= (int*) malloc(size+2*RADIUS);

	fill_vec(input,N);
	print_vec(input,N);

	start_t = clock();
	fill_stencil(input,output);
	end_t = clock();

	int zz[BLOCK_SIZE]
	memcpy(temp,in,BLOCK_SIZE+1);
	print_vec(temp,BLOCK_SIZE+1);

	total_t = (end_t - start_t);
	printf("time: %f \n",(double)(total_t)/CLOCKS_PER_SEC);

	return 0;
}