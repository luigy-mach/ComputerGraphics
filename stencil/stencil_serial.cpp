#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

//#define N 30
#define N 100000
#define BLOCK_SIZE 32
//#define BLOCK_SIZE 4
//#define NUM_BLOCKS N+BLOCK_SIZE-1/BLOCK_SIZE 
#define RADIUS 3
//#define RADIUS 2

//int temp[BLOCK_SIZE + 2 * RADIUS];

void fill_vec(int *a, int n){
	int i;
	for(i=0;i<n;i++){
		a[i] = i+1;
	}
}


void fill_vec_value(int *a, int n,int value){
	int i;
	for(i=0;i<n;i++){
		a[i] = value;
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
	int index = 0+RADIUS;
	int i;
	int j;
	int k;


//memcpy(temp,in,sizeof(int)*tamBlock);
	for( i=0 ; i<N ; i+=BLOCK_SIZE+2*RADIUS){
		//printf("aa\n");
		for( j=i+RADIUS ; j<i+BLOCK_SIZE+2*RADIUS ; j++){
			//printf("	bb\n");
			for( k=j-RADIUS ; k<j+RADIUS+1 ;k++)
				//printf("		cc %d \n",out[j]);
				out[j] += in[k];
		}

	}
}

int main()
{

	clock_t start_t, end_t, total_t;

	int* input  ;
	int* output ;
	int  size_in = N * sizeof(int);
	int  size_out = (N + 2*RADIUS) * sizeof(int);

	input 	= (int*) malloc(size_in);
	output	= (int*) malloc(size_out);

	fill_vec(input,N);
	//int tam_out=N+2*RADIUS;
	int tam_out=N+2*RADIUS;
	fill_vec_value(output,tam_out,0);

	printf("in ----------------------------- \n");
	//print_vec(input,N);
	printf("----------------------------- \n");
	printf("out ----------------------------- \n");
	//print_vec(output,tam_out);
	printf("----------------------------- \n");

	start_t = clock();
	fill_stencil(input,output);
	printf("#----------------------------- \n");
	//print_vec(output,tam_out);
	printf("#----------------------------- \n");
	end_t = clock();

	//int zz[BLOCK_SIZE];
	//memcpy(zz,input,sizeof(int)*BLOCK_SIZE);
	//print_vec(zz,BLOCK_SIZE);

	total_t = (end_t - start_t);
	printf("time: %f \n",(double)(total_t)/CLOCKS_PER_SEC);

	return 0;
}