#include<stdio.h>
#include<stdlib.h>

#define N 20
#define M 3



__global__ void add(int *a, int *b, int *c, int n) {
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if(index < n )
		c[index] = a[index] + b[index];	

}

void random_ints(int *x, int n){
	int i;
	for(i=0;i<n;i++){
		x[i]=rand()%99;
	}
}

void print(int *a, int n){
	int i;
	for(i=0;i<n;i++){
		printf(" %d ",a[i]);
	}
	printf("\n");

}



int main(void) {
	int *a, *b, *c; // host copies of a, b, c
	int *d_a, *d_b, *d_c; // device copies of a, b, c
	int size = N * sizeof(int);

	// Alloc space for device copies of a, b, c
	cudaMalloc((void **)&d_a, size);
	cudaMalloc((void **)&d_b, size);
	cudaMalloc((void **)&d_c, size);

	// Alloc space for host copies of a, b, c and setup input values
	a = (int *)malloc(size); random_ints(a, N);
	b = (int *)malloc(size); random_ints(b, N);
	c = (int *)malloc(size);

	print(a,N);
	print(b,N);
	print(c,N);

	// Copy inputs to device
	cudaMemcpy(d_a, a, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, size, cudaMemcpyHostToDevice);

	// Launch add() kernel on GPU with N blocks
	//add<<<N,1>>>(d_a, d_b, d_c);
	// Launch add() kernel on GPU with N threads
	//add<<<1,N>>>(d_a, d_b, d_c,N);
	add<<<(N + M-1) / M,M>>>(d_a, d_b, d_c, N);


	// Copy result back to host
	cudaMemcpy(c, d_c, size, cudaMemcpyDeviceToHost);

	printf("-------------------------------\n");
	print(c,N);
	// Cleanup
	free(a); free(b); free(c);
	cudaFree(d_a); cudaFree(d_b); cudaFree(d_c);
	return 0;
}
