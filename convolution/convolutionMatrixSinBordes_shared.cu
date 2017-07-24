// Luigy Machaca Arcana
// Computer science - Arequipa, Perú  2017


#include <stdlib.h>
#include <stdio.h>

#include <fstream>
#include <iostream>
#include <string>


using namespace std;
#define WIDTH_TILE 32



#define TAM_KERNEL 3
#define CONV_TILE_WIDTH 3
#define CONV_TILE_WIDTH2 CONV_TILE_WIDTH+TAM_KERNEL-1
//#define CONV_BLOCK_WIDTH 12

__global__ void convolution_shared(int** dd_mat_a, int n_rows_a, int n_cols_a ,int** dd_mat_b, int n_rows_b, int n_cols_b, int** dd_mat_c, int n_rows_c, int n_cols_c){


	int n_kernel_row = n_rows_b; //n_cols_b
	int n_kernel_col = n_cols_b; //n_cols_b
	
	//const int row_tile = CONV_TILE_WIDTH+n_kernel_row-1;
	//const int col_tile = CONV_TILE_WIDTH+n_kernel_col-1;
	__shared__ int Ns[CONV_TILE_WIDTH2][CONV_TILE_WIDTH2];

	int bx = blockIdx.x;
	int by = blockIdx.y;

	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int blockDimx = blockDim.x;
	int blockDimy = blockDim.y;


	int row = by*blockDimy + ty;
	int col = bx*blockDimx + tx;	


	if ( (row >= 0) && (row < n_rows_a) && (col >= 0) && (col < n_cols_a) ) 
	{
      Ns[ty][tx] = dd_mat_a[row][col]; 
    }
    else {
      Ns[ty][tx] = 0;
    }

    //int value;

	if( ((int)(n_kernel_row/2)-1)< row && row<(n_rows_a-(int)(n_kernel_row/2)) && 
		((int)(n_kernel_col/2)-1)< col && col<(n_cols_a-(int)(n_kernel_col/2)) 	){

		double offset = 0;
		for(int k=(n_kernel_row/2) ; k<CONV_TILE_WIDTH+CONV_TILE_WIDTH-(n_kernel_row/2); k++){
			for(int l=(n_kernel_col/2) ; l<CONV_TILE_WIDTH+CONV_TILE_WIDTH-(n_kernel_col/2) ; l++){
				double cc = dd_mat_b[k-n_kernel_row/2][l-n_kernel_col/2];
				double dd = 0;
				//dd = Ns[ty-(int)(n_kernel_row/2)+k][tx-(int)(n_kernel_col/2)+l];
				dd = Ns[k][l];
				offset += cc*dd;
			}
		}
		//dd_mat_c[row][col] = offset;
	}
		dd_mat_c[row][col] = -1;
}




__global__ void convolution(int** dd_mat_a, int n_rows_a, int n_cols_a ,int** dd_mat_b, int n_rows_b, int n_cols_b, int** dd_mat_c, int n_rows_c, int n_cols_c){

	int n_kernel_row = n_rows_b; //n_cols_b
	int n_kernel_col = n_cols_b; //n_cols_b
	
	int bx = blockIdx.x;
	int by = blockIdx.y;

	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int blockDimx = blockDim.x;
	int blockDimy = blockDim.y;


	int row = by*blockDimy + ty;
	int col = bx*blockDimx + tx;	

	if( ((int)(n_kernel_row/2)-1)< row && row<(n_rows_a-(int)(n_kernel_row/2)) && 
		((int)(n_kernel_col/2)-1)< col && col<(n_cols_a-(int)(n_kernel_col/2)) 	){

		double offset = 0;
		for(int k=0 ; k<n_kernel_row ; k++){
			for(int l=0 ; l<n_kernel_col ; l++){
				double cc = dd_mat_b[k][l];
				double dd = 0;
				dd = dd_mat_a[row-(int)(n_kernel_row/2)+k][col-(int)(n_kernel_col/2)+l];
				/*
				if( (i-(int)(p/2)+k)>=0 && (j-(int)(q/2)+l)>=0 &&
					(i-(int)(p/2)+k)<n && (j-(int)(q/2)+l)<m  ){
					dd = in.at( i-(int)(p/2)+k , j-(int)(q/2)+l );
				}
				*/
				offset += cc*dd;
			}
		}
		dd_mat_c[row][col] = offset;
		//dd_mat_c[row][col] = -1;
	}

}


__global__ void convolution_complete(int** dd_mat_a, int n_rows_a, int n_cols_a ,int** dd_mat_b, int n_rows_b, int n_cols_b, int** dd_mat_c, int n_rows_c, int n_cols_c){

	int n_kernel_row = n_rows_b; //n_cols_b
	int n_kernel_col = n_cols_b; //n_cols_b
	
	int bx = blockIdx.x;
	int by = blockIdx.y;

	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int blockDimx = blockDim.x;
	int blockDimy = blockDim.y;


	int row = by*blockDimy + ty;
	int col = bx*blockDimx + tx;	


	if( row<n_rows_a && col<n_cols_a ){

		double offset = 0;
		for(int k=0 ; k<n_kernel_row ; k++){
			for(int l=0 ; l<n_kernel_col ; l++){
				double cc = dd_mat_b[k][l];
				double dd = 0;
				//dd = dd_mat_a[row-(int)(n_kernel_row/2)+k][col-(int)(n_kernel_col/2)+l];
				if( (row-(int)(n_kernel_row/2)+k)>=0  && (row-(int)(n_kernel_row/2)+k)<n_rows_a &&
					(col-(int)(n_kernel_col/2)+l)>=0  && (col-(int)(n_kernel_col/2)+l)<n_cols_a  ){
					dd = dd_mat_a[row-(int)(n_kernel_row/2)+k][col-(int)(n_kernel_col/2)+l];
					
				}
				offset += cc*dd;
			}
		}
		dd_mat_c[row][col] = offset;
		//dd_mat_c[row][col] = -1;
	}

}



__global__ void matrix_mult_shared(int** dd_mat_a, int n_rows_a, int n_cols_a ,int** dd_mat_b, int n_rows_b, int n_cols_b, int** dd_mat_c, int n_rows_c, int n_cols_c){

	
	__shared__ int Mds[WIDTH_TILE][WIDTH_TILE];
	__shared__ int Nds[WIDTH_TILE][WIDTH_TILE];

	int bx = blockIdx.x;
	int by = blockIdx.y;

	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int value = 0;

	int row = by*WIDTH_TILE + ty;
	int col = bx*WIDTH_TILE + tx;	

	int width = n_cols_a; //n_cols_a == n_rows_b

	int k;
	for( k=0 ; k<(int)(width-1+WIDTH_TILE)/(int)WIDTH_TILE ; ++k ){
		if (k*WIDTH_TILE+tx < n_cols_a && row < n_rows_a){
			Mds[ty][tx] = dd_mat_a[row][k*WIDTH_TILE+tx];
		}
        else{
			Mds[ty][tx] = 0;
        }

        if (k*WIDTH_TILE+ty < n_rows_b && col < n_cols_b){
			Nds[ty][tx] = dd_mat_b[k*WIDTH_TILE+ty][col];
        }
        else{
			Nds[ty][tx] = 0;
        }

		__syncthreads();
		int m;
		for(m=0 ; m<WIDTH_TILE ; ++m){
			value += Mds[ty][m]*Nds[m][tx];
		}
		__syncthreads();

	}

	if(row<n_rows_c && col<n_cols_c){
		dd_mat_c[row][col]=value;
	}
	

}


__global__ void matrix_mult(int** dd_mat_a, int n_rows_a, int n_cols_a ,int** dd_mat_b, int n_rows_b, int n_cols_b, int** dd_mat_c, int n_rows_c, int n_cols_c){
	int value=0;


	int tx=threadIdx.x;
	int ty=threadIdx.y;


	int x = tx + blockIdx.x*blockDim.x;
	int y = ty + blockIdx.y*blockDim.y;

	if( y<n_rows_c && x<n_cols_c ){
		int i;
		for(i=0 ; i<n_cols_a ; i++){
			value += dd_mat_a[y][i] * dd_mat_b[i][x];
		}
		dd_mat_c[y][x]=value;
	} 
}


void fill_kernel_3x3(int** mat, int n, int m){
    mat[0][0]=2; mat[0][1]=3; mat[0][2]=2;
    mat[1][0]=2; mat[1][1]=0; mat[1][2]=2;
    mat[2][0]=2; mat[2][1]=3; mat[2][2]=2;
}



void fill(int** mat, int n, int m){
    srand(time(0));
	int i,j; 
	for(i=0; i<n ;i++){
		for(j=0; j<m ;j++)
			mat[i][j] = rand()%5;
			//mat[i][j] = 1;
	}
}


void fill_value(int** mat,int n, int m, int value=0){
	int i,j; 
	for(i=0;i<n;i++)
		for(j=0;j<m;j++)
			mat[i][j] = value;
}


void print(int** mat,int n, int m){
	int i,j; 
	for(i=0; i<n ;i++){
		for(j=0; j<m ;j++)
			printf("%d ",mat[i][j]);
		printf("\n");
	}
}



void create_copy(int**& mat, int**& d_mat, int**& dd_mat, int n_rows, int n_cols, int fillValue=-1){
	
	int i;

	int size_row = sizeof(int*) * n_rows;

	d_mat = (int**) malloc(size_row);
	cudaMalloc((void**)& d_mat[0], sizeof(int) * n_rows * n_cols );
	cudaMemcpy(  d_mat[0], mat[0], sizeof(int) * n_rows * n_cols ,cudaMemcpyHostToDevice);

	for( i=1 ; i<n_rows ; i++ ){
		d_mat[i] = (d_mat[0]+i*n_cols);
	}	
	
	cudaMalloc((void***)& dd_mat, size_row );
	cudaMemcpy( dd_mat, d_mat, size_row, cudaMemcpyHostToDevice );

}



void create(int**& mat, int**& d_mat, int**& dd_mat, int n_rows, int n_cols, int fillValue=-1){
	
	int i;
	mat 	= (int** )malloc(sizeof(int*) * n_rows 			);	
	mat[0] 	= (int*  )malloc(sizeof(int ) * n_rows * n_cols );	

	for( i=1 ; i<n_rows ; i++ ){
		mat[i] = mat[i-1]+n_cols;
	}
	if(fillValue==-1){
		fill(mat,n_rows,n_cols);	
	}
	else{
		fill_value(mat,n_rows,n_cols,fillValue);
	}

	int size_row = sizeof(int*) * n_rows;

	d_mat = (int**) malloc(size_row);
	cudaMalloc((void**)& d_mat[0], sizeof(int) * n_rows * n_cols );
	cudaMemcpy(  d_mat[0], mat[0], sizeof(int) * n_rows * n_cols ,cudaMemcpyHostToDevice);

	for( i=1 ; i<n_rows ; i++ ){
		d_mat[i] = (d_mat[0]+i*n_cols);
	}	
	
	cudaMalloc((void***)& dd_mat, size_row );
	cudaMemcpy( dd_mat, d_mat, size_row, cudaMemcpyHostToDevice );

}



void create_kernell_static(int**& mat, int**& d_mat, int**& dd_mat, int n_rows, int n_cols){
	
	int i;
	mat 	= (int** )malloc(sizeof(int*) * n_rows 			);	
	mat[0] 	= (int*  )malloc(sizeof(int ) * n_rows * n_cols );	

	for( i=1 ; i<n_rows ; i++ ){
		mat[i] = mat[i-1]+n_cols;
	}

	fill_kernel_3x3(mat,n_rows,n_cols); 

	int size_row = sizeof(int*) * n_rows;
	d_mat = (int**) malloc(size_row);
	cudaMalloc((void**)& d_mat[0], sizeof(int) * n_rows * n_cols );
	cudaMemcpy(  d_mat[0], mat[0], sizeof(int) * n_rows * n_cols ,cudaMemcpyHostToDevice);

	for( i=1 ; i<n_rows ; i++ ){
		d_mat[i] = (d_mat[0]+i*n_cols);
	}	
	
	cudaMalloc((void***)& dd_mat, size_row );
	cudaMemcpy( dd_mat, d_mat, size_row, cudaMemcpyHostToDevice );

}


void create_kernell_random(int**& mat, int**& d_mat, int**& dd_mat, int n_rows, int n_cols){
	
	int i,j; 

	mat 	= (int** )malloc(sizeof(int*) * n_rows 			);	
	mat[0] 	= (int*  )malloc(sizeof(int ) * n_rows * n_cols );	

	for( i=1 ; i<n_rows ; i++ ){
		mat[i] = mat[i-1]+n_cols;
	}

	srand(time(0));
	for(i=0; i<n_rows ;i++){
		for(j=0; j<n_cols ;j++)
			mat[i][j] = rand()%3+1;
			//mat[i][j] = 1;
	}

	int size_row = sizeof(int*) * n_rows;

	d_mat = (int**) malloc(size_row);
	cudaMalloc((void**)& d_mat[0], sizeof(int) * n_rows * n_cols );
	cudaMemcpy(  d_mat[0], mat[0], sizeof(int) * n_rows * n_cols ,cudaMemcpyHostToDevice);

	for( i=1 ; i<n_rows ; i++ ){
		d_mat[i] = (d_mat[0]+i*n_cols);
	}	
	
	cudaMalloc((void***)& dd_mat, size_row );
	cudaMemcpy( dd_mat, d_mat, size_row, cudaMemcpyHostToDevice );

}





int main(int argc, char *argv[]){

	

	/////////////////////////////////////////////////////
	

	int n_rows_a = 10;
	int n_cols_a = 10;

	int n_rows_b = 3;  //n_kernel
	int n_cols_b = 3;  //n_kernel

	int n_rows_c = 10;
	int n_cols_c = 10;



	int** mat_a; int** d_mat_a;	 int** dd_mat_a;	
	int** mat_b; int** d_mat_b;	 int** dd_mat_b;	
	int** mat_c; int** d_mat_c;	 int** dd_mat_c;	

	//create_copy( mat_a, d_mat_a, dd_mat_a, n_rows_a, n_cols_a	);
	create( mat_a, d_mat_a, dd_mat_a, n_rows_a, n_cols_a	);
	
	create_kernell_static( mat_b, d_mat_b, dd_mat_b, n_rows_b, n_cols_b 	); 
	//create_kernell_random( mat_b, d_mat_b, dd_mat_b, n_rows_b, n_cols_b 	);

	create( mat_c, d_mat_c, dd_mat_c, n_rows_c, n_cols_c, 0	);



	/////////////////////////////////////////

	dim3 blockNum(WIDTH_TILE,WIDTH_TILE,1);
	dim3 grid((int)(n_cols_c-1+blockNum.x)/blockNum.x,(int)(n_rows_c-1+blockNum.y)/blockNum.y,1);
	printf("ty: %d, tx: %d\n",(int)(n_rows_c-1+blockNum.y)/blockNum.y, (int)(n_cols_c-1+blockNum.x)/blockNum.x);
	printf("grid_row: %d, grid_col: %d\n",grid.x , grid.y );

	////////////////////////////////////////////////////
	
	//convolution<<<grid,blockNum>>>(dd_mat_a, n_rows_a, n_cols_a, dd_mat_b, n_rows_b, n_cols_b, dd_mat_c, n_rows_c, n_cols_c);
	convolution_shared<<<grid,blockNum>>>(dd_mat_a, n_rows_a, n_cols_a, dd_mat_b, n_rows_b, n_cols_b, dd_mat_c, n_rows_c, n_cols_c);

	//convolution_complete<<<grid,blockNum>>>(dd_mat_a, n_rows_a, n_cols_a, dd_mat_b, n_rows_b, n_cols_b, dd_mat_c, n_rows_c, n_cols_c);



	//matrix_mult_shared<<<grid,blockNum>>>(dd_mat_a, n_rows_a, n_cols_a, dd_mat_b, n_rows_b, n_cols_b, dd_mat_c, n_rows_c, n_cols_c);

	//matrix_mult<<<grid,blockNum>>>(dd_mat_a, n_rows_a, n_cols_a, dd_mat_b, n_rows_b, n_cols_b, dd_mat_c, n_rows_c, n_cols_c);
	

    /////////////////////////////////////////////////////

	cudaMemcpy(mat_c[0],d_mat_c[0],sizeof(int)*n_rows_c*n_cols_c,cudaMemcpyDeviceToHost);		
	
	
	printf("//////////////////\n");
	printf("//////////////////\n");
	print(mat_a,n_rows_a,n_cols_a);
	printf("//////////////////\n");
	print(mat_b,n_rows_b,n_cols_b);

	printf("//////////////////\n");
	print(mat_c,n_rows_c,n_cols_c);
	
	
  
	cudaFree(dd_mat_a);
	cudaFree(dd_mat_b);
	cudaFree(dd_mat_c);
	cudaFree(d_mat_a);
	cudaFree(d_mat_b);
	cudaFree(d_mat_c);
  	
  	free(mat_a);
  	free(mat_b);
  	free(mat_c);
  

	return 0;
}

