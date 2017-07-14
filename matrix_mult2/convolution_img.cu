// Luigy Machaca Arcana
// Computer science - Arequipa, Perú  2017


#include <stdlib.h>
#include <stdio.h>

#include <fstream>
#include <iostream>
#include <string>


using namespace std;




#define WIDTH_TILE 32



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
		offset = offset>0?offset:0;
		offset = (int)offset%254 + 10;
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
		offset = -1/256*offset;
		offset = offset>0?offset:0;
		offset = (int)offset%255 + 1;
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
    mat[0][0]=-1; mat[0][1]=0; mat[0][2]=1;
    mat[1][0]=-1; mat[1][1]=0; mat[1][2]=1;
    mat[2][0]=-1; mat[2][1]=0; mat[2][2]=1;
}

void fill_kernel_5x5(int** mat, int n, int m){
//-1/256	// 1 4 6 4 1
		// 4 16 24 16 4
		// 6 24 -476 24 6
		// 4 16 24 16 4	
		// 1 4 6 4 1
	mat[0][0]=1; mat[0][1]=4 ; mat[0][2]=6   ; mat[0][3]=4 ; mat[0][4]=1;
	mat[1][0]=4; mat[1][1]=16; mat[1][2]=24  ; mat[1][3]=16; mat[1][4]=4;
	mat[2][0]=6; mat[2][1]=24; mat[2][2]=-476; mat[2][3]=24; mat[2][4]=6;
	mat[3][0]=4; mat[3][1]=16; mat[3][2]=24  ; mat[3][3]=16; mat[3][4]=4;
	mat[4][0]=1; mat[4][1]=4 ; mat[4][2]=6   ; mat[4][3]=4 ; mat[4][4]=1;

}



void fill(int** mat, int n, int m){
    srand(time(0));
	int i,j; 
	for(i=0; i<n ;i++){
		for(j=0; j<m ;j++)
			//mat[i][j] = rand()%3+1;
			mat[i][j] = 1;
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
	//fill_kernel_5x5(mat,n_rows,n_cols); 

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

	

	string title1,title2;
	char rows[5];
	char cols[5];
	int n_rows = -1;
	int n_cols = -1;


	ifstream myReadFile;
	myReadFile.open("img/mario.pgm");

	char out_temp[10];
	
	int** mat_a;

	if (myReadFile.is_open()){

		std::getline(myReadFile,title1);
		std::getline(myReadFile,title2);

		myReadFile >> rows;
		n_rows = atoi(rows);
		//n_rows = 30;
		cout << n_rows << endl;

		myReadFile >> cols;
		n_cols = atoi(cols);
		//n_cols = 30;
		cout << n_cols << endl;

		/////////////////////////////////////////////////////////////
		mat_a 		= (int** )malloc(sizeof(int*) * n_rows 			);	
		mat_a[0] 	= (int*  )malloc(sizeof(int ) * n_rows * n_cols );	
		
		for( int i=1 ; i<n_rows ; i++ ){
			mat_a[i] = mat_a[i-1]+n_cols;
		}

		/////////////////////////////////////////////////////////////
		int n_temp;
		for(int i=0 ; i<n_rows ; i++){
			for(int j=0 ; j<n_cols ; j++){
				if(!myReadFile.eof()){
					myReadFile >> out_temp;
					n_temp = atoi(out_temp);
					mat_a[i][j] = n_temp;
					//cout << n_temp-1000 << endl;	
				}
			}
		}
		//while (!myReadFile.eof()){
		//	myReadFile >> out_temp;
		//	n_temp = atoi(out_temp);
		//	cout << n_temp-1000 << endl;
		//}
	}
	myReadFile.close();





	/////////////////////////////////////////////////////

	int n_rows_a = n_rows;
	int n_cols_a = n_cols;

	int n_rows_b = 5;  //n_kernel
	int n_cols_b = 5;  //n_kernel

	int n_rows_c = n_rows;
	int n_cols_c = n_cols;



	//int** mat_a; int** d_mat_a;	 int** dd_mat_a;	
	//int** mat_a;

				 int** d_mat_a;	 int** dd_mat_a;	
	int** mat_b; int** d_mat_b;	 int** dd_mat_b;	
	int** mat_c; int** d_mat_c;	 int** dd_mat_c;	

	create_copy( mat_a, d_mat_a, dd_mat_a, n_rows_a, n_cols_a	);
	//create( mat_a, d_mat_a, dd_mat_a, n_rows_a, n_cols_a	);
	
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

	convolution_complete<<<grid,blockNum>>>(dd_mat_a, n_rows_a, n_cols_a, dd_mat_b, n_rows_b, n_cols_b, dd_mat_c, n_rows_c, n_cols_c);



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
	
	
	//////////////////////////////////////////////
	ofstream myfile;
	myfile.open ("mario3.pgm");
	myfile << title1 <<endl;
	myfile << title2 <<endl;
	myfile << n_rows_c <<" "<< n_cols_c <<endl;

  	for(int i=0 ; i<n_rows_c ; i++){
		for(int j=0 ; j<n_cols_c ; j++){
			myfile << mat_c[i][j] <<endl;
		}
	}

	myfile.close();
	//////////////////////////////////////////////

  
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

