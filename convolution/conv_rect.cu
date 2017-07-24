#include <stdio.h>
#include <stdlib.h>

#define MAT_ROW 4
#define MAT_COL 6
#define CHANNELS 1
#define MASK_WIDTH  3
#define MASK_RADIUS MASK_WIDTH/2
#define O_TILE_WIDTH 12
#define BLOCK_WIDTH (O_TILE_WIDTH + (MASK_WIDTH-1))

void print_matrix(float* a,int n,int m)
{
        int i,j;
        for(i=0;i<n;i++)
        {
                for(j=0;j<m;j++)
                {
                        printf("%f ",a[i*m+j]);
                }
                printf("\n");
        }
}

void fill_mat(float* a,int n,int m)
{
        //srand(time(NULL));
        int i,j;
        for(i=0;i<n;i++)
        {
                for(j=0;j<m;j++)
                {
                        //a[i*n+j] = (rand()%2+1)*1.0;
                        a[i*m+j] = 1.0;
                }
        }
}

__global__ void convolution_shared(float *in, float* out,const float* __restrict__ M,int height, int width, int channels)
{

  float sum, pixel, maskVal;

  __shared__ float Ns[BLOCK_WIDTH][BLOCK_WIDTH];

  int tx = threadIdx.x;
  int ty = threadIdx.y;

  int row_o = blockIdx.y*O_TILE_WIDTH + ty;
  int col_o = blockIdx.x*O_TILE_WIDTH + tx;

  int row_i = row_o - MASK_RADIUS;
  int col_i = col_o - MASK_RADIUS;

  
  for (int c = 0; c < channels; c++) {

    if ( (row_i >= 0) && (row_i < height) &&
        (col_i >= 0) && (col_i < width) ) {
      Ns[ty][tx] = in[(row_i*width + col_i)*channels + c]; 
    }
    else {
      Ns[ty][tx] = 0.0f;
    }
    __syncthreads();

    sum = 0.0;
    if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {
      for (int y = 0; y < MASK_WIDTH; y++){
        for (int x = 0; x < MASK_WIDTH; x++){
          pixel = Ns[ty + y][tx + x];
          maskVal = M[y*MASK_WIDTH + x];
          sum += pixel*maskVal;
        }
      }
      if (row_o < height && col_o < width) {
        //out[ (row_o * width + col_o) * channels + c] = min(max(0.0f,sum),1.0f);
        out[ (row_o * width + col_o) * channels + c] = sum;
      }
    }
   // __syncthreads();
  }
}


int main()
{
        float *mat,*d_mat;
        float *mask,*d_mask;
        float *result,*d_result;

        //float elapsed_time=0;
        //cudaEvent_t start,stop;
        //cudaEventCreate(&start);
        //cudaEventCreate(&stop);

        int mat_size = MAT_ROW*MAT_COL*sizeof(float);
        int mask_size = MASK_WIDTH*MASK_WIDTH*sizeof(float);

        mat = (float*) malloc(mat_size);
        result = (float*) malloc(mat_size);
        mask = (float*) malloc(mask_size);

        fill_mat(mat,MAT_ROW,MAT_COL);
        fill_mat(mask,MASK_WIDTH,MASK_WIDTH);
        printf("Printing Matrix \n");
        print_matrix(mat,MAT_ROW,MAT_COL);
        printf("Printing Mask\n");
        print_matrix(mask,MASK_WIDTH,MASK_WIDTH);
        printf("\n");

        cudaMalloc((void** )&d_mat,mat_size);
        cudaMalloc((void** )&d_result,mat_size);
        cudaMalloc((void** )&d_mask,mask_size);

        cudaMemcpy(d_mat,mat,mat_size,cudaMemcpyHostToDevice);
        cudaMemcpy(d_mask,mask,mask_size,cudaMemcpyHostToDevice);

        dim3 my_block(BLOCK_WIDTH,BLOCK_WIDTH);
        dim3 my_grid((MAT_COL + BLOCK_WIDTH-1)/my_block.x,(MAT_ROW + BLOCK_WIDTH-1)/my_block.y);

 		convolution_shared<<<my_grid,my_block>>>(d_mat, d_result, d_mask, MAT_ROW,MAT_COL,CHANNELS);
        cudaMemcpy(result,d_result,mat_size,cudaMemcpyDeviceToHost);
        printf("Printing result\n");
        print_matrix(result,MAT_ROW,MAT_COL);

        //printf("Elapsed time: %f\n",elapsed_time);
        return 0;
}

