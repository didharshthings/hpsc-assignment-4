//siddharth
// serial matrix multiplication

#include <stdio.h>
#include <math.h>


void print_matrix(int row, int column, double matrix[row][column]){
  for(int i=0; i<row; i++)
    {
      printf("[%d] - ", i);
      for(int j=0; j<column; j++)
        {
    printf("%.1lf ", matrix[i][j]);
        }
      printf("\n");
      fflush(stdout);
    }
}

int main(int argc, char *argv[])
{
  int i,j,k;
  int size = 8;
  double matrix_a[size][size];
  double matrix_b[size][size];
  double matrix_c[size][size];


  for(i=0;i<size;i++) 
  {
    for(j=0;j<size;j++)
      {
        matrix_a[i][j] = rand()%10;
        matrix_b[i][j] = rand()%10;
      }
  }

 //naive
    // for (i = 0; i < SIZE; i++)    
    // {
    //   for (j = 0; j < SIZE; j++) 
    //   {
      // C[i][j] = 0.0;
    //     for (k = 0; k < SIZE; k++)
     //  {
     //   C[i][j] = C[i][j] + A[i][k]*B[k][j];
    //      }
    //    }
    // }


//quicker
 for (i = 0; i < size; i++)    
    {
      for (k = 0; k < size; k++) 
      {
      matrix_c[i][j] = 0.0;
        for (j = 0; j < size; j++)
      {
        matrix_c[i][j] += matrix_a[i][k]*matrix_b[k][j];
        }
       }
    }   

    print_matrix(size,size, matrix_a);
    printf("\n");
    print_matrix(size,size, matrix_b);
    printf("Product \n");
    print_matrix(size,size, matrix_c);

  return 0;
}