//siddharth
// allgather matrix mutltiplication

#include <stdio.h>
#include "mpi.h"
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
  int pid;
  int num_procs;
  int i,j,k;
  int tag = 0;
  int from,to;
  int size = 8;
  double matrix_a[size][size];
  double matrix_b[size][size];
  double matrix_c[size][size];
  MPI_Status status;

  
  MPI_Init (&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 

  if (pid==0) {
  for(i=0;i<size;i++) 
  {
    for(j=0;j<size;j++)
      {
        matrix_a[i][j] = rand()%10;
        matrix_b[i][j] = rand()%10;
      }
  }
  }

  MPI_Bcast (matrix_b, size*size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  from = pid * size/num_procs;
  to = (pid+1) * size/num_procs;
  MPI_Scatter (matrix_a, size*size/num_procs, MPI_DOUBLE, matrix_a[from], size*size/num_procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  printf("for process - %d, multiply row from %d to %d \n",pid,from, to);
  for (i=from; i<to; i++) 
    for (j=0; j<size; j++) {
      matrix_c[i][j]=0;
      for (k=0; k<size; k++)
    matrix_c[i][j] += matrix_a[i][k]*matrix_b[k][j];
    }

  MPI_Allgather (matrix_c[from], size*size/num_procs, MPI_DOUBLE, matrix_c, size*size/num_procs, MPI_DOUBLE, MPI_COMM_WORLD);
 
  if (pid ==0) {
    print_matrix(size,size, matrix_a);
    print_matrix(size,size, matrix_b);
    printf("Product \n");
    print_matrix(size,size, matrix_c);
  }

  MPI_Finalize();
  return 0;
}