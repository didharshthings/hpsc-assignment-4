//Siddharth
//matrix multiplication using Fox's Algorithm
//based on code from Pacheco

#include "mpi.h"
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

typedef struct
{
  int pid;
  int num_procs; 
  MPI_Comm grid_comm;
  MPI_Comm row_comm;
  MPI_Comm col_comm;
  int grid_size;
  int curr_row;
  int curr_col;
  int grid_pid;
}Grid;


void setup_grid(Grid* grid)
{
  int pid;
  int dimensions[2];
  int wrap_around[2];
  int coords[2];
  int local_coords[2];

  MPI_Comm_rank(MPI_COMM_WORLD,&pid);
  MPI_Comm_size(MPI_COMM_WORLD,&(grid->num_procs));

  printf("setting up the grid \n");
  //assume num_procs is a perfect square

  grid->grid_size = (int)sqrt((double)grid->num_procs);
  dimensions[0]=dimensions[1]=grid->grid_size;
  
  wrap_around[0]=wrap_around[1]=1;

  MPI_Cart_create(MPI_COMM_WORLD,2,dimensions,wrap_around,1,&(grid->grid_comm));
  MPI_Comm_rank(grid->grid_comm, &(grid->grid_pid));
  MPI_Cart_coords(grid->grid_comm,grid->grid_pid,2,coords);
  
  grid->curr_row = coords[0];
  grid->curr_col = coords[1];

  //row communicators
  local_coords[0] = 0;
  local_coords[1] = 1;
  MPI_Cart_sub(grid->grid_comm,local_coords, &(grid->row_comm));

  //column communicators
  local_coords[0] = 0;
  local_coords[1] = 0;
  MPI_Cart_sub(grid->grid_comm,local_coords,&(grid->col_comm));

  printf("grid setup complete for process - %d\n",pid);
  printf ("Number of Processes is %d\n", grid->num_procs);
  printf ("Grid Order is %d\n", grid->grid_size);
  printf ("Current Process Coordinates are (%d, %d)\n", grid->curr_row, grid->curr_col);
  printf ("Process rank in Grid is %d\n", grid->grid_pid); 
}

double**  matrix_product(double** a, double** b, int size)
{
  int i,j,k;
  printf("multiplying matrices\n");
  double **temp_c;
  temp_c = (double**) malloc(size*sizeof(double*));
  for(i=0;i<size;i++)
    {
      *(temp_c + i) = (double*) malloc(size*sizeof(double));
    }
  for(i = 0; i<size; i++)
    {
      for(j=0;j<size;j++)
      {
      temp_c[i][j] = 0.0;
    }
    }

  for(i=0;i<size;i++)
    {
      for(j=0;j<size;j++)
	{
	  for(k=0;k<size;k++)
	    {
	      temp_c[i][j] += a[i][k]*b[k][j];
	      printf("c[%d][%d] - %d \n",i,j,temp_c[i][j]);
	      printf("a[%d][%d] - %d \n",i,j,a[i][k]);
	      printf("b[%d][%d] - %d \n",i,j,b[k][j]);
	    }
	}
    }
  printf("done multiplying\n");
  return temp_c;
}

void fox_mult(Grid* grid,double** local_a,double** local_b, double** local_c,int size)
{

  int i, source, x,y ,dest, root,j,k;
  int local_size = size/grid->grid_size;
  printf("local size %d\n",local_size);
  double temp_a[local_size][local_size];
  MPI_Status status;

 
  source = (grid->curr_row + 1) % grid->grid_size;
  dest = (grid->curr_row + grid->grid_size - 1) % grid->grid_size;

  for(i=0;i < grid->grid_size;i++)
    {
      root = (grid->curr_row + i) % grid->grid_size;
      if(root == grid->curr_col)
	{
	  MPI_Bcast(local_a,local_size*local_size,MPI_DOUBLE,root,grid->row_comm);
	  local_c = matrix_product(local_a,local_b,local_size);
	}
      else
	{
	  MPI_Bcast(temp_a,local_size*local_size,MPI_DOUBLE,root,grid->row_comm);
	  local_c = matrix_product(temp_a,local_b,local_size);
	}
      printf("send_Recv source -%d, dest - %d}\n",source, dest);
      MPI_Sendrecv_replace(local_b,1,MPI_DOUBLE,dest,0,source,0,grid->col_comm,&status);
    }

}

void print_matrix(int row, int column, double** matrix)
{
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

double** read_matrix(double** matrix, Grid* grid_info, int size)
{
  int x, y, local_size, local_pid, coords[2];
  int i,j;
  double** local_matrix;
  local_size = size/grid_info->grid_size;
  local_matrix = (double**) malloc(local_size*sizeof(double*));
  for(i=0;i<local_size;i++)
    {
      *(local_matrix + i) = (double*) malloc(local_size*sizeof(double));
    }
  for(i = 0; i<local_size; i++)
    {
      for(j=0;j<local_size;j++)
	{
	  local_matrix[i][j] = 0.0;
	}
    }

  MPI_Comm_rank(grid_info->grid_comm, &local_pid);
  MPI_Cart_coords(grid_info->grid_comm, local_pid, 2, coords);
  for(x =0 ;x<local_size;x++)
    {
      for(y=0;y<local_size;y++)
	{
	  local_matrix[x][y] = matrix[local_size*coords[0]+x][local_size*coords[1]+y];
	}
    }
  printf("read matrix  - \n");
  for(i=0; i<local_size; i++)
    {
      printf("[%d] - ", i);
      for( j=0; j<local_size; j++)
	{
	  printf("%.1lf ", local_matrix[i][j]);
	}
      printf("\n");
      fflush(stdout);
    }
  return local_matrix;
}

void write_matrix(double** local_matrix, double** matrix, Grid* grid_info,int size)
{
  int x,y, local_size, local_pid, coords[2], i,j , num_procs;
  MPI_Status status;

  local_size = size/grid_info->grid_size;
  MPI_Comm_rank(grid_info->grid_comm,&local_pid);
  MPI_Cart_coords(grid_info->grid_comm,local_pid,2,coords);

  if(local_pid)
    {
      MPI_Send(&local_matrix[0][0], local_size*local_size,MPI_DOUBLE,0,0,grid_info->grid_comm);
      MPI_Send(coords,2,MPI_INT,0,0,grid_info->grid_comm);
    }
  else
    {
      for(x=0;x<local_size;x++)
	{
	  for(y=0;y<local_size;y++)
	    {
	      matrix[local_size*coords[0]+x][local_size*coords[1]+y]=local_matrix[x][y];
	    }
	}
      MPI_Comm_size(grid_info->grid_comm, &num_procs);
      for ( i=1;i<num_procs;i++)
	{
	  MPI_Recv(&local_matrix[0][0],local_size*local_size,MPI_DOUBLE,i,0,grid_info->grid_comm,&status);
	  MPI_Recv(coords,2,MPI_INT,i,0,grid_info->grid_comm,&status);
	  for(x=0;x<local_size;x++)
	    {
	      for(y=0;y<local_size;y++)
		{
		  matrix[local_size*coords[0]+x][local_size*coords[1]+y]=local_matrix[x][y];
		}
	    }
	}
    }
}

int main(int argc, char** argv)
{
  int i,j;
  int size,local_size;
  int pid;
  int num_procs; //assume its a square
  double **local_a;
  double **local_b;
  double **local_c;
  double **matrix_A;
  double **matrix_B;
  double **matrix_C;
  Grid grid_info;

  size = 16; 

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&pid);
  setup_grid(&grid_info);

  //initialize  matrix_A, matrix_B and matrix_C
  matrix_A = (double**) malloc(size*sizeof(double*));
  matrix_B = (double**) malloc(size*sizeof(double*));
  matrix_C = (double**) malloc(size*sizeof(double*));
  
  for(i=0;i<size;i++)
    {
      *(matrix_A + i) = (double*) malloc(size*sizeof(double));
      *(matrix_B + i) = (double*) malloc(size*sizeof(double));
      *(matrix_C + i) = (double*) malloc(size*sizeof(double));
    }
  
  for(i = 0; i<size; i++)
    {
      for(j=0;j<size;j++)
	{
	  matrix_A[i][j] = 1;
	  matrix_B[i][j] = 1;
	  matrix_C[i][j] = 0.0;
	}
    }
 
  
  local_size = size/grid_info.grid_size;

  MPI_Barrier(MPI_COMM_WORLD);
  // initializaing local matrices
  local_a = (double**) malloc(local_size*sizeof(double*));
  local_b = (double**) malloc(local_size*sizeof(double*));
  local_c = (double**) malloc(local_size*sizeof(double*));

  for(i=0;i<local_size;i++)
    {
      *(local_a + i) = (double*) malloc(local_size*sizeof(double));
      *(local_b + i) = (double*) malloc(local_size*sizeof(double));
      *(local_c + i) = (double*) malloc(local_size*sizeof(double));
    }

  //intialize local matrix
  for(i = 0; i<local_size; i++)
    {
      for(j=0;j<local_size;j++)
	{
	  local_a[i][j] = 0.0;
	  local_c[i][j] = 0.0;
	  local_b[i][j] = 0.0;
	}
    }

  local_a = read_matrix(matrix_A, &grid_info, size); //distributes matrix a
  local_b = read_matrix(matrix_B, &grid_info, size); //distributes matrix b
  
  fox_mult(&grid_info,local_a,local_b,local_c,size);

  write_matrix(local_c,matrix_C, &grid_info, size); //collates matrix c
  if (pid == 0)
    {
    print_matrix(size,size, matrix_A);
    print_matrix(size,size, matrix_B);
    printf("product -\n");
    print_matrix(size,size, matrix_C);
    }

  MPI_Finalize();
  return 0;
}


  
