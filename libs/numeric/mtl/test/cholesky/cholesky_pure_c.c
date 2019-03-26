
#include <math.h>
#include <time.h>
#include <stdio.h>

const int num_rows=8000;
const int num_cols=num_rows;
double matrix[num_rows][num_cols];
//double matrix[1152][1152];

void display()
{   
    int i,j; 
	for ( i=0 ; i<num_rows; i++ ){
		for( j=0; j<num_cols;  j++ ){		   
		  printf("%0.5f  ", matrix[i][j]);
		}
		printf("\n");
	}
	return;	
}
 
void fill_matrix(){
     
    double x= 1.0;     
    int i,j; 
    for( i=0;i<num_rows;i++) {
       for( j=0;j<=i;j++){
         if(i!=j){
	          matrix[i][j]=x; matrix[j][i]=x; 
	          x=x+1.0; 
	       }
       }
    }
  
    double rowsum;
    for( i=0;i< num_rows ;i++){
       rowsum=0.0;
       for( j=0;j< num_cols ;j++){
         if(i!=j){
	          rowsum += matrix[i][j]; 
	       }
       }
       matrix[i][i]=rowsum*2;
    }       
}


void verify_matrix()
{ 
	int i,j,k;	
	for(i=0; i<num_rows;  i++ ){
		for(j=0 ; j<num_rows;  j++ ){		
			double  element=0;
			for (k=0; k<num_cols; k++ )		 
				element += matrix[i][k] * matrix[j][k];		
		printf("%2.5f  ",element);		 
		}
	 	printf("\n");
	}
	return;
}


 




int main()
{	
  int i,j,k;  
  time_t starttime,endtime;
  struct tm *timeinfo;
  

  time (&starttime);
  timeinfo = localtime (&starttime);
  
  printf("---PURE C-------order = %d   -------------------->Load start: %s",
	 num_rows,  asctime (timeinfo));
  
  fill_matrix();
  //  display();
  // printf("=====================================================\n");

  time (&starttime);
  timeinfo = localtime (&starttime);
  printf("---PURE C-------order = %d   -------------------->Start date and time are: %s",  
	 num_rows, asctime (timeinfo));

  for( k=0; k<num_rows;  k++ ){
    matrix[k][k] = sqrt( matrix[k][k]);     
    
    //     for( i=k+1; i<num_rows; i++ ) matrix[i][k] /= matrix[k][k];
    
    for( i=k+1; i<num_rows;  i++){
      matrix[i][k] /= matrix[k][k];
      double d = matrix[i][k];
      for ( j=k+1; j<=i; j++ ) 
	matrix[i][j] -= d * matrix[j][k];		    
    }
  }
  
  time (&endtime);
  timeinfo = localtime (&endtime);
    printf("---PURE C-------order = %d     ------------------->End date and time are: %s", 	 num_rows,  asctime (timeinfo));

 
  printf("\nTOTAL TIME TAKEN for PURE C for order = %d  : %d  secs\n\n",num_rows,  endtime - starttime);
  

  //  printf("=====================================================\n");
  
  /* for( i=0 ; i<num_rows;  i++ ){
    for ( j=i+1; j<num_cols; j++ )			 
    matrix[i][j]=0; }*/
    //  display();
    // printf("=====================================================\n");
  
  // verify_matrix();
  return 0;
  
	
}
