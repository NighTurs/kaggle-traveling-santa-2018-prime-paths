/******************************************************************************\
*				  				 Files Manipulation							 *
\******************************************************************************/

#include "defs.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<string>
#include<cstring>
#include<fstream>

#define CHAR_LEN 100



/****************************************************************/
/* Read the problem instance (adapted from DTSP generator - Yang et al.) */
/****************************************************************/
void read_problem(const char* filename){
	int i, j, id;
	char line[CHAR_LEN], * keywords,Delimiters[] = " :=\n\t\r\f\v", * tempChar;
	
  std::ifstream fin(filename);
  while((fin.getline(line, CHAR_LEN-1))){
		if(!(keywords = strtok(line, Delimiters)))
  			continue;
		if(!strcmp(keywords, "DIMENSION")){			
  			if(!sscanf(strtok(NULL, Delimiters), "%d", &n_cities)){
				printf("DIMENSION error\n");
				exit(0);
  			}
			coord_x=aloc_vectord (n_cities);
			coord_y=aloc_vectord (n_cities);
		}
		else if(!strcmp(keywords, "EDGE_WEIGHT_TYPE")){
  			if(!(tempChar=strtok(NULL, Delimiters))){
				printf("EDGE_WEIGHT_TYPE error\n");
				exit(0);
  			}
			 if(strcmp(tempChar, "EUC_2D")){
				printf("not EUC_2D\n");
				exit(0);
			 }
		}
		else if(!strcmp(keywords, "NODE_COORD_SECTION")){
  			if(n_cities!=0){
  				for(i=0; i<n_cities; i++){
					//coordinates
					fin>>id;       
					fin>>coord_x[i]>>coord_y[i];

  				}
    		}
		}
  }
  fin.close();

  // building the weight matrix
	 if (n_cities<=15000){
		  W = aloc_matrixi(n_cities,n_cities);
		  for(i=0; i<n_cities; i++){
				for(j=i; j<n_cities; j++){
						if (i==j)
							W[i][j]=0;
						else {
							W[i][j]=round( sqrt( pow(coord_x[i]-coord_x[j],2) +pow(coord_y[i]-coord_y[j],2) ) );		// compute the distance between the cities i and j
							W[j][i]=W[i][j];            
						}
				}
		   }
	}

	
}


