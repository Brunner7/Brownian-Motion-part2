
/*
 ============================================================================
 Name        : task1.c
 Author      : ich
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#define root 4		//root of number of particles, I want a quadratic initial lattice
#define Pe 20			//Peclet number
#define t_ini 0			//t0
#define t_fin 5		//tfin
#define timestep 0.0001	//time step size
#define L 85			//box length, must be or enter manually
#define cutoff 5		//cutoff radius
#define epsilon_wca 10	//energy depth of wca potential
#define size_of_cell 5		//must be chosen int and so that L/size_of_cell=int
#define zeros_in_cell 100	//number of zeros in 3d matrix
#define savesteps ( t_fin - t_ini ) / timestep / 200 // number of iterations/100, to get points where timestep will be saved
#define N root * root		//number of particles
#define cells L / size_of_cell
#define sigma_wca cutoff * 0.890899	//sigma in wca
#define PI2 2 * M_PI
#define sqrttimestep2 sqrt(2 * timestep)
//#define sigma6 sigma_wca * sigma_wca * sigma_wca * sigma_wca * sigma_wca * sigma_wca
//#define sigma12 sigma6 * sigma6

//defining a function that generates normal-distributed random values
double
normal_distribution (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
}

int main(void) {
	int i, j, k, m, n, p, p1, p2, count=0;
	double t,  time_spent = 0.0, wca_x, wca_y, distance_inv, dis_inv7;
	double sigma_w=sigma_wca, sigma_6=sigma_w*sigma_w*sigma_w*sigma_w*sigma_w*sigma_w, sigma_12=sigma_6*sigma_6;
	double*x,*y,*phi;
	x=malloc(N*sizeof(double));
	y=malloc(N*sizeof(double));
	phi=malloc(N*sizeof(double));
	//start runtime counter
	clock_t begin = clock();
	//open file to write data to
	FILE *values;
		values = fopen("values.txt", "w");
		fprintf(values,"current time \t particle1: (x y phi) particle2: (x y phi) etc \n 0 ");
	//setup the even spaced lattice in x,y
	for( i = 0; i < N; i = i + root ){
		for( j = 0; j < root; j++ ){
			x[ i + j ] = ( j + .5 ) * L / ( 2.0 * root ) + L / 4;
			y[ i + j ] = ( i + root / 2.0 ) / root * L / ( 2.0 * root ) + L / 4;
		}
	}
	//setup the initial angles, randomly distributed
	for ( i = 0; i < N; i++ ){
		phi[i] = PI2 * normal_distribution(0,1);
		fprintf(values,"%.5f  %.5f  %.5f ",x[i],y[i],phi[i]);
	}
	//for my sanity
//	int cell_list[cells+2][cells+2][zeros_in_cell]={0}, track_list[cells+2][cells+2]={0};
//	for ( k = 0; k < N; k++){
//			i = (int) x[k] / size_of_cell + 1;
//			j = (int) y[k] / size_of_cell + 1;
//			cell_list[i][j][track_list[i][j]] = k;
//			if ( track_list[i][j] <= zeros_in_cell ){
//				track_list[i][j] += 1;
//			}
//			if (track_list[i][j]==zeros_in_cell){
//				//printf( "Error, too many particles in list!" );
//
//			}
//	}
//
//	for ( k = 0; k < cells + 2; k++){
//		i=0; while ((int)cell_list[cells][k+1][i]!=0){
//
//		cell_list[0][0][i]=cell_list[cells][k+1][i];
//		i++; printf("%d",cell_list[cells][k+1][0]);
//		}
//	}
//
//	for (i=0;i<cells+2;i++){
//		for (j=0;j<cells+2;j++){
//			printf("%d  ",cell_list[i][j][0]);
//		}
//		printf("\n");
//	}
//
//	printf("\n \n");
//	for (i=0;i<cells+2;i++){
//			for (j=0;j<cells+2;j++){
//				printf("%d  ",track_list[i][j]);
//			}
//			printf("\n");
//		}
	//loop over timestep
	for ( t = t_ini; t <= t_fin; t += timestep ){
		count += 1;
		if(count==savesteps){
			fprintf(values,"\n %.5f ",t);
			}
		//create cell list with space for ghost cells
		//int cell_list[cells+2][cells+2][zeros_in_cell]={0}, track_list[cells+2][cells+2]={0};

		int** track_list;

		// the following code allocates an array of int* of size n to p
		// p now points to an array of int* of size n
		// Each p[0 to n-1] points to int*
		track_list = (int**)malloc(cells+2 * sizeof(int*));

		// the following loop allocates an array of int of size m to p[0 to n-1]
		// every p[0 to n-1] now points to an array of int of size m
		// Each p[0 to n-1][0 to m-1] points to int
		for (i = 0; i < cells+2; ++i) {
		        track_list[i] = (int*)malloc(cells+2 * sizeof(int));
		}

		// the following code allocates an array of int** of size n to p
		// p now points to an array of int** of size n
		// Each p[0 to n-1] points to int**
		int***cell_list;
		cell_list = (int***)malloc(cells+2 * sizeof(int));

		// the following loop allocates an array of int* of size m to p[0 to n-1]
		// every p[0 to n-1] now points to an array of int* of size m
		// Each p[0 to n-1][0 to m-1] points to int*
		for (i = 0; i < cells+2; ++i) {
		        cell_list[i] = (int**)malloc(cells+2 * sizeof(int*));
		}

		// the following loop allocates an array of int of size l to p[0 to n-1][0 to m-1]
		// every p[0 to n-1][0 to m-1] nows points to an array of int of size l
		// Each p[0 to n-1][0 to m-1][0 to l-1] points to int
		for (i = 0; i < cells+2; ++i) {
		        for (j = 0; j < zeros_in_cell; ++j) {
		                cell_list[i][j] = (int*)malloc(zeros_in_cell * sizeof(int));
		        }
		}

		for (i = 0; i < cells+2; ++i) {
				for (j = 0; j < cells+2; ++j) {
					track_list[i][j] = 0;
				}
			}



		for ( k = 0; k < N; k++){  //loop over particles
			//periodic boundary conditions
			if (x[k] > L ){ x[k] = x[k] - L; }
			if (x[k] < 0 ){ x[k] = x[k] + L; }
			if (y[k] > L ){ y[k] = y[k] - L; }
			if (y[k] < 0 ){ y[k] = y[k] + L; }
			if(count==savesteps){
				fprintf(values,"%.5f  %.5f  %.5f ",x[k],y[k],phi[k]);
			}
				i = (int) x[k] / size_of_cell + 1;
				j = (int) y[k] / size_of_cell + 1;
				cell_list[i][j][track_list[i][j]] = k;
				if ( track_list[i][j] <= zeros_in_cell ){
					track_list[i][j] += 1;
				}
				if (track_list[i][j]==zeros_in_cell){
					break;
				}
		}
		if (track_list[i][j]==zeros_in_cell){
							printf( "ERROR, too many particles for 3d array! \n" );
							break;
						}
		if(count==savesteps){
			count = 0;
		}

		//fill ghost cells except corners
		for ( k = 0; k < cells ; k++){

				track_list[0][k+1]=track_list[cells][k+1];
				track_list[cells+1][k+1]=track_list[1][k+1];
				track_list[k+1][0]=track_list[k+1][cells];
				track_list[k+1][cells+1]=track_list[k+1][1];

				i=0; while ((int)cell_list[cells][k+1][i]!=0){
						cell_list[0][k+1][i]=cell_list[cells][k+1][i];i++;
				}
				i=0; while ((int)cell_list[1][k+1][i]!=0){
						cell_list[cells+1][k+1][i]=cell_list[1][k+1][i];i++;
				}
				i=0; while ((int)cell_list[k+1][cells][i]!=0){
						cell_list[k+1][0][i]=cell_list[k+1][cells][i];i++;
				}
				i=0; while ((int)cell_list[k+1][1][i]!=0){
						cell_list[k+1][cells+1][i]=cell_list[k+1][1][i];i++;
				}
			}
		//fill ghost cells in corners
			track_list[0][0]=track_list[cells][cells];
			track_list[0][cells+1]=track_list[cells][1];
			track_list[cells+1][0]=track_list[1][cells];
			track_list[cells+1][cells+1]=track_list[1][1];

			i=0; while ((int)cell_list[cells][cells][i]!=0){
					cell_list[0][0][i]=cell_list[cells][cells][i];i++;
			}
			i=0; while ((int)cell_list[cells][1][i]!=0){
					cell_list[0][cells+1][i]=cell_list[cells][1][i];i++;
			}
			i=0; while ((int)cell_list[1][cells][i]!=0){
					cell_list[cells+1][0][i]=cell_list[1][cells][i];i++;
			}
			i=0; while ((int)cell_list[1][1][i]!=0){
					cell_list[cells+1][cells+1][i]=cell_list[1][1][i];i++;
			}

		//loop over cells
		for ( i = 1; i <= cells; i++){      									//cell loop y
			for ( j = 1; j <= cells; j++) {										//cell loop x
				for (k=0; k < track_list[i][j]; k++ ){							//particles in ij-cell loop
					p1=cell_list[i][j][k];										//current particle (p1)
					wca_x = 0; wca_y=0;
//					for ( m = 0; m < 3; m++ ){									//loop over neighbour cells y
//						for ( n=0; n < 3; n++ ){								//loop over neighbour cells x
//							for ( p = 0; p < track_list[i-1+m][j-1+n]; p++ ){	//loop over particles in neighbour cell
//								p2=cell_list[i-1+m][j-1+n][p];			//pick particle 2 (p2)
//								if ( p1 != p2 ){					//ignore self
//									distance_inv = 1 / sqrt( (x[p1] - x[p2])*(x[p1] - x[p2]) + (y[p1] - y[p2])*(y[p1] - y[p2]) );
//									if ( 1 / distance_inv < cutoff){
//										dis_inv7 = distance_inv * distance_inv * distance_inv * distance_inv * distance_inv * distance_inv * distance_inv;
//										wca_x += epsilon_wca * ( 48.0 * sigma_12 * dis_inv7 * dis_inv7 - 24.0 * sigma_6 * dis_inv7 * distance_inv ) * (x[p1] - x[p2]);
//										wca_y += epsilon_wca * ( 48.0 * sigma_12 * dis_inv7 * dis_inv7 - 24.0 * sigma_6 * dis_inv7 * distance_inv ) * (y[p1] - y[p2]);
//									}
//								}
//							}
//						}
//					}
					x[p1] 		+= wca_x * timestep + Pe * cos(phi[p1]) * timestep +  normal_distribution(0,1) * sqrttimestep2;
					y[p1] 		+= wca_y * timestep + Pe * sin(phi[p1]) * timestep +  normal_distribution(0,1) * sqrttimestep2;
					phi[p1] 		+= normal_distribution(0,1) * sqrttimestep2;
				}
			}
		}




		//loop over particles
//		for ( j = 0; j < N; j++ ) {
//
//			if(count==savesteps){
//				fprintf(values,"%.5f  %.5f  %.5f ",x[j],y[j],phi[j]);
//			}
			//calculate forces within some distance
//			wca_x = 0; wca_y=0;
//			for ( k = 0; k < N && k != j; k++){
//				distance_inv = 1 / sqrt( (x[j] - x[k])*(x[j] - x[k]) + (y[j] - y[k])*(y[j] - y[k]) );
//				if ( 1 / distance_inv < cutoff){
//					dis_inv7 = distance_inv * distance_inv * distance_inv * distance_inv * distance_inv * distance_inv * distance_inv;
//					wca_x += epsilon_wca * ( 48 * sigma12 * dis_inv7 * dis_inv7 - 24 *sigma6 * dis_inv7 * distance_inv ) * (x[j] - x[k]);
//					wca_y += epsilon_wca * ( 48 * sigma12 * dis_inv7 * dis_inv7 - 24 *sigma6 * dis_inv7 * distance_inv ) * (y[j] - y[k]);
//				}
//			}

			//discretized differential equations
//			x[j] 		+= wca_x * timestep + Pe * cos(phi[j]) * timestep +  normal_distribution(0,1) * sqrttimestep2;
//			y[j] 		+= wca_y * timestep + Pe * sin(phi[j]) * timestep +  normal_distribution(0,1) * sqrttimestep2;
//			phi[j] 		+= normal_distribution(0,1) * sqrttimestep2;
//		}

	}

	fclose(values);
	//get runtime
	clock_t end = clock();
	time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
	if(time_spent<60){
	printf("Total runtime: %.1f seconds \n", time_spent);
	}
	if(time_spent>60){
		printf("Total runtime: %.1f minutes \n", time_spent/60.0);
		}

	//just for my sanity
	printf("\n \n");
	int cell_list1[cells+2][cells+2][zeros_in_cell]={0}, track_list1[cells+2][cells+2]={0};
	for ( k = 0; k < N; k++){
			i = (int) x[k] / size_of_cell + 1;
			j = (int) y[k] / size_of_cell + 1;
			cell_list1[i][j][track_list1[i][j]] = k;
			if ( track_list1[i][j] <= zeros_in_cell ){
				track_list1[i][j] += 1;
			}
			if (track_list1[i][j]==zeros_in_cell){
				//printf( "Error, too many particles in list!" );
			}
	}

	for ( k = 0; k < cells ; k++){

		track_list1[0][k+1]=track_list1[cells][k+1];
		track_list1[cells+1][k+1]=track_list1[1][k+1];
		track_list1[k+1][0]=track_list1[k+1][cells];
		track_list1[k+1][cells+1]=track_list1[k+1][1];

		i=0; while ((int)cell_list1[cells][k+1][i]!=0){
				cell_list1[0][k+1][i]=cell_list1[cells][k+1][i];i++;
		}
		i=0; while ((int)cell_list1[1][k+1][i]!=0){
				cell_list1[cells+1][k+1][i]=cell_list1[1][k+1][i];i++;
		}
		i=0; while ((int)cell_list1[k+1][cells][i]!=0){
				cell_list1[k+1][0][i]=cell_list1[k+1][cells][i];i++;
		}
		i=0; while ((int)cell_list1[k+1][1][i]!=0){
				cell_list1[k+1][cells+1][i]=cell_list1[k+1][1][i];i++;
		}
		//printf("%d",cell_list1[0][k][0]);
	}
	track_list1[0][0]=track_list1[cells][cells];
	track_list1[0][cells+1]=track_list1[cells][1];
	track_list1[cells+1][0]=track_list1[1][cells];
	track_list1[cells+1][cells+1]=track_list1[1][1];
	i=0; while ((int)cell_list1[cells][cells][i]!=0){
			cell_list1[0][0][i]=cell_list1[cells][cells][i];i++;
	}
	i=0; while ((int)cell_list1[cells][1][i]!=0){
			cell_list1[0][cells+1][i]=cell_list1[cells][1][i];i++;
	}
	i=0; while ((int)cell_list1[1][cells][i]!=0){
			cell_list1[cells+1][0][i]=cell_list1[1][cells][i];i++;
	}
	i=0; while ((int)cell_list1[1][1][i]!=0){
			cell_list1[cells+1][cells+1][i]=cell_list1[1][1][i];i++;
	}

	for (i=0;i<cells+2;i++){
		for (j=0;j<cells+2;j++){
			printf("%d  ",cell_list1[i][j][0]);
		}
		printf("\n");
	}

	printf("\n \n");
	for (i=0;i<cells+2;i++){
			for (j=0;j<cells+2;j++){
				printf("%d  ",track_list1[i][j]);
			}
			printf("\n");
		}
	free(x);
	free(y);
	free(phi);
	printf("\n %f \n %f \n %f",sigma_w,sigma_6,sigma_12);
return EXIT_SUCCESS;
}


