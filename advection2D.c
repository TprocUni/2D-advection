/*******************************************************************************
2D advection example program which advects a Gaussian u(x,y) at a fixed velocity



Outputs: initial.dat - inital values of u(x,y) 
         final.dat   - final values of u(x,y)

         The output files have three columns: x, y, u

         Compile with: gcc -o advection2D -std=c99 advection2D.c -lm

Notes: The time step is calculated using the CFL condition

********************************************************************************/

/*********************************************************************
                     Include header files 
**********************************************************************/

#include <stdio.h>
#include <math.h>

/*including openMP*/
#include <omp.h>

#include<stdlib.h>

/*calc_x_shear prototype*/
float calc_x_shear (float height);

/*********************************************************************
                      Main function
**********************************************************************/

int main(){

  /*time stuff*/
  double start_time, end_time, elapsed_time;
  start_time = omp_get_wtime();

  /* Grid properties */
  const int NX=1000;    // Number of x points
  const int NY=1000;    // Number of y points
  const float xmin=0.0; // Minimum x value
  const float xmax=30.0; // Maximum x value
  const float ymin=0.0; // Minimum y value
  const float ymax=30.0; // Maximum y value
  
  /* Parameters for the Gaussian initial conditions */
  const float x0=3.0;                    // Centre(x)
  const float y0=15.0;                    // Centre(y)
  const float sigmax=1.0;               // Width(x)
  const float sigmay=5.0;               // Width(y)
  const float sigmax2 = sigmax * sigmax; // Width(x) squared
  const float sigmay2 = sigmay * sigmay; // Width(y) squared

  /* Boundary conditions */
  const float bval_left=0.0;    // Left boudnary value
  const float bval_right=0.0;   // Right boundary value
  const float bval_lower=0.0;   // Lower boundary
  const float bval_upper=0.0;   // Upper bounary
  
  /* Time stepping parameters */
  const float CFL=0.9;   // CFL number 
  const int nsteps=800; // Number of time steps

  /* Velocity */
  const float velx=calc_x_shear(ymax); // Velocity in x direction
  const float vely=0.0; // Velocity in y direction
  
  /* Arrays to store variables. These have NX+2 elements
     to allow boundary values to be stored at both ends */
  float x[NX+2];          // x-axis values
  float y[NX+2];          // y-axis values
  float u[NX+2][NY+2];    // Array of u values
  float dudt[NX+2][NY+2]; // Rate of change of u

  /*data structure for vertically averaged distribution*/
  float *averaged = (float *)calloc(NX, sizeof(float));

  float x2;   // x squared (used to calculate iniital conditions)
  float y2;   // y squared (used to calculate iniital conditions)
  
  /* Calculate distance between points */
  float dx = (xmax-xmin) / ( (float) NX);
  float dy = (ymax-ymin) / ( (float) NY);
  
  /* Calculate time step using the CFL condition */
  /* The fabs function gives the absolute value in case the velocity is -ve */
  float dt = CFL / ( (fabs(velx) / dx) + (fabs(vely) / dy) );
  
  /*** Report information about the calculation ***/
  printf("Grid spacing dx     = %g\n", dx);
  printf("Grid spacing dy     = %g\n", dy);
  printf("CFL number          = %g\n", CFL);
  printf("Time step           = %g\n", dt);
  printf("No. of time steps   = %d\n", nsteps);
  printf("End time            = %g\n", dt*(float) nsteps);
  printf("Distance advected x = %g\n", velx*dt*(float) nsteps);
  printf("Distance advected y = %g\n", vely*dt*(float) nsteps);

  /*** Place x points in the middle of the cell ***/ 
/* LOOP 1 */
  #pragma omp parallel for shared(x, dx)
  for (int i=0; i<NX+2; i++){
    x[i] = ( (float) i - 0.5) * dx;
  }

  /*** Place y points in the middle of the cell ***/
  /* LOOP 2 */
  #pragma omp parallel for shared(y, dy)
  for (int j=0; j<NY+2; j++){
    y[j] = ( (float) j - 0.5) * dy;
  }

  /*** Set up Gaussian initial conditions ***/
  /* LOOP 3 */

  #pragma omp parallel for private(x2,y2) shared(x,y,u)
  for (int i=0; i<NX+2; i++){
    for (int j=0; j<NY+2; j++){
      x2      = (x[i]-x0) * (x[i]-x0);
      y2      = (y[j]-y0) * (y[j]-y0);
      u[i][j] = exp( -1.0 * ( (x2/(2.0*sigmax2)) + (y2/(2.0*sigmay2)) ) );
    }
  }

  /*** Write array of initial u values out to file ***/
  FILE *initialfile;
  initialfile = fopen("initial.dat", "w");

  /* LOOP 4 */

  /*Collapse directive results in program being ~0.3 seconds slower
    when using the initial values (part 2.1).
    Ordered directive results in program being ~0.05 seconds slower
    when using the initial valiues (part 2.1)
    Although further investigation could be conducted into seeing 
    how they scale, ordered will be used for now.
    It is worth noting that under the initial conditions, the 
    serialised version of the code performs best (LOOP 4 exclusive). 
    Additionally, the initial.dat files are different
                      							 */
  /*
  ANSWER:
  Since we are writing to a file, if we parallelise the loop the writes
  will be out of order. This is an output dependancy.
  */

  for (int i=0; i<NX+2; i++){
    for (int j=0; j<NY+2; j++){
      fprintf(initialfile, "%g %g %g\n", x[i], y[j], u[i][j]);
    }
  }
  fclose(initialfile);
  
  /*** Update solution by looping over time steps ***/
  /* LOOP 5 */

  /*This loop cannot be parallelised since m is a flow dependant value, as it is 
    used by some of the subsequent loops, */

  for (int m=0; m<nsteps; m++){
    
    /*** Apply boundary conditions at u[0][:] and u[NX+1][:] ***/
    /* LOOP 6 */
    #pragma omp parallel for shared(u)
    for (int j=0; j<NY+2; j++){
      u[0][j]    = bval_left;
      u[NX+1][j] = bval_right;
    }

    /*** Apply boundary conditions at u[:][0] and u[:][NY+1] ***/
    /* LOOP 7 */
    #pragma omp parallel for shared(u)
    for (int i=0; i<NX+2; i++){
      u[i][0]    = bval_lower;
      u[i][NY+1] = bval_upper;
    }
    
    /*** Calculate rate of change of u using leftward difference ***/
    /* Loop over points in the domain but not boundary values */
    /* LOOP 8 */
    #pragma omp parallel for shared(u,dudt,dy,dx)
    for (int i=1; i<NX+1; i++){
      for (int j=1; j<NY+1; j++){
	dudt[i][j] = -calc_x_shear(y[j]) * (u[i][j] - u[i-1][j]) / dx
	            - vely * (u[i][j] - u[i][j-1]) / dy;
      }
    }
    
    /*** Update u from t to t+dt ***/
    /* Loop over points in the domain but not boundary values */
    /* LOOP 9 */
    #pragma omp parallel for shared(u,dudt,dt,averaged)
    for	(int i=1; i<NX+1; i++){
      for (int j=1; j<NY+1; j++){
	u[i][j] = u[i][j] + dudt[i][j] * dt;
        /*sum all vertical values*/
        averaged[i-1] += u[i][j];
      }
      /*divide by NY for average*/
      averaged[i-1] /= NY;
    }
    
  } // time loop
  
  /*** Write array of final u values out to file ***/
  FILE *finalfile;
  finalfile = fopen("final.dat", "w");
  /* LOOP 10 */
  /*
  We cannot parallelise this loop otherwise an output dependency 
  will occur.
  */
  for (int i=0; i<NX+2; i++){
    for (int j=0; j<NY+2; j++){
      fprintf(finalfile, "%g %g %g\n", x[i], y[j], u[i][j]);
    }
  }
  fclose(finalfile);

  /*write vertically averaged distribution stuff to file*/
  FILE *vertically_averaged;
  vertically_averaged = fopen("vertically_averaged.dat", "w");
  /*loop through all values in verticallt averaged data*/
  for (int i = 0; i < NX; i++)
  {
    fprintf(vertically_averaged, "%g %g\n", x[i], averaged[i]);
  }
  fclose(vertically_averaged);
  /*release space taken up*/
  free(averaged);





  /*other time stuff*/
  end_time = omp_get_wtime();
  elapsed_time = end_time - start_time;
  printf("Elapsed time  = %f \n", elapsed_time);
  return 0;
}

/*Calc horizontal shear using logarithmic profile*/

float calc_x_shear (float height)
{

  const float frict = 0.2;	//Friction 
  const float rough = 1.0;	// Roughness length
  const float VK = 0.41;	// Von Karmans constant

  float vel_x = 0;

  //check if height is greater than the roughness length
  if (height > rough)
  {
    //do shear equation
    vel_x = (frict/VK) * logf(height/rough);
  }
  return vel_x;

}
 

/* End of file ******************************************************/
