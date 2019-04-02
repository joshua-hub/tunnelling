// phys3071 as12 melsom 42593249

/*_____________________________________________________________________________
Description: This program solves the wavefunction of a particle inside of a
harmonic potential with a central barrier. 

Inputs: The user is required to enter the Barrier beight, the lower and upper
x values, the number of spatial steps, the output time steps, and the final 
calculation time. Within variable declaration, the energy of the particle
can be modified by changing k between the values of 0 and 4 (inclusive).

Calculations: The program uses the crank-nicolson method to solve each 
time step iteration. The initial conditions call upon a fcuntion to populate
across the x positions.

outputs: The program outputs to terminal the current time, the Norm of the 
solution, the expectation value and the standard deviation. It also creates a
file (which is named in the variable declarations) which contains the 

Compiled as
        gcc as12-problem3-melsom-42593249.c complex_tridiag.c -o as12 -Wall -lm
_____________________________________________________________________________*/

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<complex.h>
#include<string.h>

// Function Prototypes ________________________________________________________

// complex diag is main diagonal, complex a is the diagonal below the main, 
// complex c is the diagonal above the main, r is the vector of known values, 
// complex x is the vector to be solved, int n is the length height of the 
// matrix that these values are from.
int complex_tridiag(double complex diag[], double complex a[], 
                    double complex c[], double complex r[], 
                    double complex x[], int n);

// Equation 4 from assignment 12 sheet
double psi_of_x (int k, double x);

// Returns the kth Hermite polynomial evaluated at x.
double hermite  (int k, double x);

// Returns k! with recursive function.
int fact(int k);

// Global Variables ___________________________________________________________


// Begin main function ________________________________________________________
int main () {
  double delta_x, x_initial, x_final; // Values for the x positions
  
  double delta_t, t_final, t_position; // Values for the time positions.
  double t_out, t_output;
  int time_count;

  int N, J; // The maximum values in spatial and temporal positions.
  int count; // Generic counter
  int scan_temp, temp; // temporary variables.
  int count2= 1;
  // Arrays for storing the position and potential at that position. 
  double *potential,*x_position;
  double complex *diag, *under_diag, *above_diag;  // Matrix tridiagonal values
  double complex *psi_old, *psi_new, *f_of_psi_old; // r, x values for tridiag
  double complex r; // Variable to store a common factor in calculations
  
  double norm; // Wavefunctin values
  double x_expect;
  double x_stand_dev;
  double psi_mod_squ;
  
  double b_nought; // potential function value
  int k= 0; // Can change which Hermite polynomial to use
  char fname[]= ".dat"; // For file output.
  char fnamenum[100];  // Beginning of file name.
  FILE *dat;
  
  printf("This program will solve the wavefunction of a particle in a harmonic"
         " potential wih an potential barrier.\n");
  printf("First I will need some values from you.");
  
  // User entry of potential barrier height.
  printf("\nPlease enter the height of the potential barrier: ");
  scanf("%lf", &b_nought);
  
  // User entry of the xmin and xmax values
  printf("Please enter the x range, x_0 and x_e as a pair of numbers: ");
  scan_temp = scanf("%lf%*[ ,;\t:]%lf",&x_initial, &x_final);
   	
	// User entry of number of spatial grid points.
  printf("Please enter the number of spatial grid positions: ");
  scanf("%i", &J);
  
  // User entry of size of time step delta_t
  printf("Please enter the size of time step, delta t: ");
  scanf("%lf", &delta_t);
  
  // User entry of output time steps.
  printf("Please enter the output time step size, delta t_out: ");
  scanf("%lf", &t_out);
  
  // User entry of final time 
  printf("Please enter the final time that you wish to calculate for: ");
  scanf("%lf", &t_final);
   
  // Checks that all inputted values are OK.
  if ((scan_temp!= 2)
   || (b_nought< 0.0)
   || (x_final< x_initial)
   || (J< 0)
   || (delta_t< 0.0)
   || (t_out < delta_t)
   || (t_final< t_out)) {
      fprintf(stderr, "ERROR: Invalid entry, Please start again.\n");
		  exit(EXIT_FAILURE);
	}
	
	// Calculate spatial stepsize, temporal final step, and time output steps.
	delta_x = (x_final- x_initial)/ (1.0* J);
  N= (int) ceil(t_final/ delta_t);
  t_output= t_out;
  
  // Reserve memory for arrays.
  potential= malloc(J* sizeof(double));
  x_position= malloc(J* sizeof(double));
  diag= malloc(J* sizeof(complex double));
  psi_old= malloc(J* sizeof(complex double));
  psi_new= malloc(J* sizeof(complex double));
  under_diag= malloc(J* sizeof(complex double));
  above_diag= malloc(J* sizeof(complex double));
  f_of_psi_old= malloc(J* sizeof(complex double));
  
  // Assign a variable that will be used often in future calculations.
  r= I* delta_t/ (delta_x* delta_x);
 
  // Populate the x positions the potential at x positions, 
  // the initial conditions and the tridiag values.
  for (count= 0; count< J; count++) {
    x_position[count]= x_initial+ count* delta_x;
    potential[count]= b_nought* exp(-16.0* pow(x_position[count], 2.0)) + 
              pow(x_position[count], 2.0);
    psi_old[count]= psi_of_x(k, (x_position[count]- 4.0));
    above_diag[count]= -r;
    under_diag[count]= -r;
    diag[count]= 
        2.0+ I* delta_t* potential[count]+ 2.0* r;
  }

  // Force boundaries to zero.
  psi_old[0]= 0.0;
  psi_old[J]= 0.0;
  
  // Iterate over each time step until t_final
  for (time_count= 0; time_count< N; time_count++) {
    
    // Calculates the RHS of the linear equation which is a function of psi_old
    for(count= 1; count< J-1; count++) {
      f_of_psi_old[count]= (2.0- I* delta_t* potential[count]- 2* r)* 
          psi_old[count]+ r* (psi_old[count-1] + psi_old[count+ 1]);
    }
    // Has to also assign a few individual values because of boundaries.
    f_of_psi_old[0]= (2.0- I* delta_t* potential[0]- 2* r)* 
                    psi_old[0]+ r* psi_old[1];
    f_of_psi_old[J-1]= (2.0- I* delta_t* potential[J-1]- 2* r)* 
                    psi_old[J-1]+ r* psi_old[J-2];
    
    // Solve for psi_new using the complex tridiag program.
    // Returns error if tridiag fails and ends program.
    if ((temp= complex_tridiag(diag, under_diag, above_diag, f_of_psi_old, 
         psi_new, J)) != 0) {
      fprintf(stderr, "ERROR: failed matrix inversion.\n");
      exit(EXIT_FAILURE); 
    }
    
    // Replace psi_old with psi_new for the next step.
    for (count= 0; count< J; count++) {
      psi_old[count]= psi_new[count];
    }
        
    // Initialise the t_output and it advances every time an output is reached.
    t_position= time_count* delta_t;
    // Output the values to terminal.
    if (t_position >= t_output) {
      // Calculate the norm, expectation of x and the standard deviation.
      // First initialise values and then add the sums.
      norm= 0.0;
      x_expect= 0.0;
      x_stand_dev= 0.0;
      for (count= 0; count< J; count++) {
        // Calculate the sum Norm(t_n)
        norm += delta_x* (pow(creal(psi_new[count]), 2.0)+ 
                pow(cimag(psi_new[count]), 2.0));
        // Calculate the sum <x> (t_n)     
        x_expect += delta_x* (x_position[count])* 
                (pow(creal(psi_new[count]), 2.0)+ 
                pow(cimag(psi_new[count]), 2.0));
        // Calculates the sum <x^2>, need to finish off calc outside of loop    
        x_stand_dev += delta_x* pow(x_position[count], 2.0)* 
                (pow(creal(psi_new[count]), 2.0)+ 
                pow(cimag(psi_new[count]), 2.0));        
      }
      
      // Finally calculate the end value of standard deviation.
      x_stand_dev = sqrt(x_stand_dev- pow(x_expect, 2.0));
           
      // output to terminal
      printf("%lf\t%lf\t%lf\t%lf\n", t_position, norm, x_expect, x_stand_dev);
      

      sprintf(fnamenum,"%03d",count2);  // Assign integer to char - middle of file name.
      strcat(fnamenum, fname);// Concatenate first 2 of 3 strings
      // Open file to write to 
      if ((dat=fopen(fnamenum,"w")) == NULL) {
        fprintf(stderr, "ERROR: Could not open file name: %s\n", fnamenum);
        exit(EXIT_FAILURE); 
      }
      for (count= 0; count< J; count++) {
        // output to file
        psi_mod_squ= pow(creal(psi_new[count]), 2.0)+ 
                   pow(cimag(psi_new[count]), 2.0); 
        fprintf(dat, "%lf\t%lf\n", x_position[count], psi_mod_squ);
      }    
      fclose(dat);
      
      // Increment to next t_out value.
      t_output += t_out;
      count2++;
    }
  }
  
  return (EXIT_SUCCESS);
}

// Functions __________________________________________________________________

// Solves equation 4 from assignment 12 sheets.
double psi_of_x (int k, double x) {
  return ((1.0/(sqrt(pow(2.0, k)* 1.0* fact(k)* sqrt(M_PI))))* 
          exp(-0.5* pow(x, 2.0))* hermite(k, x));
}

// returns the kth hermite polynomial evaluated at x.
double hermite(int k,  double x) {
  double H;  
  // Allows for the first 5 Hermite polynomials to be calculated.
  if (k== 0) (H= 1.0);
  else if (k== 1) (H= 2.0* x);
  else if (k== 2) (H= 4.0* pow(x, 2.0)- 2);
  else if (k== 3) (H= 8.0* pow(x, 3.0)- 12* x);
  else if (k== 4) (H= 16.0* pow(x, 4.0) - 48.0* pow(x,2.0)+ 12.0);
  
  return (H);
}

// Return k! using recursive function.
int fact(int k) {
  
  if (k == 0) 
    return (1);
  else 
    return (k* fact(k- 1));
}
