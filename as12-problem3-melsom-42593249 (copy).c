// phys3071 as12 melsom 42593249

/*_____________________________________________________________________________
Description: 

Inputs: 

Calculations: 

outputs: 

Compiled as
        gcc as12-problem3-melsom-42593249.c complex_tridiag.c -o as12 -Wall -lm
_____________________________________________________________________________*/

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<complex.h>

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
  double delta_x;
  double delta_t;
  double t_final;
  double x_initial;
  double x_final;
  double b_nought;
  double t_out;
  double complex r;
  int N, J;
  int count;
  int scan_temp;
  double *potential;
  double complex *diag;
  double complex *under_diag;
  double complex *above_diag;
  double complex *psi_old;
  double complex *psi_new;
  double complex *f_of_psi_old;
  int k= 2;
  double x_position;
  int temp;
  
  char fname[]= "temp.dat";
  FILE *dat;
  
  printf("This program will solve the wavefunction of a particle in a harmonic"
         " potential wih an potential barrier.\n");
  printf("First I will need some values from you.");
  
  /*
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
  */
  
  // START TEMP COMMENTS
  // Temporary hard coded value assignments.
  printf("\n");
  x_initial = -20.0;
  x_final= 20.0;
  t_final = 10.0;
  delta_t = 0.0025;
  t_out= 0.05;
  J = 2000;
  scan_temp= 2;
  b_nought= 0.0;
  // END TEMP COMMENTS
  
  // Checks that all imputted values are OK.
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
	
	delta_x = (x_final- x_initial)/ (1.0* J);
  N= (int) ceil(t_final/ delta_t);
  
  potential= malloc(J* sizeof(double));
  diag= malloc(J* sizeof(complex double));
  psi_old= malloc(J* sizeof(complex double));
  psi_new= malloc(J* sizeof(complex double));
  under_diag= malloc(J* sizeof(complex double));
  above_diag= malloc(J* sizeof(complex double));
  f_of_psi_old= malloc(J* sizeof(complex double));
  
  // Assign a variable that will be used often in future calculations.
  r= I* delta_t/ (delta_x* delta_x);
  
  //NEED TO LOOP OVER TIME FROM HERE>!$%&!$%&!$%&!$%&!$%&!$%!&!$%&!$%&!$!%&!$&%!
  
  // Open file to write to 
  if ((dat=fopen(fname,"w")) == NULL) {
    fprintf(stderr, "ERROR: Could not open file name: %s\n", fname);
    exit(EXIT_FAILURE); 
  }
  
  // Populate arrays
  for (count= 0; count< J; count++) {
    potential[count]= pow((x_initial+ count* delta_x), 2.0);
    x_position= x_initial+ count* delta_x;
    psi_old[count]= psi_of_x(k, x_position);
    above_diag[count]= -r;
    under_diag[count]= -r;
    diag[count]= 
        2.0+ I* delta_t* potential[count]+ 2.0* r;
  }

  // Force boundaries to zero.
  psi_old[0]= 0.0;
  psi_old[J]= 0.0;

  // Calculates the RHS of the linear equation which is a function of psi_old
  for(count= 1; count< J-1; count++) {
    f_of_psi_old[count]= (2.0- I* delta_t* potential[count]- 2* r)* 
        psi_old[count]+ r* (psi_old[count-1] + psi_old[count+ 1]);
  }
  // Has to also assign a few individual values because of boundaries.
  f_of_psi_old[0]= (2.0- I* delta_t* potential[0]- 2* r)* 
                  psi_old[count]+ r* psi_old[1];
  f_of_psi_old[J-1]= (2.0- I* delta_t* potential[J-1]- 2* r)* 
                  psi_old[J-1]+ r* psi_old[J-2];
  
  // Solve for psi_new using the complex tridiag program.
  if ((temp= complex_tridiag(diag, under_diag, above_diag, f_of_psi_old, 
       psi_new, J)) != 0) {
    fprintf(stderr, "ERROR: failed matrix inversion.\n");
    exit(EXIT_FAILURE); 
  }
  
  // Replace psi_old with psi_new for the next step.
  for (count= 0; count< J; count++) {
    psi_old[count]= psi_new[count];
  }
  
  fclose(dat);
  
  //NEED TO LOOP OVER TIME TO HERE>!$%&!$%&!$%&!$%&!$%&!$%!&!$%&!$%&!$!%&!$&%!!$
  
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
