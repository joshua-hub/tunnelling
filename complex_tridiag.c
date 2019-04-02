#include<complex.h>
#include<stdlib.h>

int complex_tridiag(double complex diag[], double complex a[], double complex c[], double complex r[], double complex x[], int n)
{
        int j;
        double complex bet,*gam;

        gam = malloc(n*sizeof(double complex));
        if (cabs(diag[0]) == 0.0) return -1;    // Error
        x[0]=r[0]/(bet=diag[0]);
        for (j=1;j<n;j++) {
             gam[j]=c[j-1]/bet;
             bet=diag[j]-a[j]*gam[j];
             if (cabs(bet) == 0.0) return -1;   // Error
             x[j]=(r[j]-a[j]*x[j-1])/bet;
        }
        for (j=n-2;j>=0;j--)
             x[j] -= gam[j+1]*x[j+1];
        free(gam); 
        return 0;
}
