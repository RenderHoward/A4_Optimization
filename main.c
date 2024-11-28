#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

// gcc -lm test_j.c -o test_j

// Calculated constants for duration of program
double pi, fptildemin, aX, gX, alphamin, gammamin, sigma_a, sigma_b;

// Compile time constants for all time.
double a  = 0.0081;
double b  = 0.6;
double g  = 9.807;
double saC = 0.0547;
double sbC = 0.0783;
double sbX = 0.16;
double saX = 0.32;
double gC = 5.87;
double aC = 0.0317;


void set_constants()
{
    pi = 4.*atan(1.);

    fptildemin = (1.0/2.0/pi) * pow((4.0 * b / 5.0), (1.0/4.0));

    aX  = (log(a)-log(aC))/log(fptildemin);

    gX  = -log(gC)/log(fptildemin);

    alphamin  = aC  * pow(fptildemin, aX);
    gammamin  = gC  * pow(fptildemin, gX);

    sigma_a = saC * pow(fptildemin, saX);
    sigma_b = sbC * pow(fptildemin, sbX);
}

double function_j(double f, double fp, double fptilde) {

   double alpha, gamma, exp1arg, sigma;

   exp1arg = -1.25 * pow((f/fp),-4);

   double fpt = MAX(fptilde, fptildemin);

   if( fptilde > fptildemin )
   {
       alpha   = aC  * pow(fpt, aX);
       gamma   = gC  * pow(fpt, gX);

       if (f <= fp)
           sigma = saC * pow(fpt, saX);
       else
           sigma = sbC * pow(fpt, sbX);
   }


   double exp2arg = -0.5 * pow((f-fp)/(sigma*fp), 2);

   double S = alpha * pow(g, 2) * pow((2*pi), -4) * pow(f,-5) * exp(exp1arg) * pow(gamma, exp(exp2arg));

   return S;
}

int main(int argc, char** argv )
{
    if( argc > 1 )
        step = atof(argv[1]);

    double S, f, fp, fptilde,
            accum = 0.0;  // Just for preventing auto optimization from deleting everything.

    set_constants();
    for (f = -5.; f <= 5.; f += step) {
      for (fp = 0.; fp <= 10.; fp += step) {
        for (fptilde = 0.; fptilde <= 10.; fptilde += step) {
          S = function_j(f, fp, fptilde);
          accum += S;
          printf("%.15f \n", S);
        }
      }
    }

    printf( "%f", accum );

    return 0;
}

