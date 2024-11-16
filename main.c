#include <stdio.h>
#include <math.h>
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

// gcc -lm test_j.c -o test_j

double function_j(double f, double fp, double fptilde) {

   double a  = 0.0081;
   double b  = 0.6;
   double g  = 9.807;
   double pi = 4.*atan(1.);

   double fptildemin = (1.0/2.0/pi) * pow((4.0 * b / 5.0), (1.0/4.0));

   double gC = 5.87;
   double aC = 0.0317;

   double aX  = (log(a)-log(aC))/log(fptildemin);
   double gX  = -log(gC)/log(fptildemin);

   double saC = 0.0547;
   double saX = 0.32;

   double sbC = 0.0783;
   double sbX = 0.16;

   double fpt = MAX(fptilde, fptildemin);

   double alpha   = aC  * pow(fpt, aX);
   double gamma   = gC  * pow(fpt, gX);
   double sigma_a = saC * pow(fpt, saX);
   double sigma_b = sbC * pow(fpt, sbX);
   double exp1arg = -1.25 * pow((f/fp),-4);
   double sigma   = (f <= fp) * sigma_a + (f > fp) * sigma_b;

   double exp2arg = -0.5 * pow((f-fp)/(sigma*fp), 2);

   double S = alpha * pow(g, 2) * pow((2*pi), -4) * pow(f,-5) * exp(exp1arg) * pow(gamma, exp(exp2arg));

   return S;
}

main()
{
    double S, f, fp, fptilde;

    for (f = -5.; f <= 5.; f += 0.01) {
      for (fp = 0.; fp <= 10.; fp += 0.01) {
        for (fptilde = 0.; fptilde <= 10.; fptilde += 0.01) {
          S = function_j(f, fp, fptilde);
          printf("%.15f \n", S);
        }
      }
    }
}

