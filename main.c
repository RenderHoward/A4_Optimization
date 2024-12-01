#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "jstate.h"
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
double g2div16pi4;

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

    g2div16pi4 = pow(g, 2) * pow((2*pi), -4);
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

double step = 0.1;

int main(int argc, char** argv )
{
    if( argc > 1 )
        step = atof(argv[1]);

    double S, f, fp, fptilde,
            accum = 0.0;  // Just for preventing auto optimization from deleting everything.

    set_constants();

    JState state;

    init_state( &state, -5, 0, 0);

    for  (f = -5.; f <= 5.; f += step){
      for (fp = 0.; fp <= 10.; fp += step) {
        for (fptilde = 0.; fptilde <= 10.; fptilde += step){

         set_parameters(&state, f, fp, fptilde);

         S = func_j(&state);

         accum += S;
         // printf("%.15f \n", S);
        }
      }
    }

    printf( "%f", accum );

    return 0;
}


double func_j( JState *state )
{
    return state->alpha * pow(g, 2) * pow((2*pi), -4) * pow(state->f,-5) * exp(state->exp1arg) * pow(state->gamma, exp(state->exp2arg));
}

void update_exp1(JState *state)  {  state->exp1arg = -1.25 * pow((state->f/state->fp),-4);  }

void update_exp2(JState *state)
{
    state->exp2arg = -0.5 * pow((state->f - state->fp)/(state->sigma*state->fp), 2);
}

void update_sigma(JState *state)
{
    if( state->fptilde > fptildemin )
    {
        if (state->f_lte_fp)
            state->sigma = saC * pow(state->fptilde, saX);
        else
            state->sigma = sbC * pow(state->fptilde, sbX);
    }
    else
        state->sigma = (state->f_lte_fp) ? sigma_a : sigma_b;
}

void update_ag(JState *state)
{
    if( state->fptilde > fptildemin )
    {
        state->alpha = aC * pow(state->fptilde, aX);
        state->gamma = gC * pow(state->fptilde, gX);
    }
    else
    {
        state->alpha = alphamin;
        state->gamma = gammamin;
    }
}

void set_parameters( JState *state,  double f, double fp, double fptilde )
{
    bool ag_dirty, sigma_dirty, exp1_dirty, exp2_dirty;

    ag_dirty = sigma_dirty = exp1_dirty = exp2_dirty = false;

    if( state->f != f )
        exp2_dirty = exp1_dirty = true;
    if( state->fp != fp )
        exp2_dirty = exp1_dirty = true;
    if( state->f_lte_fp != (f <= fp) )
        sigma_dirty = exp2_dirty = true;
    if( state->fptilde != fptilde )
        ag_dirty = sigma_dirty = exp2_dirty = true;

    state->f = f;
    state->fp = fp;
    state->fptilde = fptilde;
    state->f_lte_fp = (f <= fp);

    if(ag_dirty)
        update_ag(state);
    if(exp1_dirty)
        update_exp1(state);
    if(sigma_dirty)
        update_sigma(state);
    if(exp2_dirty)
        update_exp2(state);
}

void init_state( JState *state, double f, double fp, double fptilde )
{
    state->f = f, state->fp = fp, state->fptilde = fptilde;

    state->f_lte_fp = f <= fp;

    update_ag(state);

    update_sigma(state);

    update_exp1(state);

    update_exp2(state);
}

