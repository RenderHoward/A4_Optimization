#include <glib.h>
#include <math.h>
#include "jstate.h"

// gcc -lm test_j.c -o test_j

// Calculated constants for duration of program
double pi, fptildemin, aX, gX, log_alphamin, log_gammamin, sigma_a, sigma_b, log_aC, log_gC;

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

    log_alphamin  = log(aC  * pow(fptildemin, aX));
    log_gammamin  = log(gC  * pow(fptildemin, gX));

    sigma_a = saC * pow(fptildemin, saX);
    sigma_b = sbC * pow(fptildemin, sbX);

    g2div16pi4 = pow(g, 2) * pow((2*pi), -4);

    log_aC = log(aC);

    log_gC = log(gC);
}

double step = 0.1;

AGvals *cache;
int ind_fptilde;

int main(int argc, char** argv )
{
    if( argc > 1 )
        step = atof(argv[1]);

    double S = 0, f, fp, fptilde,
            accum = 0.0;

    set_constants();

    JState state;

    init_state( &state, -5, 0, 0);

    int count = (int)floor(10/step);

    //  Cache of intermediate values for re-use
    cache = malloc(sizeof(AGvals) * (count +1));

    for( int i=0; i<=count; i++)
        cache[i].valid = FALSE;

    f = -5;
    for  (int ind_f = 0 ; ind_f <= count; ind_f++){  fp = 0;
      for (int ind_fp = 0; ind_fp <= count; ind_fp++){  fptilde = 0;
          for (ind_fptilde = 0; ind_fptilde <= count; ind_fptilde++){

            S = fast_function_j( &state, f, fp, fptilde );

            // Just for preventing auto optimization from deleting everything.
            accum += S;
#ifdef TEST // Dump subsampled points to verify faster code still matches output
            g_print("%.15f \n", S);
#endif
            fptilde += step;
       }fp += step;
      }f += step;
    }

    free(cache);

    return 0;
}

double fast_function_j( JState *state, double f, double fp, double fptilde )
{
    set_parameters(state, f, fp, fptilde);

    return func_j(state);
}

double func_j( JState *state )
{
    double f5 = state->f * state->f * state->f * state->f * state->f;
    double ret =
            g2div16pi4 *
            exp(
                state->log_alpha +
                state->exp1arg +
                state->log_gamma * state->exp_exp2arg
                )/f5;
    return ret;
}

void update_exp1(JState *state)  {  state->exp1arg = -1.25 * pow((state->f/state->fp),-4);  }

void update_exp2(JState *state)
{
    state->exp_exp2arg = exp( -0.5 * pow((state->f - state->fp)/(state->sigma*state->fp), 2) );
}

void update_sigma(JState *state)
{
    if( state->fptilde > fptildemin )
    {
        AGvals *vals = cache + ind_fptilde;

        if( !vals->valid )
        {
            vals->log_alpha = log_aC + log(state->fptilde) * aX;
            vals->log_gamma = log_gC + log(state->fptilde) * gX;
            vals->sigma_a = saC * pow(state->fptilde, saX);
            vals->sigma_b = sbC * pow(state->fptilde, sbX);
            vals->valid = TRUE;
        }

        if (state->f_lte_fp)
            state->sigma = vals->sigma_a;
        else
            state->sigma = vals->sigma_b;
    }
    else
        state->sigma = (state->f_lte_fp) ? sigma_a : sigma_b;
}

void AddToCache(JState *state)
{
    AGvals *vals = cache + ind_fptilde;

    vals->log_alpha = state->log_alpha;
    vals->log_gamma = state->log_gamma;
    vals->sigma_a = saC * pow(state->fptilde, saX);
    vals->sigma_b = sbC * pow(state->fptilde, sbX);
    vals->valid = TRUE;
}

void update_ag(JState *state)
{
    if( state->fptilde > fptildemin )
    {
        AGvals *vals = cache + ind_fptilde;

        if( vals->valid )
        {
            state->log_alpha = vals->log_alpha;
            state->log_gamma = vals->log_gamma;
        }
        else
        {
            state->log_alpha = log_aC + log(state->fptilde) * aX;
            state->log_gamma = log_gC + log(state->fptilde) * gX;

            AddToCache(state);
        }
    }
    else
    {
        state->log_alpha = log_alphamin;
        state->log_gamma = log_gammamin;
    }
}

void set_parameters( JState *state,  double f, double fp, double fptilde )
{
    gboolean ag_dirty, sigma_dirty, exp1_dirty, exp2_dirty;

    ag_dirty = sigma_dirty = exp1_dirty = exp2_dirty = FALSE;

    if( state->f != f )
        exp2_dirty = exp1_dirty = TRUE;
    if( state->fp != fp )
        exp2_dirty = exp1_dirty = TRUE;
    if( state->f_lte_fp != (f <= fp) )
        sigma_dirty = exp2_dirty = TRUE;
    if( state->fptilde != fptilde )
        ag_dirty = sigma_dirty = exp2_dirty = TRUE;

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

