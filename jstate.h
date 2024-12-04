#ifndef JSTATE_H
#define JSTATE_H
#include <glib.h>
typedef struct
{
    double f, fp, fptilde;

    gboolean f_lte_fp;

    double log_alpha, log_gamma, sigma;

    int fpt_ind;

    double exp1arg, exp_exp2arg;
} JState;

typedef struct
{
    gboolean valid;
    double log_alpha, log_gamma, sigma_a, sigma_b, fptilde;
} AGvals;

double func_j( JState *state );

void update_exp1(JState *state);

void update_exp2(JState *state);

void update_sigma(JState *state);

void update_ag(JState *state);

void set_parameters( JState *state,  double f, double fp, double fptilde, int ind_fpt );

void init_state( JState *state, double f, double fp, double fptilde );

double fast_function_j( JState *state, double f, double fp, double fptilde );

#endif // JSTATE_H
