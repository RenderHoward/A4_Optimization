#ifndef JSTATE_H
#define JSTATE_H
#include <glib.h>
typedef struct
{
    double f, fp, fptilde;

    gboolean f_lte_fp;

    double alpha, gamma, sigma;

    double exp1arg, exp2arg;
} JState;

typedef struct
{
    gboolean valid;
    double alpha, gamma, sigma_a, sigma_b;
} AGvals;

double func_j( JState *state );

void update_exp1(JState *state);

void update_exp2(JState *state);

void update_sigma(JState *state);

void update_ag(JState *state);

void set_parameters( JState *state,  double f, double fp, double fptilde );

void init_state( JState *state, double f, double fp, double fptilde );

#endif // JSTATE_H
