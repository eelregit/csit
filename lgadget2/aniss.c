#include <stdio.h>
#include <gsl/gsl_spline.h>
#include "allvars.h"
#include "proto.h"

/* adding ANIsotroptic Super-Sample modes (functions of time) to the simulation background
 *
 * the input should be a text file tabulating (a, delta_a_x, delta_a_y, delta_a_z)
 * preceded by the number of entries, where delta_a = da/a
 *
 * change velocity conversion factor All.Time^1.5 -> All.TimeAni[]^1.5
 * mods are done in io.c and init.c, as well as in 2lptic
 *
 * in run.c following move_particles() (drift), update All.TimeAni[] right after All.Time
 * eval_aniss() in the middle of drift in predict.c
 *
 * TreePM forces modified in pm_periodic.c and forcetree.c, with kick unchanged
 *
 * in forcetree.c the splitting and softening scales are rendered
 * isotropic in physical coordinates (thus anisotropic in comoving sizes)
 */

static gsl_spline *splx, *sply, *splz;
static gsl_interp_accel *accx, *accy, *accz;

int init_aniss()
{
    FILE *fp;
    if(!(fp = fopen(All.AniFile, "r")))
    {
        printf("error in opening file '%s'\n", All.AniFile);
        endrun(870510);
    }
    int numa, i;
    fscanf(fp, "%d\n", &numa);
    double *a = mymalloc(numa * sizeof(double));
    double *dax = mymalloc(numa * sizeof(double));
    double *day = mymalloc(numa * sizeof(double));
    double *daz = mymalloc(numa * sizeof(double));
    for(i=0; i<numa; ++i)
    {
        fscanf(fp, "%lf %lf %lf %lf\n", a+i, dax+i, day+i, daz+i);
    }
    fclose(fp);

    splx = gsl_spline_alloc(gsl_interp_cspline, numa);
    sply = gsl_spline_alloc(gsl_interp_cspline, numa);
    splz = gsl_spline_alloc(gsl_interp_cspline, numa);
    accx = gsl_interp_accel_alloc();
    accy = gsl_interp_accel_alloc();
    accz = gsl_interp_accel_alloc();
    gsl_spline_init(splx, a, dax, numa);
    gsl_spline_init(sply, a, day, numa);
    gsl_spline_init(splz, a, daz, numa);

    return numa;
}

void eval_aniss(double Time, double* TimeAni)
{
    TimeAni[0] = Time * (1. + gsl_spline_eval(splx, Time, accx));
    TimeAni[1] = Time * (1. + gsl_spline_eval(sply, Time, accy));
    TimeAni[2] = Time * (1. + gsl_spline_eval(splz, Time, accz));
}
