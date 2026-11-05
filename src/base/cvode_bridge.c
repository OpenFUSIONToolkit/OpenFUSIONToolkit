/*-----------------------------------------------------------------------------
* Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
*
* SPDX-License-Identifier: LGPL-3.0-only
*------------------------------------------------------------------------------
* Minimal CVODE (SUNDIALS) interface for non-stiff ODE integration
*
* Wraps serial CVODE using the Adams-Moulton method with functional iteration
* (the equivalent of LSODE MF=10) behind a small opaque-handle API that can be
* called from Fortran. Supports SUNDIALS 6.0 and newer.
*----------------------------------------------------------------------------*/
#ifdef HAVE_CVODE
#include <stdlib.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#if SUNDIALS_VER_MAJOR >= 7
#define OFT_SUN_COMM_NULL SUN_COMM_NULL
#else
#define OFT_SUN_COMM_NULL NULL
#endif

/* Right-hand-side callback implemented in Fortran: ydot = f(t,y) */
typedef int (*oft_cvode_rhs_c)(int neq, double t, const double *y, double *ydot, void *udata);

typedef struct {
  void *mem;                 /* CVODE memory */
  N_Vector y;                /* Solution vector */
  SUNNonlinearSolver nls;    /* Functional iteration solver */
  SUNContext ctx;            /* SUNDIALS context (one per integrator) */
  int neq;                   /* Number of equations */
  oft_cvode_rhs_c f;         /* Fortran callback */
  void *udata;               /* Fortran callback data */
} oft_cvode_t;

/* Adapter between CVODE callback signature and Fortran callback */
static int oft_cvode_rhs_adapter(double t, N_Vector y, N_Vector ydot, void *user_data)
{
  oft_cvode_t *self = (oft_cvode_t *) user_data;
  return self->f(self->neq, t, N_VGetArrayPointer(y), N_VGetArrayPointer(ydot), self->udata);
}

static void oft_cvode_destroy(oft_cvode_t *self)
{
  if (self->mem != NULL) CVodeFree(&self->mem);
  if (self->nls != NULL) SUNNonlinSolFree(self->nls);
  if (self->y != NULL) N_VDestroy(self->y);
  if (self->ctx != NULL) SUNContext_Free(&self->ctx);
  free(self);
}

/* Create integrator and set initial condition: returns 0 on success */
int oft_cvode_create_c(int neq, double t0, const double *y0, double rtol, double atol,
                       oft_cvode_rhs_c f, void *udata, oft_cvode_t **handle)
{
  int i, flag;
  double *yd;
  oft_cvode_t *self = (oft_cvode_t *) calloc(1, sizeof(oft_cvode_t));
  *handle = NULL;
  if (self == NULL) return -1;
  self->neq = neq;
  self->f = f;
  self->udata = udata;
  if (SUNContext_Create(OFT_SUN_COMM_NULL, &self->ctx) != 0) { free(self); return -2; }
  self->y = N_VNew_Serial((sunindextype) neq, self->ctx);
  if (self->y == NULL) { oft_cvode_destroy(self); return -3; }
  yd = N_VGetArrayPointer(self->y);
  for (i = 0; i < neq; i++) yd[i] = y0[i];
  self->mem = CVodeCreate(CV_ADAMS, self->ctx);
  if (self->mem == NULL) { oft_cvode_destroy(self); return -4; }
  flag = CVodeInit(self->mem, oft_cvode_rhs_adapter, t0, self->y);
  if (flag != CV_SUCCESS) { oft_cvode_destroy(self); return -5; }
  flag = CVodeSetUserData(self->mem, self);
  if (flag != CV_SUCCESS) { oft_cvode_destroy(self); return -6; }
  flag = CVodeSStolerances(self->mem, rtol, atol);
  if (flag != CV_SUCCESS) { oft_cvode_destroy(self); return -7; }
  /* Use functional iteration (non-stiff) so no linear solver is required */
  self->nls = SUNNonlinSol_FixedPoint(self->y, 0, self->ctx);
  if (self->nls == NULL) { oft_cvode_destroy(self); return -8; }
  flag = CVodeSetNonlinearSolver(self->mem, self->nls);
  if (flag != CV_SUCCESS) { oft_cvode_destroy(self); return -9; }
  *handle = self;
  return 0;
}

/* Restart integration from (t0,y0), optionally with initial step size h0 (<=0 for default) */
int oft_cvode_reinit_c(oft_cvode_t *self, double t0, const double *y0, double h0)
{
  int i, flag;
  double *yd = N_VGetArrayPointer(self->y);
  for (i = 0; i < self->neq; i++) yd[i] = y0[i];
  flag = CVodeReInit(self->mem, t0, self->y);
  if (flag != CV_SUCCESS) return flag;
  return CVodeSetInitStep(self->mem, (h0 > 0.0) ? h0 : 0.0);
}

/* Take one internal step: returns CVODE flag (0=success, 1=reached tstop, <0=failure) */
int oft_cvode_step_c(oft_cvode_t *self, double tout, int use_tstop, double tstop,
                     double *y, double *t, double *hlast)
{
  int i, flag;
  double *yd = N_VGetArrayPointer(self->y);
  if (use_tstop) {
    flag = CVodeSetStopTime(self->mem, tstop);
    if (flag != CV_SUCCESS) return flag;
  }
  flag = CVode(self->mem, tout, self->y, t, CV_ONE_STEP);
  for (i = 0; i < self->neq; i++) y[i] = yd[i];
  if (flag >= 0) CVodeGetLastStep(self->mem, hlast);
  return flag;
}

/* Destroy integrator and release all memory */
void oft_cvode_free_c(oft_cvode_t **handle)
{
  if (*handle != NULL) oft_cvode_destroy(*handle);
  *handle = NULL;
}
#endif
