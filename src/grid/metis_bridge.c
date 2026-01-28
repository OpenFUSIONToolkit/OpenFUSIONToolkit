/*-----------------------------------------------------------------------------
* Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
*
* SPDX-License-Identifier: LGPL-3.0-only
*------------------------------------------------------------------------------
* METIS interface for the Open FUSION Toolkit
*----------------------------------------------------------------------------*/
#include <stdlib.h>
#include "metis.h"

void
oft_metis_partMesh(int *nc, int *np, int *ncp, int *lc, int *npart, int *cpart, int *info)
{
	idx_t *vwgt;
	idx_t *vsize;
	idx_t *eptr;
	idx_t *eind;
	idx_t *ppart;
	real_t *tpwgts;
	idx_t options[METIS_NOPTIONS];
	idx_t objval;
	idx_t ncommon=3;
	if (*ncp == 8) { ncommon=4; }
	// Allocate temporary partitioning variables
	vwgt = (idx_t *) malloc(*nc * sizeof(idx_t));
	vsize = (idx_t *) malloc(*nc * sizeof(idx_t));
	eptr = (idx_t *) malloc((*nc +1) * sizeof(idx_t));
	eind = (idx_t *) malloc((*ncp)*(*nc) * sizeof(idx_t));
	ppart = (idx_t *) malloc(*np * sizeof(idx_t));
	tpwgts = (real_t *) malloc(*npart * sizeof(real_t));
	// Initialize constraints and sizes for partitioning
	real_t tmp_wgt = 1.0f / *npart;
	int i = 0;
	for (i = 0; i < *npart; i++) { tpwgts[i] = tmp_wgt; }
	for (i = 0; i < *nc; i++) {
		vwgt[i] = 1;
		vsize[i] = 1;
		eptr[i] = (*ncp)*i;
	}
	eptr[*nc] = (*ncp)*(*nc);
	for (i = 0; i < (*ncp)*(*nc); i++) { eind[i] = lc[i]-1; }
	// Perform partitioning with default options
	*info = METIS_SetDefaultOptions(options);
	*info = METIS_PartMeshDual(nc,np,eptr,eind,vwgt,vsize,&ncommon,
			                   npart,tpwgts,options,&objval,cpart,ppart);
	// Free temporary variables
	free(vwgt);
	free(vsize);
	free(eptr);
	free(eind);
	free(ppart);
	free(tpwgts);
	// Update partitioning to 1-based indexing
	for (i = 0; i < *nc; i++) { cpart[i]++; }
	// Finished
	return;
}

void
oft_metis_partGraph(int *nr, int *nnz, int *kr, int *lc, int *npart, int *part, int *type, int *info)
{
	idx_t ncon = 1;
	idx_t *vwgt;
	idx_t *vsize;
	idx_t *adjwgt;
	real_t *tpwgts;
	real_t *ubvec;
	idx_t options[METIS_NOPTIONS];
	idx_t objval;
	// Allocate temporary partitioning variables
	vwgt = (idx_t *) malloc(*nr * sizeof(idx_t));
	vsize = (idx_t *) malloc(*nr * sizeof(idx_t));
	adjwgt = (idx_t *) malloc(*nnz * sizeof(idx_t));
	tpwgts = (real_t *) malloc(*npart * sizeof(real_t));
	ubvec = (real_t *) malloc(ncon * sizeof(real_t));
	// Initialize constraints and sizes for partitioning
	real_t tmp_wgt = 1.0f / *npart;
	int i = 0;
	for (i = 0; i < *npart; i++) { tpwgts[i] = tmp_wgt; }
	for (i = 0; i < ncon; i++) { ubvec[i] = (real_t) 1.01; }
	for (i = 0; i < *nr; i++) {
		vwgt[i] = 1;
		vsize[i] = 1;
	}
	for (i = 0; i < *nnz; i++) { adjwgt[i] = 1; }
	// Update graph to 0-based indexing
	for (i = 0; i <= *nr; i++) { kr[i]--; }
	for (i = 0; i < *nnz; i++) { lc[i]--; }
	// Perform partitioning with default options
	*info = METIS_SetDefaultOptions(options);
	if (*type == 1) {
		*info = METIS_PartGraphRecursive(nr,&ncon,kr,lc,vwgt,vsize,adjwgt,
				                         npart,tpwgts,ubvec,options,&objval,part);
	} else if (*type == 2) {
		*info = METIS_PartGraphKway(nr,&ncon,kr,lc,vwgt,vsize,adjwgt,
				                    npart,tpwgts,ubvec,options,&objval,part);
	} else {
		*info = -99;
	}
	// Free temporary variables
	free(vwgt);
	free(vsize);
	free(adjwgt);
	free(tpwgts);
	free(ubvec);
	// Update partitioning to 1-based indexing
	for (i = 0; i < *nr; i++) { part[i]++; }
	// Reset graph to 1-based indexing
	for (i = 0; i <= *nr; i++) { kr[i]++; }
	for (i = 0; i < *nnz; i++) { lc[i]++; }
	// Finished
	return;
}
