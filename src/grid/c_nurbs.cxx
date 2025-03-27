/*-----------------------------------------------------------------------------
* Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
*
* SPDX-License-Identifier: LGPL-3.0-only
*------------------------------------------------------------------------------
* OpenNURBS interface for the Open FUSION Toolkit
*----------------------------------------------------------------------------*/
#ifdef HAVE_ONURBS
#include "opennurbs.h"

// Setup Logging
ON_TextLog dump_to_stderr(stderr); //!< Log variable
ON_TextLog* dump = &dump_to_stderr; //!< Log variable
bool bVerboseTextDump = true;

// Model container
ONX_Model model;

// CAD wireframe objects
int ncurves = -1;
int nsurfs = -1;
ON_ClassArray<ONX_Model_Object> curves;
ON_ClassArray<ONX_Model_Object> surfs;

int nurbs_curve_name( int ind, char * name)
{
	// Validate index
	if ( ind < 0 || ind > ncurves-1 ) { return -1; }
	//
	int len_s;
	len_s = strlen(name);
	//
	const wchar_t* tmp_name = curves[ind].m_attributes.m_name;
	int nlength = curves[ind].m_attributes.m_name.Length();
	//printf("Curve %i %i, %S\n",ind,len_s,tmp_name);
	if ( len_s < nlength ) { return -2; }
	for (int i=0; i<nlength; i++) {
		name[i]=tmp_name[i];
	}
	//
	return 0;
}

int nurbs_eval_curve( int ind, double u, double *r )
{
	// Validate index
	if ( ind < 0 || ind > ncurves-1 ) { return -1; }
	// Get curve reference
	const ON_Object* m_object = curves[ind].m_object;
	// Cast to curve object
	const ON_Curve* pCurve = ON_Curve::Cast(m_object);
	// Create point and evaluate
	ON_3dPoint pt;
	ON_BOOL32 chk = pCurve->EvPoint(u,pt);
	if ( !chk ) { return -2; }
	// Copy point to vector
	r[0] = pt.x;
	r[1] = pt.y;
	r[2] = pt.z;
	// Return
	return 0;
}

int nurbs_curve_span_size( int ind, int *size )
{
	// Validate index
	if ( ind < 0 || ind > ncurves-1 ) { return -1; }
	// Get curve reference
	const ON_Object* m_object = curves[ind].m_object;
	// Cast to curve object
	const ON_Curve* pCurve = ON_Curve::Cast(m_object);
	// Get parameter span size
	*size = pCurve->SpanCount();
	// Return
	return 0;
}

int nurbs_curve_span( int ind, double *span )
{
	// Validate index
	if ( ind < 0 || ind > ncurves-1 ) { return -1; }
	// Get curve reference
	const ON_Object* m_object = curves[ind].m_object;
	// Cast to curve object
	const ON_Curve* pCurve = ON_Curve::Cast(m_object);
	// Get parameter span size
	ON_BOOL32 chk = pCurve->GetSpanVector(span);
	if ( !chk ) { return -2; }
	// Return
	return 0;
}

int nurbs_curve_domain( int ind, double *domain )
{
	// Validate index
	if ( ind < 0 || ind > ncurves-1 ) { return -1; }
	// Get curve reference
	const ON_Object* m_object = curves[ind].m_object;
	// Cast to curve object
	const ON_Curve* pCurve = ON_Curve::Cast(m_object);
	// Get parameter span size
	double t1,t2;
	ON_BOOL32 chk;
	chk = pCurve->GetDomain(&t1,&t2);
	if ( !chk ) { return -2; }
	domain[0] = t1;
	domain[1] = t2;
	// Return
	return 0;
}

int nurbs_curve_periodic( int ind, int *p)
{
	// Validate index
	if ( ind < 0 || ind > ncurves-1 ) { return -1; }
	// Get curve reference
	const ON_Object* m_object = curves[ind].m_object;
	// Cast to curve object
	const ON_Curve* pCurve = ON_Curve::Cast(m_object);
	// Get parameter span size
	ON_BOOL32 chk;
	//
	*p = 0;
	chk = pCurve->IsPeriodic();
	if ( chk ) { *p = 1; }
	//
	return 0;
}

int nurbs_curve_linear(int ind, int *p)
{
  *p=0;
  // Validate index
  if ( ind < 0 || ind > ncurves-1 ) { return -1; }
  // Get curve reference
  const ON_Object* m_object = curves[ind].m_object;
  // Cast to curve object
  const ON_Curve* pCurve = ON_Curve::Cast(m_object);
  // Get parameter span size
  ON_BOOL32 chk;
  // Check if planar
  chk = pCurve->IsLinear(ON_ZERO_TOLERANCE);
  // Set return
  if (chk) { *p = 1; }
  //
  return 0;
}

int nurbs_surf_name( int ind, char * name)
{
	// Validate index
	if ( ind < 0 || ind > nsurfs-1 ) { return -1; }
	//
	int len_s;
	len_s = strlen(name);
	//
	const wchar_t* tmp_name = surfs[ind].m_attributes.m_name;
	int nlength = surfs[ind].m_attributes.m_name.Length();
	//
	if ( len_s < nlength ) { return -2; }
	for (int i=0; i<nlength; i++) {
		name[i]=tmp_name[i];
	}
	//
	return 0;
}

int nurbs_eval_surf( int ind, double u, double v, double *r )
{
	// Validate index
	if ( ind < 0 || ind > nsurfs-1 ) { return -1; }
	// Get surface reference
	const ON_Object* m_object = surfs[ind].m_object;
	// Cast to surface object
	const ON_Surface* pSurface = ON_Surface::Cast(m_object);
	// Create point and evaluate
	ON_3dPoint pt;
	ON_BOOL32 chk = pSurface->EvPoint(u,v,pt);
	if ( !chk ) { return -2; }
	// Copy point to vector
	r[0] = pt.x;
	r[1] = pt.y;
	r[2] = pt.z;
	// Return
	return 0;
}

int nurbs_surf_span_size( int ind, int *size1, int *size2 )
{
	// Validate index
	if ( ind < 0 || ind > nsurfs-1 ) { return -1; }
	// Get curve reference
	const ON_Object* m_object = surfs[ind].m_object;
	// Cast to curve object
	const ON_Surface* pSurface = ON_Surface::Cast(m_object);
	// Get parameter span size
	*size1 = pSurface->SpanCount(0);
	*size2 = pSurface->SpanCount(1);
	// Return
	return 0;
}

int nurbs_surf_span( int ind, double *span1, double *span2 )
{
	// Validate index
	if ( ind < 0 || ind > nsurfs-1 ) { return -1; }
	// Get curve reference
	const ON_Object* m_object = surfs[ind].m_object;
	// Cast to curve object
	const ON_Surface* pSurface = ON_Surface::Cast(m_object);
	// Get parameter span size
	double t1,t2;
	ON_BOOL32 chk;
	chk = pSurface->GetSpanVector(0,span1);
	if ( !chk ) { return -2; }
	chk = pSurface->GetSpanVector(1,span2);
	if ( !chk ) { return -2; }
	chk = pSurface->GetDomain(0,&t1,&t2);
	//printf("Surface domain %14.10f, %14.10f",t1,t2);
	chk = pSurface->GetDomain(1,&t1,&t2);
	//printf("Surface domain %14.10f, %14.10f",t1,t2);
	// Return
	return 0;
}

int nurbs_surf_domain( int ind, double *domain1, double *domain2 )
{
	// Validate index
	if ( ind < 0 || ind > nsurfs-1 ) { return -1; }
	// Get curve reference
	const ON_Object* m_object = surfs[ind].m_object;
	// Cast to curve object
	const ON_Surface* pSurface = ON_Surface::Cast(m_object);
	// Get parameter span size
	double t1,t2;
	ON_BOOL32 chk;
	chk = pSurface->GetDomain(0,&t1,&t2);
	if ( !chk ) { return -2; }
	domain1[0] = t1;
	domain1[1] = t2;
	//
	chk = pSurface->GetDomain(1,&t1,&t2);
	if ( !chk ) { return -2; }
	domain2[0] = t1;
	domain2[1] = t2;
	// Return
	return 0;
}

int nurbs_surf_periodic( int ind, int *p1, int *p2)
{
	// Validate index
	if ( ind < 0 || ind > nsurfs-1 ) { return -1; }
	// Get curve reference
	const ON_Object* m_object = surfs[ind].m_object;
	// Cast to curve object
	const ON_Surface* pSurface = ON_Surface::Cast(m_object);
	// Get parameter span size
	ON_BOOL32 chk;
	//
	*p1 = 0;
	chk = pSurface->IsPeriodic(0);
	if ( chk ) { *p1 = 1; }
	//
	*p2 = 0;
	chk = pSurface->IsPeriodic(1);
	if ( chk ) { *p2 = 1; }
	//
	return 0;
}

int nurbs_surf_singular( int ind, int *p)
{
	// Validate index
	if ( ind < 0 || ind > nsurfs-1 ) { return -1; }
	// Get curve reference
	const ON_Object* m_object = surfs[ind].m_object;
	// Cast to curve object
	const ON_Surface* pSurface = ON_Surface::Cast(m_object);
	// Get parameter span size
	ON_BOOL32 chk;
	//
	for ( int i=0; i<4; i++ )
	{
		p[i] = 0;
		chk = pSurface->IsSingular(i);
		if ( chk ) { p[i] = 1; }
	}
	//
	return 0;
}

int nurbs_surf_planar(int ind, int *p)
{
  *p = 0;
  // Validate index
  if ( ind < 0 || ind > nsurfs-1 ) { return -1; }
  // Get curve reference
  const ON_Object* m_object = surfs[ind].m_object;
  // Cast to curve object
  const ON_Surface* pSurface = ON_Surface::Cast(m_object);
  // Get parameter span size
  ON_BOOL32 chk;
  // Check if planar
  chk = pSurface->IsPlanar(NULL,ON_ZERO_TOLERANCE);
  // Set return
  if (chk) { *p = 1; }
  //
  return 0;
}

int nurbs_read_in( char *filename )
{
	// open file containing opennurbs archive
	FILE* archive_fp = ON::OpenFile( filename, "rb");
	// Error on open
	if ( !archive_fp ) { return -1; }
	// create achive object from file pointer
	ON_BinaryFile archive( ON::read3dm, archive_fp );
	// read the contents of the file into "model"
	bool rc = model.Read( archive, dump );
	// close the file
	ON::CloseFile( archive_fp );
	// Error on read
	if ( !rc ) { return -2; }
	// Error in model
	if ( !model.IsValid(dump) ) { return -3; }
	//
	const ON_Object* m_object;
	const ON_Geometry* pGeometry;
	//
	ncurves = 0;
	nsurfs = 0;
	for ( int i = 0; i < model.m_object_table.Count(); i++ )
	{
		// Object definition
		m_object = model.m_object_table[i].m_object;
		pGeometry = ON_Geometry::Cast(m_object);
		if ( pGeometry )
		{
			// m_object is some type of geometric object
			if ( ON_Extrusion::Cast(m_object) )
			{
				// m_object is derived from ON_Extrusion
				//  Note:
				//   ON_Extrusion::BrepForm() will return a brep form
				//   if you don't want to explicitly handle extrusions.
				//const ON_Extrusion* extrusion = ON_Extrusion::Cast(m_object);
				dump->Print("Extrusion found.\n");
				return -4;
			}
			else if ( ON_Brep::Cast(m_object) )
			{
				// m_object is derived from ON_Brep
				//const ON_Brep* brep = ON_Brep::Cast(m_object);
				dump->Print("Brep found.\n");
				return -4;
			}
			else if ( ON_Curve::Cast(m_object) )
			{
				//const ON_Curve* pCurve = ON_Curve::Cast(m_object);
				ncurves++;
				ONX_Model_Object& mo = curves.AppendNew();
				mo.m_object = m_object;
				mo.m_attributes = model.m_object_table[i].m_attributes;
			}
			else if ( ON_Surface::Cast(m_object) )
			{
				//const ON_Surface* pSurface = ON_Surface::Cast(m_object);
				nsurfs++;
				ONX_Model_Object& mo = surfs.AppendNew();
				mo.m_object = m_object;
				mo.m_attributes = model.m_object_table[i].m_attributes;
			}
			else if ( ON_Mesh::Cast(m_object) )
			{
				//const ON_Mesh* pMesh = ON_Mesh::Cast(m_object);
				dump->Print("Mesh found.\n");
				return -4;
			}
			else
			{
				dump->Print("Invalid object.\n");
				return -4;
			}
		}
	}
	// Return
	return 0;
}

int nurbs_init()
{
	// Initialize OpenNURBS library
	ON::Begin();
	// Return
	return 0;
}

int nurbs_finalize()
{
	// destroy this model
	model.Destroy();
	// Finalize OpenNURBS library
	ON::End();
	// Return
	return 0;
}

extern "C" {
	//
	void nurbs_init(int *ierr) { *ierr = nurbs_init(); }
	void nurbs_finalize(int *ierr) { *ierr = nurbs_finalize(); }
	void nurbs_read_in(char *filename, int *ierr) { *ierr = nurbs_read_in( filename ); }
	void nurbs_get_count(int *nc, int *ns, int *ierr) {
		if ( ncurves < 0 || nsurfs < 0 ) {
			*ierr = -1;
		} else {
			*nc = ncurves;
			*ns = nsurfs;
			*ierr = 0;
		}
	}
	//
	void nurbs_curve_name(int *ind, char *name, int *ierr) { *ierr = nurbs_curve_name( *ind, name ); }
	void nurbs_curve_span_size(int *ind, int *size, int *ierr) { *ierr = nurbs_curve_span_size( *ind, size ); }
	void nurbs_curve_span(int *ind, double *span, int *ierr) { *ierr = nurbs_curve_span( *ind, span ); }
	void nurbs_curve_domain(int *ind, double *domain, int *ierr) { *ierr = nurbs_curve_domain( *ind, domain ); }
	void nurbs_curve_periodic(int *ind, int *p, int *ierr) { *ierr = nurbs_curve_periodic( *ind, p); }
	void nurbs_curve_linear(int *ind, int *p, int *ierr) { *ierr = nurbs_curve_linear( *ind, p); }
	//
	void nurbs_surf_name(int *ind, char *name, int *ierr) { *ierr = nurbs_surf_name( *ind, name ); }
	void nurbs_surf_span_size(int *ind, int *size1, int *size2, int *ierr) { *ierr = nurbs_surf_span_size( *ind, size1, size2 ); }
	void nurbs_surf_span(int *ind, double *span1, double *span2, int *ierr) { *ierr = nurbs_surf_span( *ind, span1, span2 ); }
	void nurbs_surf_domain(int *ind, double *domain1, double *domain2, int *ierr) { *ierr = nurbs_surf_domain( *ind, domain1, domain2 ); }
	void nurbs_surf_periodic(int *ind, int *p1, int *p2, int *ierr) { *ierr = nurbs_surf_periodic( *ind, p1, p2); }
	void nurbs_surf_singular(int *ind, int *p, int *ierr) { *ierr = nurbs_surf_singular( *ind, p); }
	void nurbs_surf_planar(int *ind, int *p, int *ierr) { *ierr = nurbs_surf_planar( *ind, p); }
	//
	void nurbs_eval_curve(int *ind, double *u, double *r, int *ierr) { *ierr = nurbs_eval_curve( *ind, *u, r ); }
	void nurbs_eval_surf(int *ind, double *u, double *v, double *r, int *ierr) { *ierr = nurbs_eval_surf( *ind, *u, *v, r); }
}
#endif
