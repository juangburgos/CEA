/*
 * mydll.cpp
 *
 *  Created on: Apr 26, 2013
 *      Author: USER
 */

#include <stdlib.h>
#include "stdafx.h"
#define EXPORTING_DLL

// to avoid C++ incompatibilities (after all, it is expecting C)
extern "C" {               // make sure these are visible in directory lookup
#include "f2c.h"           // our solution: copy them to /usr/include
#include "clapack.h"
}

#include "mydll.h"

extern "C" void __declspec(dllexport) __cdecl dgesvd(integer *m, integer *n, doublereal *a,
													 doublereal *s, doublereal *u, doublereal *vt)
{

	  char jobu      = 'A';
	  char jobvt     = 'A';
	  integer info   = 0.0;
	  integer minmn  = 0;

	  integer lwork = 5*(*m)*(*n);
	  doublereal *work = (doublereal*) calloc(lwork, sizeof(doublereal));

	  minmn = min(*m,*n);
	  doublereal *sigma = (doublereal*) calloc(minmn, sizeof(doublereal));

	  dgesvd_(&jobu, &jobvt, m, n, a,
			  m, sigma, u, m, vt, n,
			  work, &lwork, &info);

	  // fill diagonal of s with singular values
	  for (int i = 0; i < minmn; i++)
	  {
			s[i+i*(*m)] = sigma[i];
	  }

	  free (work);

}
