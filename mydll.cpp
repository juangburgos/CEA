/*
 * mydll.cpp
 *
 *  Created on: Apr 26, 2013
 *      Author: USER
 */

#include <stdlib.h>
#include <float.h>
#include "stdafx.h"
#define EXPORTING_DLL

// to avoid C++ incompatibilities (after all, it is expecting C)
extern "C" {               // make sure these are visible in directory lookup
#include "f2c.h"           // our solution: copy them to /usr/include
#include "clapack.h"
}

#include "mydll.h"

extern "C" void __declspec(dllexport) __cdecl mypinv(integer *m, integer *n, doublereal *a, doublereal *mytol, doublereal *ainv)
{
	  // moore–penrose pseudoinverse
	  if (*n > *m)
	  {
		  // allocate at matrix
		  doublereal *at = (doublereal*) calloc((*n)*(*m), sizeof(doublereal));
		  // transpose of matrix a
		  for (int i = 0; i < *m; i++) //stored in column major
		  {
			  for(int j = 0; j < *n; j++)
			  {
				  at[(*n)*(i)+j] = a[j*(*m)+i];
			  }
		  }
		  // apply same algorithm to a^t and transpose the result
		  monpen(n,m,at,mytol,ainv);
		  // free at matrix
		  free (at);
	  }
	  else
	  {
		  // allocate ainvt matrix
		  doublereal *ainvt = (doublereal*) calloc((*m)*(*n), sizeof(doublereal));
		  // apply normal algorithm to a
		  monpen(m,n,a,mytol,ainvt);
		  // transpose of matrix ainvt
		  for (int i = 0; i < *m; i++) //stored in column major
		  {
			  for(int j = 0; j < *n; j++)
			  {
				  ainv[(*n)*(i)+j] = ainvt[j*(*m)+i];
			  }
		  }
		  // free ainvt matrix
		  free (ainvt);
	  }
}

void monpen(integer *m, integer *n, doublereal *a, doublereal *mytol, doublereal *ainvt)
{
	// init tolerance
    doublereal tol = 0.0;
    // define tolerance
//    if (*mytol <= 0)
//    {
//    	// linear search for machine tolerance (**eliminate this when possible**)
//		doublereal machEps = 1.0;
//		do {
//		   machEps /= 2.0;
//		}
//		while ((doublereal)(1.0 + (machEps/2.0)) != 1.0);
//		// tolerance
//		tol = (*m) * machEps;
//    }
//    else
//    {
    	tol = *mytol;
//    }
    // allocate work matrix
    integer lwork = 5*(*m)*(*n);
    doublereal *work = (doublereal*) calloc(lwork, sizeof(doublereal));
    // allocate s matrix
    doublereal *s = (doublereal*) calloc(*n, sizeof(doublereal)); // *m > *n for sure
    // allocate u matrix
    doublereal *u = (doublereal*) calloc((*m)*(*m), sizeof(doublereal));
    // allocate vt matrix
    doublereal *vt = (doublereal*) calloc((*n)*(*n), sizeof(doublereal));
    // dgesvd function output
    integer info  = 0.0;
    // singular value decomposition
    char jobu      = 'A';
    char jobvt     = 'A';
    dgesvd_(&jobu, &jobvt, m, n, a,
		    m, s, u, m, vt, n,
		    work, &lwork, &info);
    // allocate sfull matrix
    doublereal *sfull = (doublereal*) calloc((*m)*(*n), sizeof(doublereal));
    // fill sfull with zeros
    for (int i = 0; i < (*m)*(*n); i++)
    {
    	sfull[i] = 0.0;
    }
    // fill diagonal of sfull with inverse of singular values
    for (int i = 0; i < *n; i++)
    {
		if (s[i] > tol)
		{
			sfull[i+i*(*m)] = 1.0/(s[i]);
		}
    }
    // allocate usfull matrix
    doublereal *usfull = (doublereal*) calloc((*m)*(*n), sizeof(doublereal));
    // multiply u*sfull
    char transa      = 'N';
	char transb      = 'N';
	doublereal alpha = 1.0;
	doublereal beta  = 0.0;
	dgemm_(&transa, &transb,
	       m, n, m, &alpha, u, m,
		   sfull, m, &beta, usfull, m);
	// multiply usfull*vt
	dgemm_(&transa, &transb,
		   m, n, n, &alpha, usfull, m,
		   vt, n, &beta, ainvt, m);

	free (work);
    free (s);
    free (u);
    free (vt);
    free (sfull);
    free (usfull);
}
