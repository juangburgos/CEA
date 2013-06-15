/*
 * mydll.h
 *
 *  Created on: Apr 26, 2013
 *      Author: USER
 */

#ifndef MYDLL_H_
#define MYDLL_H_

#ifdef EXPORTING_DLL
extern "C" void __declspec(dllexport) __cdecl dgesvd(integer *m, integer *n, doublereal *a, doublereal *s, doublereal *u, doublereal *vt);
#else
extern "C" void __declspec(dllimport) __cdecl dgesvd(integer *m, integer *n, doublereal *a, doublereal *s, doublereal *u, doublereal *vt);
#endif

#endif /* MYDLL_H_ */
