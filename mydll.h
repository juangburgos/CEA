/*
 * mydll.h
 *
 *  Created on: Apr 26, 2013
 *      Author: USER
 */

#ifndef MYDLL_H_
#define MYDLL_H_

#ifdef EXPORTING_DLL
extern "C" void __declspec(dllexport) __cdecl mypinv(integer *m, integer *n, doublereal *a, doublereal *mytol, doublereal *ainv);
#else
extern "C" void __declspec(dllimport) __cdecl mypinv(integer *m, integer *n, doublereal *a, doublereal *mytol, doublereal *ainv);
#endif

void monpen(integer *m, integer *n, doublereal *a, doublereal *mytol, doublereal *ainvt);

#endif /* MYDLL_H_ */
