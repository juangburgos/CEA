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
extern "C" void __declspec(dllexport) __cdecl sqltest(integer *m, integer *n, doublereal *mymat);
extern "C" void __declspec(dllexport) __cdecl matobserver(doublereal *u, doublereal *ym, doublereal *y);
#else
extern "C" void __declspec(dllimport) __cdecl mypinv(integer *m, integer *n, doublereal *a, doublereal *mytol, doublereal *ainv);
extern "C" void __declspec(dllimport) __cdecl sqltest(integer *m, integer *n, doublereal *mymat);
extern "C" void __declspec(dllimport) __cdecl matobserver(doublereal *u, doublereal *ym, doublereal *y);
#endif

void monpen(integer *m, integer *n, doublereal *a, doublereal *mytol, doublereal *ainvt);

// Matrix Structure
typedef struct mat{

    int m;
    int n;
    doublereal * mat;

}mat;

// Observer Functions
void zeros( int m, int n, mat *matrix );

// SQL Functions
void readmatsql( sqlite3 *db, char *matname, mat *matrix );

#endif /* MYDLL_H_ */
