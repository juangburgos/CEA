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
extern "C" void __declspec(dllexport) __cdecl matcontroller(doublereal *u, doublereal *y, doublereal *ym, doublereal *wayPointX, doublereal *wayPointY, int numWaypoints, doublereal *radio, doublereal *seconds);
extern "C" void __declspec(dllexport) __cdecl Control(double *position, double *velocity, double *action, int numAxis, double *wayPointX, double *wayPointY, int numWaypoints, double *actualWayPoint, double *param, int numParam, double *seconds);
#else
extern "C" void __declspec(dllimport) __cdecl mypinv(integer *m, integer *n, doublereal *a, doublereal *mytol, doublereal *ainv);
extern "C" void __declspec(dllimport) __cdecl sqltest(integer *m, integer *n, doublereal *mymat);
extern "C" void __declspec(dllimport) __cdecl matobserver(doublereal *u, doublereal *ym, doublereal *y);
extern "C" void __declspec(dllimport) __cdecl matcontroller(doublereal *u, doublereal *y, doublereal *ym, doublereal *wayPointX, doublereal *wayPointY, int numWaypoints, doublereal *radio, doublereal *seconds);
extern "C" void __declspec(dllimport) __cdecl Control(double *position, double *velocity, double *action, int numAxis, double *wayPointX, double *wayPointY, int numWaypoints, double *actualWayPoint, double *param, int numParam, double *seconds);
#endif

void monpen(integer *m, integer *n, doublereal *a, doublereal *mytol, doublereal *ainvt);
void pinv(integer *m, integer *n, doublereal *a, doublereal *mytol, doublereal *ainv);

void monpennew(integer *m, integer *n, doublereal *a, doublereal *mytol, doublereal *ainvt,
		doublereal *work, doublereal *s, doublereal *u, doublereal *vt, doublereal *sfull, doublereal *usfull);
void pinvnew(integer *m, integer *n, doublereal *a, doublereal *mytol, doublereal *ainv,
		doublereal *work, doublereal *s, doublereal *u, doublereal *vt, doublereal *sfull, doublereal *usfull, doublereal *at, doublereal *ainvt);

// Matrix Structure
typedef struct mat{

    integer m;
    integer n;
    doublereal * mat;

}mat;

// Observer Functions
void zeros( mat *matrix );
void multmat( mat *A, mat *B, mat*C );

// Least Squares Function
//void myls( mat *R, mat *F, mat*x );

// SQL Functions
void readmatsql( sqlite3 *db, char *matname, mat *matrix );

// Fitting Function
doublereal sinfit( doublereal u, mat *fit_params );

// Third Order Setpoint Function (allocates memory for tx and xp !!! must be free before updating)
void thirdord( doublereal p, doublereal v, doublereal a, doublereal j, doublereal Ts, mat *tx, mat *xp );

// Linear Interpolation Function
void interpola ( mat *tr_x, mat *pr_x, mat *my_tr, mat *my_pr );

#endif /* MYDLL_H_ */
