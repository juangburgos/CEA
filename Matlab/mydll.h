/*
 * mydll.h
 *
 *  Created on: Apr 26, 2013
 *      Author: USER
 */

#ifndef MYDLL_H_
#define MYDLL_H_

#ifdef EXPORTING_DLL
extern "C" void __declspec(dllexport) __cdecl matobserver(double *u, double *ym, double *y);
extern "C" void __declspec(dllexport) __cdecl Control(double *position, double *velocity, double *action, int numAxis, double *wayPointX, double *wayPointY, int numWaypoints, double *actualWayPoint, double *param, int numParam);
#else
extern "C" void __declspec(dllimport) __cdecl matobserver(double *u, double *ym, double *y);
extern "C" void __declspec(dllimport) __cdecl Control(double *position, double *velocity, double *action, int numAxis, double *wayPointX, double *wayPointY, int numWaypoints, double *actualWayPoint, double *param, int numParam);
#endif

#endif /* MYDLL_H_ */
