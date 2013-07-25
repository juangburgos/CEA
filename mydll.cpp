/*
 * mydll.cpp
 *
 *  Created on: Apr 26, 2013
 *      Author: USER
 */
#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define EXPORTING_DLL

// to avoid C++ incompatibilities
extern "C" {
#include "./lapack/f2c.h"
#include "./lapack/clapack.h"
#include "./sqlite/sqlite3.h"
}

// Matrix Structure
typedef struct mat{

    integer m;
    integer n;
    doublereal * mat;

}mat;

// Matrix Functions
void zeros( mat *matrix );
void multmat( mat *A, mat *B, mat*C );
// Least Squares Function
void myls( mat *R, mat *F);
// SQL Functions
void readmatsql( sqlite3 *db, char *matname, mat *matrix );
// Fitting Function
doublereal sinfit( doublereal u, mat *fit_params );
// Third Order Setpoint Function (allocates memory for tx and xp)
void thirdord( doublereal p, doublereal v, doublereal a, doublereal j, doublereal Ts, mat *tx, mat *xp );
// Linear Interpolation Function
void interpola ( mat *tr_x, mat *pr_x, mat *my_tr, mat *my_pr );

#include "mydll.h"

doublereal Ts = 0.06;
int N  = 50;
doublereal radio = 0.05;

extern "C" void __declspec(dllexport) __cdecl matobserver(doublereal *u, doublereal *ym, doublereal *y)
{
	// % Controller Status Variable (Init = 0, Run = 1)
	static int cStatus = 0;
	// % Internal System Model State
	static mat x;
	// % Internal System Model A Matrix
	static mat sysA;
	// % Internal System Model B Matrix
	static mat sysB;
	// % Internal System Model C Matrix
	static mat sysC;
	// % Internal Disturbance Model State
	static mat xd;
	//% Internal Disturbance Model A Matrix
	static mat disA;
	// % Internal Disturbance Model B Matrix
	static mat disB;
	// % Internal Disturbance Model C Matrix
	static mat disC;
	// % Internal Disturbance Model D Matrix
	static mat disD;
	// % Roll Scheduling Parameters
	static mat sRoll;
	// % Pitch Scheduling Parameters
	static mat sPitch;

	// Input Input in Matrix Form
	static mat m_u;
	// Measured Output in Matrix Form
	static mat m_ym;
	// Estimated Output in Matrix Form
	static mat m_y;

	// Temporal Matrices
	static mat sysCinv;
	static mat Sched;
	static mat Bsch;
	static mat stEq1;
	static mat stEq2;
	static mat stEq3;
	static mat oErr;
	static mat dtEq1;
	static mat suEq1;
	static mat duEq1;

	// % If Observer Initialization
	if (cStatus == 0)
	{
		// % Initlialize Parameters and States
		sqlite3_initialize();
		sqlite3 *db;
		// Open Database
		int error = sqlite3_open_v2( "test.db", &db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL );
		if (error)
		{
			sqlite3_close(db);
			system("pause");
		}
		readmatsql(db, "Ap", &sysA);
		readmatsql(db, "Bp", &sysB);
		readmatsql(db, "Cp", &sysC);
		readmatsql(db, "Aj", &disA);
		readmatsql(db, "Bj", &disB);
		readmatsql(db, "Cj", &disC);
		readmatsql(db, "Dj", &disD);
		readmatsql(db, "rollparams", &sRoll);
		readmatsql(db, "pitchparams", &sPitch);
		// Close Database
		sqlite3_close(db);
		// Init Inputs, States and Outputs in Matrix Form
		m_u.m    = sysB.n;
		m_u.n    = 1;
		//m_u.mat  = (doublereal*) calloc(m_u.m*m_u.n, sizeof(doublereal));
		m_ym.m   = sysC.m;
		m_ym.n   = 1;
		//m_ym.mat = (doublereal*) calloc(m_ym.m*m_ym.n, sizeof(doublereal));
		m_y.m    = sysC.m;
		m_y.n    = 1;
		//m_y.mat  = (doublereal*) calloc(m_y.m*m_y.n, sizeof(doublereal));
		x.m      = sysA.m;
		x.n      = 1;
		x.mat    = (doublereal*) calloc(x.m*x.n, sizeof(doublereal));
		xd.m     = disA.m;
		xd.n     = 1;
		xd.mat   = (doublereal*) calloc(xd.m*xd.n, sizeof(doublereal));
		// Init Temporal Matrices
		sysCinv.m   = sysC.n;
		sysCinv.n   = sysC.m;
		sysCinv.mat = (doublereal*) calloc(sysCinv.m*sysCinv.n, sizeof(doublereal));
		Bsch.m      = sysB.m;
		Bsch.n      = sysB.n;
		Bsch.mat    = (doublereal*) calloc(Bsch.m*Bsch.n, sizeof(doublereal));
		Sched.m     = sysB.n;
		Sched.n     = sysB.n;
		Sched.mat   = (doublereal*) calloc(Sched.m*Sched.n, sizeof(doublereal));
		stEq1.m     = sysA.m;
		stEq1.n     = x.n;
		stEq1.mat   = (doublereal*) calloc(stEq1.m*stEq1.n, sizeof(doublereal));
		stEq2.m     = disC.m;
		stEq2.n     = xd.n;
		stEq2.mat   = (doublereal*) calloc(stEq2.m*stEq2.n, sizeof(doublereal));
		stEq3.m     = Bsch.m;
		stEq3.n     = m_u.n;
		stEq3.mat   = (doublereal*) calloc(stEq3.m*stEq3.n, sizeof(doublereal));
		oErr.m      = m_y.m;
		oErr.n      = m_y.n;
		oErr.mat    = (doublereal*) calloc(oErr.m*oErr.n, sizeof(doublereal));
		dtEq1.m     = disA.m;
		dtEq1.n     = xd.n;
		dtEq1.mat   = (doublereal*) calloc(dtEq1.m*dtEq1.n, sizeof(doublereal));
		suEq1.m     = disD.m;
		suEq1.n     = oErr.n;
		suEq1.mat   = (doublereal*) calloc(suEq1.m*suEq1.n, sizeof(doublereal));
		duEq1.m     = disB.m;
		duEq1.n     = oErr.n;
		duEq1.mat   = (doublereal*) calloc(duEq1.m*duEq1.n, sizeof(doublereal));
		// 0) Input/Output Pointers
		m_u.mat  = u;
		m_ym.mat = ym;
		m_y.mat  = y;
		// % Initial Conditions
		// Initial Plant States
		//   Calculate C inverse
//		pinv(&(sysC.m), &(sysC.n), sysC.mat, &mytol, sysCinv.mat);
//		//   Multiply x = Cinv * ym
//		multmat(&sysCinv, &m_ym, &x);
//		//   Init xd = zeros
//		zeros(&xd);
		// % 1) PREDICT
		// % Scheduled B matrix
		zeros(&Sched);
		Sched.mat[0] = sinfit(m_u.mat[0], &sRoll);
		Sched.mat[3] = sinfit(m_u.mat[1], &sPitch);
		multmat(&sysB, &Sched, &Bsch);
		// % Nominal Model State Equation
		multmat(&sysA, &x, &stEq1);
		multmat(&disC, &xd, &stEq2);
		multmat(&Bsch, &m_u, &stEq3);
		for(int i = 0; i < x.m; i++)
		{
			x.mat[i] = stEq1.mat[i] + stEq2.mat[i] + stEq3.mat[i];
		}
		// % Disturbance Model State Equation
		multmat(&disA, &xd, &dtEq1);
		for(int i = 0; i < xd.m; i++)
		{
			xd.mat[i] = dtEq1.mat[i];
		}
//		//fprintf (pFile, "Step 1 successfully %d \n",5);
//		//fflush (pFile);
		// % Output Equation
		multmat(&sysC, &x, &m_y);
		// % 2) CORRECT
		// Observer Error
		for(int i = 0; i < oErr.m; i++)
		{
			oErr.mat[i] = m_ym.mat[i] - m_y.mat[i];
		}
		// % Update Nominal Model States
		multmat(&disD, &oErr, &suEq1);
		for(int i = 0; i < x.m; i++)
		{
			x.mat[i] = x.mat[i] + suEq1.mat[i];
		}
		// % Update Disturbance Model States
		multmat(&disB, &oErr, &duEq1);
		for(int i = 0; i < xd.m; i++)
		{
			xd.mat[i] = xd.mat[i] + duEq1.mat[i];
		}
		// % Update Output
		multmat(&sysC, &x, &m_y);
//		//fprintf (pFile, "Step 2 successfully %d \n",6);
//		//fflush (pFile);
		// % Initialization Finished
		cStatus = 1;

	}
	else
	{
//		myCount++;
//		//fprintf (pFile, "RUN %d \n",0);
//		//fflush (pFile);
		// 0) Input/Output Pointers
		m_u.mat  = u;
		m_ym.mat = ym;
		m_y.mat  = y;
		// % 1) PREDICT
		// % Scheduled B matrix
		zeros(&Sched);
		Sched.mat[0] = sinfit(m_u.mat[0], &sRoll);
		Sched.mat[3] = sinfit(m_u.mat[1], &sPitch);
		multmat(&sysB, &Sched, &Bsch);
		// % Nominal Model State Equation
		multmat(&sysA, &x, &stEq1);
		multmat(&disC, &xd, &stEq2);
		multmat(&Bsch, &m_u, &stEq3);
		for(int i = 0; i < x.m; i++)
		{
			x.mat[i] = stEq1.mat[i] + stEq2.mat[i] + stEq3.mat[i];
		}
		// % Disturbance Model State Equation
		multmat(&disA, &xd, &dtEq1);
		for(int i = 0; i < xd.m; i++)
		{
			xd.mat[i] = dtEq1.mat[i];
		}
		// % Output Equation
		multmat(&sysC, &x, &m_y);
//		//fprintf (pFile, "Step 1 successfully %6.6f \n",m_y.mat[1]);
//		//fflush (pFile);
		// % 2) CORRECT
		// Observer Error
		for(int i = 0; i < oErr.m; i++)
		{
			oErr.mat[i] = m_ym.mat[i] - m_y.mat[i];
		}
		// % Update Nominal Model States
		multmat(&disD, &oErr, &suEq1);
		for(int i = 0; i < x.m; i++)
		{
			x.mat[i] = x.mat[i] + suEq1.mat[i];
		}
		// % Update Disturbance Model States
		multmat(&disB, &oErr, &duEq1);
		for(int i = 0; i < xd.m; i++)
		{
			xd.mat[i] = xd.mat[i] + duEq1.mat[i];
		}
		// % Update Output
		multmat(&sysC, &x, &m_y);
//		//fprintf (pFile, "Step 2 successfully %6.6f \n",m_y.mat[1]);
//		//fprintf (pFile, "Finished Round %d \n",myCount);
//		//fflush (pFile);
	}

}

extern "C" void __declspec(dllexport) __cdecl Control(double *position, double *velocity, double *action, int numAxis, double *wayPointX, double *wayPointY, int numWaypoints, double *actualWayPoint, double *param, int numParam)
{
	// % Controller Status Variable (Init = 0, Run = 1)
	static int cStatus = 0;
	// % Old Control Input
	static mat ukm1;
	// % Internal System Model State
	static mat x;
	// % Internal System Model A Matrix
	static mat sysA;
	// % Internal System Model B Matrix
	static mat sysB;
	// % Internal System Model C Matrix
	static mat sysC;
	// % Internal Disturbance Model State
	static mat xd;
	//% Internal Disturbance Model A Matrix
	static mat disA;
	// % Internal Disturbance Model B Matrix
	static mat disB;
	// % Internal Disturbance Model C Matrix
	static mat disC;
	// % Internal Disturbance Model D Matrix
	static mat disD;
	// % Roll Scheduling Parameters
	static mat sRoll;
	// % Pitch Scheduling Parameters
	static mat sPitch;
	// % Internal System Model C Position Matrix
	static mat sysCpos;
	// % Prediction Matrix Phi
	static mat pPhi;
	// % Prediction Matrix Jota
	static mat pJota;
	// % Prediction Matrix Theta
	static mat pTheta;
	// % Prediction Matrix Gamma
	static mat pGamma;
	// % Prediction Matrix Omega
	static mat pOmega;
	// % Prediction Matrix Psi
	static mat pPsi;
	// % Internal Time
	static doublereal iTime;
	// % Reference State
	static int res;
	// % Time when the last Waypoint was achieved
	static doublereal lastT;
	// % References
	static mat tr;
	static mat pr_x;
	static mat pr_y;
	// Measured Output in Matrix Form
	static mat m_ym;
	// Estimated Output in Matrix Form
	static mat m_y;


	// Temporal Matrices
	static mat sysCinv;
	static mat Sched;
	static mat Bsch;
	static mat stEq1;
	static mat stEq2;
	static mat stEq3;
	static mat oErr;
	static mat dtEq1;
	static mat suEq1;
	static mat duEq1;
	static mat rtmpx;
	static mat rtmpy;
	static mat mytime;
	static mat ref;
	static mat be;
	static mat ntheta;
	static mat ngamma;
	static mat ngammat;
	static mat G;
	static mat Ginv;
	static mat gEq1;
	static mat gEq2;
	static mat fEq1;
	static mat fEq2;
	static mat fEq3;
	static mat fEq4;
	static mat F;
	static mat du;

	// % Load Params
	doublereal dx, dy, vymax, aymax, jymax, vxmax, axmax, jxmax;
	vymax = param[0];
	aymax = param[1];
	jymax = param[2];
	vxmax = param[3];
	axmax = param[4];
	jxmax = param[5];

	// % If Observer Initialization
	if (cStatus == 0)
	{
		// % Initlialize Parameters and States
		sqlite3_initialize();
		sqlite3 *db;
		// Open Database
		int error = sqlite3_open_v2( "test.db", &db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL );
		if (error)
		{
			sqlite3_close(db);
			system("pause");
		}
		readmatsql(db, "Ap", &sysA);
		readmatsql(db, "Bp", &sysB);
		readmatsql(db, "Cp", &sysC);
		readmatsql(db, "Aj", &disA);
		readmatsql(db, "Bj", &disB);
		readmatsql(db, "Cj", &disC);
		readmatsql(db, "Dj", &disD);
		readmatsql(db, "rollparams", &sRoll);
		readmatsql(db, "pitchparams", &sPitch);

		readmatsql(db, "Cpos", &sysCpos);
		readmatsql(db, "phi", &pPhi);
		readmatsql(db, "jota", &pJota);
		readmatsql(db, "theta", &pTheta);
		readmatsql(db, "gama", &pGamma);
		readmatsql(db, "omega", &pOmega);
		readmatsql(db, "pssi", &pPsi);
		iTime = 0.0;
		res   = 0;
		lastT = 0.0;
		// Close Database
		sqlite3_close(db);
		// Init Inputs, States, Outputs and wayPoints in Matrix Form
		ukm1.m    = sysB.n;
		ukm1.n    = 1;
		ukm1.mat  = (doublereal*) calloc(ukm1.m*ukm1.n, sizeof(doublereal));
		m_ym.m   = sysC.m;
		m_ym.n   = 1;
		m_ym.mat = (doublereal*) calloc(m_ym.m*m_ym.n, sizeof(doublereal));
		m_y.m    = sysC.m;
		m_y.n    = 1;
		m_y.mat  = (doublereal*) calloc(m_y.m*m_y.n, sizeof(doublereal));
		x.m      = sysA.m;
		x.n      = 1;
		x.mat    = (doublereal*) calloc(x.m*x.n, sizeof(doublereal));
		xd.m     = disA.m;
		xd.n     = 1;
		xd.mat   = (doublereal*) calloc(xd.m*xd.n, sizeof(doublereal));
		// Init Temporal Matrices
		sysCinv.m   = sysC.n;
		sysCinv.n   = sysC.m;
		sysCinv.mat = (doublereal*) calloc(sysCinv.m*sysCinv.n, sizeof(doublereal));
		Bsch.m      = sysB.m;
		Bsch.n      = sysB.n;
		Bsch.mat    = (doublereal*) calloc(Bsch.m*Bsch.n, sizeof(doublereal));
		Sched.m     = sysB.n;
		Sched.n     = sysB.n;
		Sched.mat   = (doublereal*) calloc(Sched.m*Sched.n, sizeof(doublereal));
		stEq1.m     = sysA.m;
		stEq1.n     = x.n;
		stEq1.mat   = (doublereal*) calloc(stEq1.m*stEq1.n, sizeof(doublereal));
		stEq2.m     = disC.m;
		stEq2.n     = xd.n;
		stEq2.mat   = (doublereal*) calloc(stEq2.m*stEq2.n, sizeof(doublereal));
		stEq3.m     = Bsch.m;
		stEq3.n     = ukm1.n;
		stEq3.mat   = (doublereal*) calloc(stEq3.m*stEq3.n, sizeof(doublereal));
		oErr.m      = m_y.m;
		oErr.n      = m_y.n;
		oErr.mat    = (doublereal*) calloc(oErr.m*oErr.n, sizeof(doublereal));
		dtEq1.m     = disA.m;
		dtEq1.n     = xd.n;
		dtEq1.mat   = (doublereal*) calloc(dtEq1.m*dtEq1.n, sizeof(doublereal));
		suEq1.m     = disD.m;
		suEq1.n     = oErr.n;
		suEq1.mat   = (doublereal*) calloc(suEq1.m*suEq1.n, sizeof(doublereal));
		duEq1.m     = disB.m;
		duEq1.n     = oErr.n;
		duEq1.mat   = (doublereal*) calloc(duEq1.m*duEq1.n, sizeof(doublereal));
		rtmpx.m     = N;
		rtmpx.n     = 1;
		rtmpx.mat   = (doublereal*) calloc(rtmpx.m*rtmpx.n, sizeof(doublereal));
		rtmpy.m     = N;
		rtmpy.n     = 1;
		rtmpy.mat   = (doublereal*) calloc(rtmpy.m*rtmpy.n, sizeof(doublereal));
		mytime.m    = N;
		mytime.n    = 1;
		mytime.mat  = (doublereal*) calloc(mytime.m*mytime.n, sizeof(doublereal));
		ref.m       = 2*N;
		ref.n       = 1;
		ref.mat     = (doublereal*) calloc(ref.m*ref.n, sizeof(doublereal));
		be.m        = N*Bsch.m;
		be.n        = N*Bsch.n;
		be.mat      = (doublereal*) calloc(be.m*be.n, sizeof(doublereal));
		ntheta.m    = pTheta.m;
		ntheta.n    = Bsch.n;
		ntheta.mat  = (doublereal*) calloc(ntheta.m*ntheta.n, sizeof(doublereal));
		ngamma.m    = pGamma.m;
		ngamma.n    = be.n;
		ngamma.mat  = (doublereal*) calloc(ngamma.m*ngamma.n, sizeof(doublereal));
		ngammat.m   = ngamma.n;
		ngammat.n   = ngamma.m;
		ngammat.mat = (doublereal*) calloc(ngammat.m*ngammat.n, sizeof(doublereal));
		gEq1.m      = pOmega.m;
		gEq1.n      = ngamma.n;
		gEq1.mat    = (doublereal*) calloc(gEq1.m*gEq1.n, sizeof(doublereal));
		gEq2.m      = ngammat.m;
		gEq2.n      = gEq1.n;
		gEq2.mat    = (doublereal*) calloc(gEq2.m*gEq2.n, sizeof(doublereal));
		G.m         = ngammat.m;
		G.n         = ngamma.n;
		G.mat       = (doublereal*) calloc(G.m*G.n, sizeof(doublereal));
		Ginv.m      = G.n;
		Ginv.n      = G.m;
		Ginv.mat    = (doublereal*) calloc(Ginv.m*Ginv.n, sizeof(doublereal));
		fEq1.m      = ntheta.m;
		fEq1.n      = ukm1.n;
		fEq1.mat    = (doublereal*) calloc(fEq1.m*fEq1.n, sizeof(doublereal));
		fEq2.m      = pJota.m;
		fEq2.n      = xd.n;
		fEq2.mat    = (doublereal*) calloc(fEq2.m*fEq2.n, sizeof(doublereal));
		fEq3.m      = pPhi.m;
		fEq3.n      = x.n;
		fEq3.mat    = (doublereal*) calloc(fEq3.m*fEq3.n, sizeof(doublereal));
		fEq4.m      = pOmega.m;
		fEq4.n      = ref.n;
		fEq4.mat    = (doublereal*) calloc(fEq4.m*fEq4.n, sizeof(doublereal));
		du.m        = ukm1.m;
		du.n        = ukm1.n;
		du.mat      = (doublereal*) calloc(du.m*du.n, sizeof(doublereal));
		F.m         = ngammat.m;
		F.n         = fEq4.n;
		F.mat       = (doublereal*) calloc(F.m*F.n, sizeof(doublereal));

		// % Read Meas (Here just in Init, needed to estimate initial states)
		m_ym.mat[0] = position[0];
		m_ym.mat[1] = velocity[0];
		m_ym.mat[2] = position[1];
		m_ym.mat[3] = velocity[1];
//		// 0) Input/Output Pointers
//		zeros(&ukm1);
//		// Calculate C inverse
//		pinv(&(sysC.m), &(sysC.n), sysC.mat, &mytol, sysCinv.mat);
//		//   Multiply x = Cinv * ym
//		multmat(&sysCinv, &m_ym, &x);
//		//   Init xd = zeros
//		zeros(&xd);

		// % <----------- START CONTROL ----------------------------------------
		// % SUPERVISORY CONTROL
		// % Reference Design
		dx = wayPointX[res+1] - wayPointX[res];
		dy = wayPointY[res+1] - wayPointY[res];
		if ( abs(dy) > abs(0.8*(vymax/vxmax)*dx) )
		{
			thirdord( dy, vymax, aymax, jymax, Ts, &tr, &pr_y );
			if ( dy < 0 )
			{
				for (int i = 0; i < pr_y.m; i++)
				{
					pr_y.mat[i] = -pr_y.mat[i];
				}
			}
			pr_x.m   = pr_y.m;
			pr_x.n   = pr_y.n;
			pr_x.mat = (doublereal*) calloc((pr_x.m)*(pr_x.n), sizeof(doublereal));
			for (int i = 0; i < pr_y.m; i++)
			{
				pr_x.mat[i] = (dx/dy)*pr_y.mat[i];
			}
		}
		else if ( dx == 0.0 && dy == 0.0 )
		{
			tr.m = 2;
			tr.n = 1;
			tr.mat = (doublereal*) calloc((tr.m)*(tr.n), sizeof(doublereal));
			tr.mat[0] = 0.0;
			tr.mat[1] = Ts;
			pr_x.m = 2;
			pr_x.n = 1;
			pr_x.mat = (doublereal*) calloc((pr_x.m)*(pr_x.n), sizeof(doublereal));
			pr_x.mat[0] = 0.0;
			pr_x.mat[1] = 0.0;
			pr_y.m = 2;
			pr_y.n = 1;
			pr_y.mat = (doublereal*) calloc((pr_y.m)*(pr_y.n), sizeof(doublereal));
			pr_y.mat[0] = 0.0;
			pr_y.mat[1] = 0.0;
		}
		else
		{
			thirdord( dx, vxmax, axmax, jxmax, Ts, &tr, &pr_x );
			if ( dx < 0 )
			{
				for (int i = 0; i < pr_x.m; i++)
				{
					pr_x.mat[i] = -pr_x.mat[i];
				}
			}
			pr_y.m   = pr_x.m;
			pr_y.n   = pr_x.n;
			pr_y.mat = (doublereal*) calloc((pr_y.m)*(pr_y.n), sizeof(doublereal));
			for (int i = 0; i < pr_y.m; i++)
			{
				pr_y.mat[i] = (dy/dx)*pr_x.mat[i];
			}
		}
		for (int i = 0; i < tr.m; i++)
		{
			tr.mat[i] = tr.mat[i] + lastT;
			pr_x.mat[i] = pr_x.mat[i] + wayPointX[res];
			pr_y.mat[i] = pr_y.mat[i] + wayPointY[res];
		}
		// % 1) PREDICT
		// % Scheduled B matrix
		zeros(&Sched);
		Sched.mat[0] = sinfit(ukm1.mat[0], &sRoll);
		Sched.mat[3] = sinfit(ukm1.mat[1], &sPitch);
		multmat(&sysB, &Sched, &Bsch);
		// % Nominal Model State Equation
		multmat(&sysA, &x, &stEq1);
		multmat(&disC, &xd, &stEq2);
		multmat(&Bsch, &ukm1, &stEq3);
		for(int i = 0; i < x.m; i++)
		{
			x.mat[i] = stEq1.mat[i] + stEq2.mat[i] + stEq3.mat[i];
		}
		// % Disturbance Model State Equation
		multmat(&disA, &xd, &dtEq1);
		for(int i = 0; i < xd.m; i++)
		{
			xd.mat[i] = dtEq1.mat[i];
		}
		// % Output Equation
		multmat(&sysC, &x, &m_y);
		// % 2) CORRECT
		// Observer Error
		for(int i = 0; i < oErr.m; i++)
		{
			oErr.mat[i] = m_ym.mat[i] - m_y.mat[i];
		}
		// % Update Nominal Model States
		multmat(&disD, &oErr, &suEq1);
		for(int i = 0; i < x.m; i++)
		{
			x.mat[i] = x.mat[i] + suEq1.mat[i];
		}
		// % Update Disturbance Model States
		multmat(&disB, &oErr, &duEq1);
		for(int i = 0; i < xd.m; i++)
		{
			xd.mat[i] = xd.mat[i] + duEq1.mat[i];
		}
		// % Update Output
		multmat(&sysC, &x, &m_y);
		// % STATE FEEDBACK
		// % 1) Compute Reference Vector
		for ( int i = 0; i < N; i++ )
		{
			mytime.mat[i] = iTime + i*Ts;
		}
		interpola( &tr, &pr_x, &mytime, &rtmpx );
		interpola( &tr, &pr_y, &mytime, &rtmpy );
		for ( int i = 0; i < N; i++ )
		{
			ref.mat[2*i]     = rtmpx.mat[i];
			ref.mat[(2*i)+1] = rtmpy.mat[i];
		}
		// % 2) Compute Unconstrained MPC Input
		int filas    = Bsch.m;
		int columnas = Bsch.n;
		for ( int i = 0; i < N; i++ )
		{
			for ( int j = 0; j < columnas; j++ )
			{
				for ( int k = 0; k < filas; k++ )
				{
					be.mat[i*(N*filas*columnas+filas) + j*N*filas + k] = Bsch.mat[j*filas + k];
				}
			}
		}
		multmat(&pTheta, &Bsch, &ntheta);
		multmat(&pGamma, &be, &ngamma);
		// transpose of matrix ngamma
	    for (int i = 0; i < ngamma.m; i++) //stored in column major
	    {
	  	    for(int j = 0; j < ngamma.n; j++)
		    {
			    ngammat.mat[(ngamma.n)*(i)+j] = ngamma.mat[j*(ngamma.m)+i];
		    }
	    }
	    // G matrix
	    multmat(&pOmega, &ngamma, &gEq1);
	    multmat(&ngammat, &gEq1, &gEq2);
	    for (int i = 0; i < G.m*G.n; i++ )
	    {
			G.mat[i] = 2.0*(pPsi.mat[i] + gEq2.mat[i]);
	    }
	    // F matrix
	    multmat(&ntheta, &ukm1, &fEq1);
	    multmat(&pJota, &xd, &fEq2);
	    multmat(&pPhi, &x, &fEq3);
	    for (int i = 0; i < fEq1.m*fEq1.n; i++ )
	    {
	    	fEq1.mat[i] = fEq1.mat[i] + fEq2.mat[i] + fEq3.mat[i] - ref.mat[i];
	    }
	    multmat(&pOmega, &fEq1, &fEq4);
	    multmat(&ngammat, &fEq4, &F);
	    for ( int i = 0; i < F.m*F.n; i++ )
	    {
	    	F.mat[i] = 2.0*F.mat[i];
	    }
	    // Solution
	    myls(&G, &F);
		for ( int i = 0; i < sysB.n; i++ )
		{
			// Keep du in case we want to limit it
			du.mat[i] = F.mat[i];
			ukm1.mat[i] = ukm1.mat[i] - du.mat[i];
		}
	    // % 3) Limit the Inputs
	    ukm1.mat[0] = min(ukm1.mat[0],1.0);
	    ukm1.mat[0] = max(ukm1.mat[0],-1.0);
	    ukm1.mat[1] = min(ukm1.mat[1],1.0);
		ukm1.mat[1] = max(ukm1.mat[1],-1.0);
		// % 4) Write Input
		for ( int i = 0; i < sysB.n; i++ )
		{
			action[i] = ukm1.mat[i];
		}
		// Next internal time
	    iTime = iTime + Ts;

		// % Initialization Finished
		cStatus = 1;
	}
	else if (cStatus == 1) // <---------------- RUNNING ---------------------------------
	{
		// SUPERVISORY CONTROL
		// Check weather next point has been reached
		if ( sqrt( pow((position[0]-wayPointX[res+1]), 2.0) + pow((position[1]-wayPointY[res+1]), 2.0) ) < radio )
		{
			if ( (res+2) < numWaypoints ) // red+2 still because later res++ and access res+1
			{
				// Free previous reference
				free(tr.mat);
				free(pr_x.mat);
				free(pr_y.mat);
				// % Update Reference State and Last Time
				res   = res + 1;
				lastT = iTime;
				// % Reference Design
				dx = wayPointX[res+1] - wayPointX[res];
				dy = wayPointY[res+1] - wayPointY[res];
				if ( abs(dy) > abs(0.8*(vymax/vxmax)*dx) )
				{
					thirdord( dy, vymax, aymax, jymax, Ts, &tr, &pr_y );
					if ( dy < 0 )
					{
						for (int i = 0; i < pr_y.m; i++)
						{
							pr_y.mat[i] = -pr_y.mat[i];
						}
					}
					pr_x.m   = pr_y.m;
					pr_x.n   = pr_y.n;
					pr_x.mat = (doublereal*) calloc((pr_x.m)*(pr_x.n), sizeof(doublereal));
					for (int i = 0; i < pr_y.m; i++)
					{
						pr_x.mat[i] = (dx/dy)*pr_y.mat[i];
					}
				}
				else if ( dx == 0.0 && dy == 0.0 )
				{
					tr.m = 2;
					tr.n = 1;
					tr.mat = (doublereal*) calloc((tr.m)*(tr.n), sizeof(doublereal));
					tr.mat[0] = 0.0;
					tr.mat[1] = Ts;
					pr_x.m = 2;
					pr_x.n = 1;
					pr_x.mat = (doublereal*) calloc((pr_x.m)*(pr_x.n), sizeof(doublereal));
					pr_x.mat[0] = 0.0;
					pr_x.mat[1] = 0.0;
					pr_y.m = 2;
					pr_y.n = 1;
					pr_y.mat = (doublereal*) calloc((pr_y.m)*(pr_y.n), sizeof(doublereal));
					pr_y.mat[0] = 0.0;
					pr_y.mat[1] = 0.0;
				}
				else
				{
					thirdord( dx, vxmax, axmax, jxmax, Ts, &tr, &pr_x );
					if ( dx < 0 )
					{
						for (int i = 0; i < pr_x.m; i++)
						{
							pr_x.mat[i] = -pr_x.mat[i];
						}
					}
					pr_y.m   = pr_x.m;
					pr_y.n   = pr_x.n;
					pr_y.mat = (doublereal*) calloc((pr_y.m)*(pr_y.n), sizeof(doublereal));
					for (int i = 0; i < pr_y.m; i++)
					{
						pr_y.mat[i] = (dy/dx)*pr_x.mat[i];
					}
				}
				for (int i = 0; i < tr.m; i++)
				{
					tr.mat[i] = tr.mat[i] + lastT;
					pr_x.mat[i] = pr_x.mat[i] + wayPointX[res];
					pr_y.mat[i] = pr_y.mat[i] + wayPointY[res];
				}
			}
			else // Ruta Completed !
			{
				// Free previous reference
				free(tr.mat     );
				free(pr_x.mat   );
				free(pr_y.mat   );
				// Free controller matrices
				free(ukm1.mat   );
				free(m_ym.mat   );
				free(m_y.mat    );
				free(x.mat      );
				free(xd.mat     );
				free(sysCinv.mat);
				free(Bsch.mat   );
				free(Sched.mat  );
				free(stEq1.mat  );
				free(stEq2.mat  );
				free(stEq3.mat  );
				free(oErr.mat   );
				free(dtEq1.mat  );
				free(suEq1.mat  );
				free(duEq1.mat  );
				free(rtmpx.mat  );
				free(rtmpy.mat  );
				free(mytime.mat );
				free(ref.mat    );
				free(be.mat     );
				free(ntheta.mat );
				free(ngamma.mat );
				free(ngammat.mat);
				free(gEq1.mat   );
				free(gEq2.mat   );
				free(G.mat      );
				free(Ginv.mat   );
				free(fEq1.mat   );
				free(fEq2.mat   );
				free(fEq3.mat   );
				free(fEq4.mat   );
				free(du.mat     );
				free(F.mat      );
				// % Status End
				cStatus = -1;
				// % Write Input
				for ( int i = 0; i < sysB.n; i++ )
				{
					action[i] = 0.0;
				}
				// Quick Return
				return;
			}
		}

		// % Read Meas
		m_ym.mat[0] = position[0];
		m_ym.mat[1] = velocity[0];
		m_ym.mat[2] = position[1];
		m_ym.mat[3] = velocity[1];
		// % 1) PREDICT
		// % Scheduled B matrix
		zeros(&Sched);
		Sched.mat[0] = sinfit(ukm1.mat[0], &sRoll);
		Sched.mat[3] = sinfit(ukm1.mat[1], &sPitch);
		multmat(&sysB, &Sched, &Bsch);
		// % Nominal Model State Equation
		multmat(&sysA, &x, &stEq1);
		multmat(&disC, &xd, &stEq2);
		multmat(&Bsch, &ukm1, &stEq3);
		for(int i = 0; i < x.m; i++)
		{
			x.mat[i] = stEq1.mat[i] + stEq2.mat[i] + stEq3.mat[i];
		}
		// % Disturbance Model State Equation
		multmat(&disA, &xd, &dtEq1);
		for(int i = 0; i < xd.m; i++)
		{
			xd.mat[i] = dtEq1.mat[i];
		}
		// % Output Equation
		multmat(&sysC, &x, &m_y);

		// % 2) CORRECT
		// Observer Error
		for(int i = 0; i < oErr.m; i++)
		{
			oErr.mat[i] = m_ym.mat[i] - m_y.mat[i];
		}
		// % Update Nominal Model States
		multmat(&disD, &oErr, &suEq1);
		for(int i = 0; i < x.m; i++)
		{
			x.mat[i] = x.mat[i] + suEq1.mat[i];
		}
		// % Update Disturbance Model States
		multmat(&disB, &oErr, &duEq1);
		for(int i = 0; i < xd.m; i++)
		{
			xd.mat[i] = xd.mat[i] + duEq1.mat[i];
		}
		// % Update Output
		multmat(&sysC, &x, &m_y);
		// % STATE FEEDBACK
		// % 1) Compute Reference Vector
		for ( int i = 0; i < N; i++ )
		{
			mytime.mat[i] = iTime + i*Ts;
		}
		interpola( &tr, &pr_x, &mytime, &rtmpx );
		interpola( &tr, &pr_y, &mytime, &rtmpy );
		for ( int i = 0; i < N; i++ )
		{
			ref.mat[2*i]     = rtmpx.mat[i];
			ref.mat[(2*i)+1] = rtmpy.mat[i];
		}
		// % 2) Compute Unconstrained MPC Input
		int filas    = Bsch.m;
		int columnas = Bsch.n;
		for ( int i = 0; i < N; i++ )
		{
			for ( int j = 0; j < columnas; j++ )
			{
				for ( int k = 0; k < filas; k++ )
				{
					be.mat[i*(N*filas*columnas+filas) + j*N*filas + k] = Bsch.mat[j*filas + k];
				}
			}
		}
		multmat(&pTheta, &Bsch, &ntheta);
		multmat(&pGamma, &be, &ngamma);
		// transpose of matrix ngamma
		for (int i = 0; i < ngamma.m; i++) //stored in column major
		{
			for(int j = 0; j < ngamma.n; j++)
			{
				ngammat.mat[(ngamma.n)*(i)+j] = ngamma.mat[j*(ngamma.m)+i];
			}
		}
		// G matrix
		multmat(&pOmega, &ngamma, &gEq1);
		multmat(&ngammat, &gEq1, &gEq2);
		for (int i = 0; i < G.m*G.n; i++ )
		{
			G.mat[i] = 2.0*(pPsi.mat[i] + gEq2.mat[i]);
		}
		// F matrix
		multmat(&ntheta, &ukm1, &fEq1);
		multmat(&pJota, &xd, &fEq2);
		multmat(&pPhi, &x, &fEq3);
		for (int i = 0; i < fEq1.m*fEq1.n; i++ )
		{
			fEq1.mat[i] = fEq1.mat[i] + fEq2.mat[i] + fEq3.mat[i] - ref.mat[i];
		}
		multmat(&pOmega, &fEq1, &fEq4);
		multmat(&ngammat, &fEq4, &F);
		for ( int i = 0; i < F.m*F.n; i++ )
		{
			F.mat[i] = 2.0*F.mat[i];
		}
		// Solution
		myls(&G, &F);
		for ( int i = 0; i < sysB.n; i++ )
		{
			// Keep du in case we want to limit it
			du.mat[i] = F.mat[i];
			ukm1.mat[i] = ukm1.mat[i] - du.mat[i];
		}
		// 3) Limit the Inputs
		ukm1.mat[0] = min(ukm1.mat[0],1.0);
		ukm1.mat[0] = max(ukm1.mat[0],-1.0);
		ukm1.mat[1] = min(ukm1.mat[1],1.0);
		ukm1.mat[1] = max(ukm1.mat[1],-1.0);
		// % 4) Write Input
		for ( int i = 0; i < sysB.n; i++ )
		{
			action[i] = ukm1.mat[i];
		}
		// Next internal time
		iTime = iTime + Ts;
	}
	else // Ruta Completed and Memory was Released Already
	{
		// Write Input
		for ( int i = 0; i < sysB.n; i++ )
		{
			action[i] = 0.0;
		}
	}

}

void zeros( mat *matrix )
{
	for (int i = 0; i < (matrix->m)*(matrix->n); i++)
	    {
			matrix->mat[i] = 0.0;
	    }
}

void readmatsql( sqlite3 *db, char *matname, mat *matrix )
{
	// Read Matrix Sizes from 'info' table
	sqlite3_stmt *res1;
	char querybegin1[] = "SELECT * FROM info WHERE matname = \"";
	char *query1;
	query1 = (char*) malloc(strlen(querybegin1)+strlen(matname)+2);
	strcpy(query1, querybegin1);
	strcat(query1, matname);
	strcat(query1, "\"");
	int error = sqlite3_prepare_v2(db,query1,-1,&res1,NULL);
	if (error != SQLITE_OK)
	{
		sqlite3_close(db);
		puts("query 1 prepare error");
		puts(query1);
		system("pause");
	}
	sqlite3_step(res1);
	matrix->m = sqlite3_column_int(res1,1);
	matrix->n = sqlite3_column_int(res1,2);
	sqlite3_finalize(res1);
	// Allocate Matrix Memory
	matrix->mat = (doublereal*) calloc(matrix->m*matrix->n, sizeof(doublereal));
	// Read Matrix Values
	sqlite3_stmt *res2;
	char querybegin2[] = "SELECT * FROM \"";
	char *query2;
	query2 = (char*) malloc(strlen(querybegin2)+strlen(matname)+2);
	strcpy(query2, querybegin2);
	strcat(query2, matname);
	strcat(query2, "\"");
	error = sqlite3_prepare_v2(db,query2,-1,&res2,NULL);
	if (error != SQLITE_OK)
	{
		sqlite3_close(db);
		puts("query 2 prepare error");
		puts(query2);
		system("pause");
	}
	int i = 0;
	while (sqlite3_step(res2) == SQLITE_ROW)
	{
		for(int j = 0; j < matrix->n; j++)
		{
			matrix->mat[j*(matrix->m)+i] = sqlite3_column_double(res2,j);
		}
		i = i + 1;
	}
	sqlite3_finalize(res2);

}

void monpennew(integer *m, integer *n, doublereal *a, doublereal *mytol, doublereal *ainvt,
		doublereal *work, doublereal *s, doublereal *u, doublereal *vt, doublereal *sfull, doublereal *usfull)
{
	// init tolerance
    doublereal tol = 0.0;
    // define tolerance
//    if (*mytol <= 0)
//    {
//		// tolerance
//		tol = (*m) * machEps;
//    }
//    else
//    {
    	tol = *mytol;
//    }
    // allocate work matrix
    integer lwork = 5*(*m)*(*n);
//    doublereal *work = (doublereal*) calloc(5*(*m)*(*n), sizeof(doublereal));
    // allocate s matrix
//    doublereal *s = (doublereal*) calloc(*n, sizeof(doublereal)); // *m > *n for sure
    // allocate u matrix
//    doublereal *u = (doublereal*) calloc((*m)*(*m), sizeof(doublereal));
    // allocate vt matrix
//    doublereal *vt = (doublereal*) calloc((*n)*(*n), sizeof(doublereal));
    // dgesvd function output
    integer info  = 0.0;
    // singular value decomposition
    char jobu      = 'A';
    char jobvt     = 'A';
    dgesvd_(&jobu, &jobvt, m, n, a,
		    m, s, u, m, vt, n,
		    work, &lwork, &info);
    // allocate sfull matrix
//    doublereal *sfull = (doublereal*) calloc((*m)*(*n), sizeof(doublereal));
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
//    doublereal *usfull = (doublereal*) calloc((*m)*(*n), sizeof(doublereal));
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

//	free (work);
//    free (s);
//    free (u);
//    free (vt);
//    free (sfull);
//    free (usfull);
}

void pinvnew(integer *m, integer *n, doublereal *a, doublereal *mytol, doublereal *ainv,
		doublereal *work, doublereal *s, doublereal *u, doublereal *vt, doublereal *sfull, doublereal *usfull, doublereal *at, doublereal *ainvt)
{
	  // moore–penrose pseudoinverse
	  if (*n > *m)
	  {
		  // allocate at matrix
//		  doublereal *at = (doublereal*) calloc((*n)*(*m), sizeof(doublereal));
		  // transpose of matrix a
		  for (int i = 0; i < *m; i++) //stored in column major
		  {
			  for(int j = 0; j < *n; j++)
			  {
				  at[(*n)*(i)+j] = a[j*(*m)+i];
			  }
		  }
		  // apply same algorithm to a^t and transpose the result
		  monpennew(n,m,at,mytol,ainv,work,s,u,vt,sfull,usfull);
		  // free at matrix
//		  free (at);
	  }
	  else
	  {
		  // allocate ainvt matrix
//		  doublereal *ainvt = (doublereal*) calloc((*m)*(*n), sizeof(doublereal));
		  // apply normal algorithm to a
		  monpennew(m,n,a,mytol,ainvt,work,s,u,vt,sfull,usfull);
		  // transpose of matrix ainvt
		  for (int i = 0; i < *m; i++) //stored in column major
		  {
			  for(int j = 0; j < *n; j++)
			  {
				  ainv[(*n)*(i)+j] = ainvt[j*(*m)+i];
			  }
		  }
		  // free ainvt matrix
//		  free (ainvt);
	  }
}

void multmat( mat *A, mat *B, mat *C )
{
	char transa      = 'N';
    char transb      = 'N';
    doublereal alpha = 1.0;
    doublereal beta  = 0.0;

    dgemm_(&transa, &transb,
	  	   &A->m, &B->n, &A->n, &alpha, A->mat, &A->m,
		   B->mat, &A->n, &beta, C->mat, &A->m);

//	/*           Form  C := A*B */
//
//if (C->m > C->n) // size(A) > size(B) -> A might contain more zeros
//{
//	for (int j = 0; j < A->m; j++) //m
//		{
//			for (int i = 0; i < C->n; i++)
//			{
//				C->mat[(i * C->m) + j] = 0.0; // Clean column first
//			}
//			for (int l = 0; l < A->n; l++) //k
//			{
//				if (A->mat[(l * A->m) + j] != 0.0)
//				{
//					for (int i = 0; i < B->n; i++) //m
//					{
//						C->mat[(i * C->m) + j] += B->mat[(i * B->m) + l] * A->mat[(l * A->m) + j];
//					}
//				}
//			}
//		}
//}
//else // Apparently not working
//{
//		for (int j = 0; j < B->n; j++) //n
//		{
//			for (int i = 0; i < C->m; i++)
//			{
//				C->mat[i + (j * C->m)] = 0.0; // Clean column first
//			}
//			for (int l = 0; l < A->n; l++) //k
//			{
//				if (B->mat[l + j * B->m] != 0.0)
//				{
//					for (int i = 0; i < A->m; i++) //m
//					{
//						C->mat[i + (j * C->m)] += A->mat[i + (l * A->m)] * B->mat[l + (j * B->m)];
//					}
//				}
//			}
//		}
//}

}

void myls( mat *R, mat *F)
{
	integer info;
	char trans      = 'N';
	// allocate work matrix
	integer lwork = 5*(R->m)*(R->n);
	doublereal *work = (doublereal*) calloc(lwork, sizeof(doublereal));

	dgels_(&trans, &R->m, &R->n, &F->n, R->mat, &R->m, F->mat, &F->m, work, &lwork, &info);

	free(work);
}

doublereal sinfit( doublereal u, mat *fit_params )
{
	int i = 0;
	doublereal out = 0.0;
	while( i < max(fit_params->m,fit_params->n) )
	{
		out = out + (fit_params->mat[i])*sin( fit_params->mat[i+1]*u + fit_params->mat[i+2] );
		i = i + 3;
	}
	return out;
}

void thirdord( doublereal p, doublereal v, doublereal a, doublereal j, doublereal Ts, mat *tx, mat *xp )
{
	// Init Variables
	doublereal t1;
	doublereal jd;
	doublereal t2;
	doublereal t3;
	mat t;
	t.m   = 1;
	t.n   = 3;
	t.mat = (doublereal*) calloc((t.m)*(t.n), sizeof(doublereal));
	zeros(&t);
	mat tconst;
	tconst.m   = 3;
	tconst.n   = 8;
//	tconst.mat = (doublereal*) calloc((tconst.m)*(tconst.n), sizeof(doublereal));
	doublereal t_const[] = {0,0,0, 1,0,0, 1,1,0, 2,1,0, 2,1,1, 3,1,1, 3,2,1, 4,2,1};
	tconst.mat = &t_const[0];
	mat tt;
	tt.m = t.m;
	tt.n = tconst.n;
	tt.mat = (doublereal*) calloc((tt.m)*(tt.n), sizeof(doublereal));
	zeros(&tt);
	mat ttest;
	ttest.m = tt.m;
	ttest.n = tt.n + 1;
	ttest.mat = (doublereal*) calloc((ttest.m)*(ttest.n), sizeof(doublereal));
	zeros(&ttest);
	int len;
	mat xj;
	xj.n = 1;
	mat xa;
	xa.n = 1;
	mat xv;
	xv.n = 1;

	//fprintf (pFile, "THIRDORD: Init %d \n",66);
	//fflush (pFile);

	// % PART 1
	p = abs(p);
	v = abs(v);
	a = abs(a);
	j = abs(j);

	//fprintf (pFile, "THIRDORD: p = %6.6f, v = %6.6f, a = %6.6f, j = %6.6f, Ts = %6.6f \n",p,v,a,j,Ts);
	//fflush (pFile);

	// % Calculation t1
    t1 = pow((p/(2.0*j)),(1.0/3.0)) ; // % largest t1 with bound on jerk
    //fprintf (pFile, "THIRDORD: t1 = %6.6f \n",t1);
	//fflush (pFile);
    t1 = ceil(t1/Ts)*Ts;
    //fprintf (pFile, "THIRDORD: t1 = %6.6f \n",t1);
	//fflush (pFile);
    jd = 1.0/2.0*p/(pow(t1,3.0));
    // % velocity test
	if ( v < jd*(pow(t1,2.0)) )         // % v bound violated ?
	{
	   t1 = pow((v/j),(1.0/2.0)) ;  // % t1 with bound on velocity not violated
	   //fprintf (pFile, "THIRDORD: t1 = %6.6f \n",t1);
	   	//fflush (pFile);
	   t1 = ceil(t1/Ts)*Ts;
	   //fprintf (pFile, "THIRDORD: t1 = %6.6f \n",t1);
	   	//fflush (pFile);
	   jd = v/(pow(t1,2.0));
	}
	// % acceleration test
	if (a < jd*t1)     // % a bound violated ?
	{
	   t1 = a/j ;    // % t1 with bound on acceleration not violated
	   //fprintf (pFile, "THIRDORD: t1 = %6.6f \n",t1);
	   	//fflush (pFile);
	   t1 = ceil(t1/Ts)*Ts;
	   //fprintf (pFile, "THIRDORD: t1 = %6.6f \n",t1);
	   	//fflush (pFile);
	   jd = a/t1;
	}
	j = jd;  // % as t1 is now fixed, jd is the new bound on jerk
	//fprintf (pFile, "THIRDORD: t1 = %6.6f, jd = %6.6f \n",t1,jd);
	//fflush (pFile);

	// % Calculation t2
	t2 = pow((pow(t1,2.0)/4.0+p/j/t1),(1.0/2.0)) - (3.0/2.0)*t1 ;   // % largest t2 with bound on acceleration
	//fprintf (pFile, "THIRDORD: t2 = %6.6f \n",t2);
	//fflush (pFile);
	t2 = ceil(t2/Ts)*Ts;
	//fprintf (pFile, "THIRDORD: t2 = %6.6f \n",t2);
	//fflush (pFile);
	jd = p/( 2.0*(pow(t1,3.0)) + 3.0*(pow(t1,2.0))*t2 + t1*(pow(t2,2.0)) );

	// % velocity test
	if ( v < (jd*pow(t1,2.0) + jd*t1*t2) )   // % v bound violated ?
	{
	   t2 = v/(j*t1) - t1 ;       // % t2 with bound on velocity not violated
	   //fprintf (pFile, "THIRDORD: t2 = %6.6f \n",t2);
	   //fflush (pFile);
	   t2 = ceil(t2/Ts)*Ts;
	   jd = v/( pow(t1,2.0) + t1*t2 );
	}
	j = jd;  // % as t2 is now fixed, jd is the new bound on jerk
	//fprintf (pFile, "THIRDORD: t2 = %6.6f, jd = %6.6f \n",t2,jd);
	//fflush (pFile);

	// % Calculation t3
	t3 = ( p - 2.0*j*(pow(t1,3.0)) - 3.0*j*(pow(t1,2.0))*t2 - j*t1*(pow(t2,2.0)) )/v ; // % t3 with bound on velocity
	//fprintf (pFile, "THIRDORD: t3 = %6.6f \n",t3);
	   //fflush (pFile);
	t3 = ceil(t3/Ts)*Ts;
	//fprintf (pFile, "THIRDORD: t3 = %6.6f \n",t3);
   //fflush (pFile);
	jd = p/( 2.0*(pow(t1,3.0)) + 3.0*(pow(t1,2.0))*t2 + t1*(pow(t2,2.0)) + (pow(t1,2.0))*t3 + t1*t2*t3 );

	// % All time intervals are now calculated
	t.mat[0] = t1;
	t.mat[1] = t2;
	t.mat[2] = t3;

	// % PART 2
	//fprintf (pFile, "THIRDORD: t1 = %6.6f, t2 = %6.6f, t3 = %6.6f, jd = %6.6f \n",t1,t2,t3,jd);
	//fprintf (pFile, "THIRDORD: Part 1 %d \n",66);
	//fflush (pFile);

	multmat(&t, &tconst, &tt);

	for ( int i = 0; i < tt.n; i++ )
	{
		ttest.mat[i] = tt.mat[i];
		//fprintf (pFile, "%6.6f, ",ttest.mat[i]);
	}
	ttest.mat[tt.n] = 1.5*tt.mat[7];
	//fprintf (pFile, "%6.6f \n ",ttest.mat[tt.n]);
	//fprintf (pFile, "THIRDORD: tt_end: %6.6f \n",tt.mat[7]);
	//fflush (pFile);

	len = (int) lround(1.2*tt.mat[7]/Ts + 1.0);
	//fprintf (pFile, "THIRDORD: len: %d \n",len);
	//fflush (pFile);
	xj.m = len;
	xj.mat = (doublereal*) calloc((xj.m)*(xj.n), sizeof(doublereal));
	zeros(&xj);
	xa.m = len;
	xa.mat = (doublereal*) calloc((xa.m)*(xa.n), sizeof(doublereal));
	zeros(&xa);
	xv.m = len;
	xv.mat = (doublereal*) calloc((xv.m)*(xv.n), sizeof(doublereal));
	zeros(&xv);
	xj.mat[0] = jd;
	tx->m = len;
	tx->n = 1;
	tx->mat = (doublereal*) calloc((tx->m)*(tx->n), sizeof(doublereal));
	xp->m = len;
	xp->n = 1;
	xp->mat = (doublereal*) calloc((xp->m)*(xp->n), sizeof(doublereal));
	for ( int i = 0; i < tx->m; i++ )
	{
		tx->mat[i] = i*Ts;
		//fprintf (pFile, "%6.6f, ",tx->mat[i]);
	}
	//fprintf (pFile, "fin %d \n",666);
	//fflush (pFile);


	int i;
	for ( int k = 0; k < (len-1); k++ )
	{
		i = -1;
		for ( int j = 0; j < ttest.n; j++ )
		{
			if ( Ts*((k+1) + 0.5) <= ttest.mat[j] )
			{
				i = j - 1;
				break;
			}
		}
		if ( i == 0 || i == 6 )
		{
			xj.mat[k+1] = jd;
		}
		else if ( i == 2 || i == 4 )
		{
			xj.mat[k+1] = -jd;
		}
		else
		{
			xj.mat[k+1] = 0;
		}
		xa.mat[k+1] = xa.mat[k] + xj.mat[k]*Ts;
		xv.mat[k+1] = xv.mat[k] + xa.mat[k]*Ts;
		xp->mat[k+1] = xp->mat[k] + xv.mat[k]*Ts;
		//fprintf (pFile, "%6.6f, ",xp->mat[k+1]);
	}
	//fprintf (pFile, " fin %d \n",666);
	//fflush (pFile);

	//fprintf (pFile, "THIRDORD: Part 2 %d \n",66);
	//fflush (pFile);

	free(xj.mat);
	free(xa.mat);
	free(xv.mat);
	free(t.mat);
	free(tt.mat);
	free(ttest.mat);

	//fprintf (pFile, "THIRDORD: Free Stuff %d \n",66);
	//fflush (pFile);

}

void interpola ( mat *tr_x, mat *pr_x, mat *my_tr, mat *my_pr )
{
	int tr_mem;
	int my_flag;
	doublereal m;
	doublereal b;

	zeros(my_pr);
	tr_mem = 0;
	for ( int i = 0; i < my_tr->m; i++ )
	{
		my_flag = 0;
		for ( int j = tr_mem; j < tr_x->m; j++ )
		{
			if ( my_tr->mat[i] < tr_x->mat[j] )
			{
				m = (pr_x->mat[j]-pr_x->mat[j-1])/(tr_x->mat[j]-tr_x->mat[j-1]);
				b = pr_x->mat[j] - m*(tr_x->mat[j]);
				my_pr->mat[i] = m*(my_tr->mat[i]) + b;
				my_flag = 1;
				tr_mem = j;
				break;
			}
			else if ( my_tr->mat[i] == tr_x->mat[j] )
			{
				my_pr->mat[i] = pr_x->mat[j];
				my_flag = 1;
				tr_mem = j;
				break;
			}
			else
			{
				continue;
			}
		}
		if ( my_flag == 1 )
		{
			continue;
		}
		else
		{
			my_pr->mat[i] = pr_x->mat[pr_x->m - 1];
		}
	}
}
