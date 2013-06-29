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

#define EXPORTING_DLL

// to avoid C++ incompatibilities (after all, it is expecting C)
extern "C" {               // make sure these are visible in directory lookup
#include "./lapack/f2c.h"           // our solution: copy them to /usr/include
#include "./lapack/clapack.h"
#include "./sqlite/sqlite3.h"
}

#include "mydll.h"

//doublereal machEps = 1.0;

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

extern "C" void __declspec(dllexport) __cdecl sqltest(integer *m, integer *n, doublereal *mymat)
{
	sqlite3_initialize();

	sqlite3 *db;
	sqlite3_stmt *res;

	const char *tail;

	int error = sqlite3_open_v2( "test.db", &db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL );

	if (error)
	{
		sqlite3_close(db);
		system("pause");
	}

	char query[] = "SELECT * FROM Aj";

	error = sqlite3_prepare_v2(db,&query[0],sizeof(query),&res,&tail);

	if (error != SQLITE_OK)
	{
		sqlite3_close(db);
		system("pause");
	}

	int i = 0;
	while (sqlite3_step(res) == SQLITE_ROW)
	{
		for(int j = 0; j < *n; j++)
	    {
		    mymat[j*(*m)+i] = sqlite3_column_double(res,j);
	    }
		i = i + 1;
	}

	sqlite3_finalize(res);

	sqlite3_close(db);

}

extern "C" void __declspec(dllexport) __cdecl matobserver(doublereal *u, doublereal *ym, doublereal *y)
{
	// Debugging Log
//	static FILE * pFile;
//	static int myCount = 0;
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

	// Calculate Machine Epsilon
	// linear search for machine tolerance (**eliminate this when possible**)
//	do {
//	   machEps /= 2.0;
//	}
//	while ((doublereal)(1.0 + (machEps/2.0)) != 1.0);

	// % If Observer Initialization
	if (cStatus == 0)
	{
//		pFile = fopen ("log.txt","w");
//		fprintf (pFile, "INIT %d \n",0);
//		fflush (pFile);
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
//		fprintf (pFile, "Matrices loaded successfully %d \n",1);
//		fflush (pFile);
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
//		fprintf (pFile, "Matrices inited successfully %d \n",2);
//		fflush (pFile);
		// 0) Input/Output Pointers
		m_u.mat  = u;
		m_ym.mat = ym;
		m_y.mat  = y;
//		fprintf (pFile, "Step 0 successfully %d \n",4);
//		fflush (pFile);
		// % Initial Conditions
		// Initial Plant States
		//   Calculate C inverse
		// doublereal mytol = -1.0;
		doublereal mytol = 0.0;
		pinv(&(sysC.m), &(sysC.n), sysC.mat, &mytol, sysCinv.mat);
		//   Multiply x = Cinv * ym
		multmat(&sysCinv, &m_ym, &x);
		//   Init xd = zeros
		zeros(&xd);
//		fprintf (pFile, "Initial Conditions successfully %d \n",3);
//		fflush (pFile);
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
//		fprintf (pFile, "Step 1 successfully %d \n",5);
//		fflush (pFile);
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
//		fprintf (pFile, "Step 2 successfully %d \n",6);
//		fflush (pFile);
		// % Initialization Finished
		cStatus = 1;

	}
	else
	{
//		myCount++;
//		fprintf (pFile, "RUN %d \n",0);
//		fflush (pFile);
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
//		fprintf (pFile, "Step 1 successfully %6.6f \n",m_y.mat[1]);
//		fflush (pFile);
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
//		fprintf (pFile, "Step 2 successfully %6.6f \n",m_y.mat[1]);
//		fprintf (pFile, "Finished Round %d \n",myCount);
//		fflush (pFile);
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


void pinv(integer *m, integer *n, doublereal *a, doublereal *mytol, doublereal *ainv)
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

void multmat( mat *A, mat *B, mat *C )
{
	char transa      = 'N';
    char transb      = 'N';
    doublereal alpha = 1.0;
    doublereal beta  = 0.0;

    dgemm_(&transa, &transb,
	  	   &A->m, &B->n, &A->n, &alpha, A->mat, &A->m,
		   B->mat, &A->n, &beta, C->mat, &A->m);
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
