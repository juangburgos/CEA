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

#define EXPORTING_DLL

// to avoid C++ incompatibilities (after all, it is expecting C)
extern "C" {               // make sure these are visible in directory lookup
#include "./lapack/f2c.h"           // our solution: copy them to /usr/include
#include "./lapack/clapack.h"
#include "./sqlite/sqlite3.h"
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

		int i;
		for(i = 0; i < sysA.m*sysA.n; i++)
		{
			u[i]  = sysA.mat[i];
		}
		for(i = 0; i < disA.m*disA.n; i++)
		{
			ym[i]  = disA.mat[i];
		}
		for(i = 0; i < sRoll.m*sRoll.n; i++)
		{
			y[i]  = sRoll.mat[i];
		}

		// % Initialization Finished
		cStatus = 1;
		// Close Database
		sqlite3_close(db);
	}
	else
	{

	}

}

void zeros( int m, int n, mat *matrix )
{
	matrix->mat = (doublereal*) calloc(m*n, sizeof(doublereal));
	for (int i = 0; i < (m)*(n); i++)
	    {
			matrix->mat[i] = 0.0;
	    }
	matrix->m = m;
	matrix->n = n;

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
