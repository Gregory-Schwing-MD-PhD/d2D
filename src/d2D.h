#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#ifdef LIB
#include "globals.h"
#endif

#ifndef M_PI
#define M_PI      3.14159265358979323846
#endif

#ifndef LIB
#define MAXLENGTH 270000
#endif

#define NSIGMA    5.0

enum {
  A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y, X, Z, NR
};

enum {
  CA, CB, CO, HN, HA, NH, NSHIFT
};

#ifndef LIB
int debug=0;
#endif
int shifty=0;
int dbformat=0;

char* shtostring(int ashift)
{
  switch(ashift) {
    case CA: return "CA"; break;
    case CB: return "CB"; break;
    case CO: return "CO"; break;
    case HN: return "HN"; break;
    case HA: return "HA"; break;
    case NH: return "NH"; break;
  }
  return NULL;
}

double v_dot_p3(const double x[3], const double y[3])
{
    double res=0.0; int i;
    for(i=0;i<3;i++) res += x[i]*y[i];
    return res;
}

void m_v_m3(const double mat[3][3], const double vec[3], double *result)
{ 
    int i;
    for (i=0;i<3;i++) result[i] = v_dot_p3(mat[i], vec);
}

double v_dot_p4(const double x[4], const double y[4])
{
    double res=0.0; int i;
    for (i=0;i<4;i++) res += x[i]*y[i];
    return res;
}

void m_v_m4(const double mat[4][4], const double vec[4], double *result)
{ 
    int i;
    for (i=0;i<4;i++) result[i] = v_dot_p4(mat[i], vec);
}

double v_dot_p5(const double x[5], const double y[5])
{
    double res=0.0; int i;
    for (i=0;i<5;i++) res += x[i]*y[i];
    return res;
}

void m_v_m5(const double mat[5][5], const double vec[5], double *result)
{ 
    int i;
    for (i=0;i<5;i++) result[i] = v_dot_p5(mat[i], vec);
}

double v_dot_p6(const double x[6], const double y[6])
{
    double res=0.0; int i;
    for (i=0;i<6;i++) res += x[i]*y[i];
    return res;
}

void m_v_m6(const double mat[6][6], const double vec[6], double *result)
{ 
    int i;
    for (i=0;i<6;i++) result[i] = v_dot_p6(mat[i], vec);
}

void m_inv3(double mat[3][3], double id[3][3]) 
{
	const double Epsilon=0.0000001;
	int   i,j,s=0,pindex,pivot=1;                    
	double f,Maximum,h;

	for(i=0;i<3;i++) for(j=0;j<3;j++) {
          if(i==j) id[i][j] = 1.0;
          else id[i][j] = 0.;
        }
	
	do {
	  Maximum = fabs(mat[s][s]);
          /* find the highest in the coulumn */
	  if(pivot) {
	    pindex = s ; 
	    for(i=s+1;i<3;i++) 
              if(fabs(mat[i][s]) > Maximum) {
		Maximum = fabs(mat[i][s]) ;
		pindex = i;
	    }
	  }
	  if(Maximum < Epsilon) break;
          if(pivot) {
	    if(pindex!=s) { 
              for(j=s;j<3;j++) {
		h = mat[s][j];
		mat[s][j] = mat[pindex][j];
		mat[pindex][j]= h;
	      }
	      for(j=0;j<3;j++) {
		h = id[s][j];
		id[s][j] = id[pindex][j];
		id[pindex][j]= h;						
	      }
	    }
	  }
	  f = mat[s][s];
	  for(j=s;j<3;j++) mat[s][j] = mat[s][j]/f;
	  for(j=0;j<3;j++) id[s][j] = id[s][j]/f;
	  for(i=0;i<3;i++) {
	    if(i!=s) {
	      f = -mat[i][s];                 
	      for(j=s;j<3;j++) mat[i][j] += f*mat[s][j];
	      for(j=0;j<3;j++) id[i][j] += f*id[s][j];
	    }
	  }
	  s++;
	} while (s<3);
}

void m_inv4(double mat[4][4], double id[4][4]) 
{
	const double Epsilon=0.0000001;
	int   i,j,s=0,pindex,pivot=1;                    
	double f,Maximum,h;

	for(i=0;i<4;i++) for(j=0;j<4;j++) {
          if(i==j) id[i][j] = 1.0;
          else id[i][j] = 0.;
        }
	
	do {
	  Maximum = fabs(mat[s][s]);
          /* find the highest in the coulumn */
	  if(pivot) {
	    pindex = s ; 
	    for(i=s+1;i<4;i++) 
              if(fabs(mat[i][s]) > Maximum) {
		Maximum = fabs(mat[i][s]) ;
		pindex = i;
	    }
	  }
	  if(Maximum < Epsilon) break;
          if(pivot) {
	    if(pindex!=s) { 
              for(j=s;j<4;j++) {
		h = mat[s][j];
		mat[s][j] = mat[pindex][j];
		mat[pindex][j]= h;
	      }
	      for(j=0;j<4;j++) {
		h = id[s][j];
		id[s][j] = id[pindex][j];
		id[pindex][j]= h;						
	      }
	    }
	  }
	  f = mat[s][s];
	  for(j=s;j<4;j++) mat[s][j] = mat[s][j]/f;
	  for(j=0;j<4;j++) id[s][j] = id[s][j]/f;
	  for(i=0;i<4;i++) {
	    if(i!=s) {
	      f = -mat[i][s];                 
	      for(j=s;j<4;j++) mat[i][j] += f*mat[s][j];
	      for(j=0;j<4;j++) id[i][j] += f*id[s][j];
	    }
	  }
	  s++;
	} while (s<4);
}

void m_inv5(double mat[5][5], double id[5][5]) 
{
	const double Epsilon=0.0000001;
	int   i,j,s=0,pindex,pivot=1;                    
	double f,Maximum,h;

	for(i=0;i<5;i++) for(j=0;j<5;j++) {
          if(i==j) id[i][j] = 1.0;
          else id[i][j] = 0.;
        }
	
	do {
	  Maximum = fabs(mat[s][s]);
          /* find the highest in the coulumn */
	  if(pivot) {
	    pindex = s ; 
	    for(i=s+1;i<5;i++) 
              if(fabs(mat[i][s]) > Maximum) {
		Maximum = fabs(mat[i][s]) ;
		pindex = i;
	    }
	  }
	  if(Maximum < Epsilon) break;
          if(pivot) {
	    if(pindex!=s) { 
              for(j=s;j<5;j++) {
		h = mat[s][j];
		mat[s][j] = mat[pindex][j];
		mat[pindex][j]= h;
	      }
	      for(j=0;j<5;j++) {
		h = id[s][j];
		id[s][j] = id[pindex][j];
		id[pindex][j]= h;						
	      }
	    }
	  }
	  f = mat[s][s];
	  for(j=s;j<5;j++) mat[s][j] = mat[s][j]/f;
	  for(j=0;j<5;j++) id[s][j] = id[s][j]/f;
	  for(i=0;i<5;i++) {
	    if(i!=s) {
	      f = -mat[i][s];                 
	      for(j=s;j<5;j++) mat[i][j] += f*mat[s][j];
	      for(j=0;j<5;j++) id[i][j] += f*id[s][j];
	    }
	  }
	  s++;
	} while (s<5);
}

void m_inv6(double mat[6][6], double id[6][6]) 
{
	const double Epsilon=0.0000001;
	int   i,j,s=0,pindex,pivot=1;                    
	double f,Maximum,h;

	for(i=0;i<6;i++) for(j=0;j<6;j++) {
          if(i==j) id[i][j] = 1.0;
          else id[i][j] = 0.;
        }
	
	do {
	  Maximum = fabs(mat[s][s]);
          /* find the highest in the coulumn */
	  if(pivot) {
	    pindex = s ; 
	    for(i=s+1;i<6;i++) 
              if(fabs(mat[i][s]) > Maximum) {
		Maximum = fabs(mat[i][s]) ;
		pindex = i;
	    }
	  }
	  if(Maximum < Epsilon) break;
          if(pivot) {
	    if(pindex!=s) { 
              for(j=s;j<6;j++) {
		h = mat[s][j];
		mat[s][j] = mat[pindex][j];
		mat[pindex][j]= h;
	      }
	      for(j=0;j<6;j++) {
		h = id[s][j];
		id[s][j] = id[pindex][j];
		id[pindex][j]= h;						
	      }
	    }
	  }
	  f = mat[s][s];
	  for(j=s;j<6;j++) mat[s][j] = mat[s][j]/f;
	  for(j=0;j<6;j++) id[s][j] = id[s][j]/f;
	  for(i=0;i<6;i++) {
	    if(i!=s) {
	      f = -mat[i][s];                 
	      for(j=s;j<6;j++) mat[i][j] += f*mat[s][j];
	      for(j=0;j<6;j++) id[i][j] += f*id[s][j];
	    }
	  }
	  s++;
	} while (s<6);
}

int cmp(const void *x, const void *y)
{
  double xx = *(double*)x, yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
}
