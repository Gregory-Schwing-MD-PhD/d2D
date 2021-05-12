/* d2D - v. 1.4.0
 * (c) Carlo Camilloni
 * Camilloni C., De Simone A., Vranken W., and Vendruscolo M.
 * Determination of Secondary Structure Populations in Disordered States of Proteins using NMR Chemical Shifts
 * Biochemistry 2012, 51: 2224-2231 
 */

/* 
 * compile with
 * gcc  predictor.c -lm -O3 -Wall -Wextra -ansi -pedantic -o d2D.x -I. -DLORENTZ 
 */
 
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
#define M_PI            3.14159265358979323846
#endif

#ifdef LIB
#define MAXLENGTH 12800
#else
#define MAXLENGTH 220000
#endif

enum {
  A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y, X, Z, NR
};

enum {
  CA, CB, CO, HN, HA, NH, NSHIFT
};

#ifndef LIB
int debug=0;
#endif
int ppii=1;
int shifty=0;
int force=0;

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

void do_averaging(double *pp[], double *spp[], int *nseq, int sql, double wc[6][6], int *star, char ss[MAXLENGTH])
{
  int    i, j, MAX;
  double dcount, maxi, norm;

  for(i=0;i<sql;i++)
  {
    if(pp[i][0]!=0.||pp[i][1]!=0.||pp[i][2]!=0.) 
    {
      if(star[i]<3) {
        if(((i>0)&&(i<sql-1))&&star[i-1]>2&&star[i+1]>2&&
           (pp[i-1][0]!=0.||pp[i-1][1]!=0.||pp[i-1][2]!=0.)&&
           (pp[i+1][0]!=0.||pp[i+1][1]!=0.||pp[i+1][2]!=0.)){
          spp[i][0] = (pp[i+1][0]+pp[i-1][0])/2.;
          spp[i][1] = (pp[i+1][1]+pp[i-1][1])/2.;
          spp[i][2] = (pp[i+1][2]+pp[i-1][2])/2.;
          spp[i][3] = (pp[i+1][3]+pp[i-1][3])/2.;
          dcount=1.0;
        } else {
          spp[i][0] = 0.;
          spp[i][1] = 0.;
          spp[i][2] = 0.;
          spp[i][3] = 0.;
          ss[i]=' ';
          continue;
        }
      } else {
        spp[i][0] = pp[i][0];
        spp[i][1] = pp[i][1];
        spp[i][2] = pp[i][2];
        spp[i][3] = pp[i][3];
        dcount=1.0;

        if(i>2) if(nseq[i-1]!=Z&&nseq[i-2]!=Z&&(pp[i-2][0]!=0.||pp[i-2][1]!=0.||pp[i-2][2]!=0.)&&star[i-2]>2) {
          spp[i][0] += wc[5][0]*pp[i-2][0];
          spp[i][1] += wc[5][0]*pp[i-2][1];
          spp[i][2] += wc[5][0]*pp[i-2][2];
          spp[i][3] += wc[5][0]*pp[i-2][3];
          dcount+=wc[5][0];
        }

        if(i>1) if(nseq[i-1]!=Z&&(pp[i-1][0]!=0.||pp[i-1][1]!=0.||pp[i-1][2]!=0.)&&star[i-1]>2) {
          spp[i][0] += wc[5][1]*pp[i-1][0];
          spp[i][1] += wc[5][1]*pp[i-1][1];
          spp[i][2] += wc[5][1]*pp[i-1][2];
          spp[i][3] += wc[5][1]*pp[i-1][3];
          dcount+=wc[5][1];
        }
        if(i<sql-2) if(nseq[i+1]!=Z&&(pp[i+1][0]!=0.||pp[i+1][1]!=0.||pp[i+1][2]!=0.)&&star[i+1]>2) {
          spp[i][0] += wc[5][1]*pp[i+1][0];
          spp[i][1] += wc[5][1]*pp[i+1][1];
          spp[i][2] += wc[5][1]*pp[i+1][2];
          spp[i][3] += wc[5][1]*pp[i+1][3];
          dcount+=wc[5][1];
        }

        if(i<sql-3) if(nseq[i+1]!=Z&&nseq[i+2]!=Z&&(pp[i+2][0]!=0.||pp[i+2][1]!=0.||pp[i+2][2]!=0.)&&star[i+2]>2) {
          spp[i][0] += wc[5][0]*pp[i+2][0];
          spp[i][1] += wc[5][0]*pp[i+2][1];
          spp[i][2] += wc[5][0]*pp[i+2][2];
          spp[i][3] += wc[5][0]*pp[i+2][3];
          dcount+=wc[5][0];
        }
      }
      spp[i][0] /= dcount;
      spp[i][1] /= dcount;
      spp[i][2] /= dcount;
      spp[i][3] /= dcount;

      if(ppii) {
        norm = spp[i][0] + spp[i][1] + spp[i][2] + spp[i][3];
        spp[i][3] /= norm;
      } else norm = spp[i][0] + spp[i][1] + spp[i][2];
      spp[i][0] /= norm;
      spp[i][1] /= norm;
      spp[i][2] /= norm;

      maxi = spp[i][0]; MAX=0;
      if(ppii) {
        for(j=1;j<4;j++) if( maxi < spp[i][j] ) {maxi = spp[i][j]; MAX=j;}
      } else {
        for(j=1;j<3;j++) if( maxi < spp[i][j] ) {maxi = spp[i][j]; MAX=j;}
      }

      switch(MAX)
      {
        case 0: ss[i]='H'; break;
        case 1: ss[i]='E'; break;
        case 2: ss[i]='C'; break;
        case 3: ss[i]='P'; break;
      }
    } else ss[i]=' ';
  }
} 

int do_alpha(int ashift, double *tcs, int *nseq, int sql, int ph, char *dbpath)
{
  char   n1[25], n2[25], n3[25], n4[25], n5[25], n6[25], *format, str[MAXLENGTH];
  double mean[21];
  double s1[21][21], d3[21][21], d1[21][21], d2[21][21], d4[21][21];
  double as1, ad1, ad2, ad3, ad4;
  int    i, err;
  FILE   *m=NULL, *uno, *due, *tre, *qua, *cin, *par;

  switch(ashift)
  {
    case CA: strcpy(n1,"/helix-db/CA-uno.dat"); 
             strcpy(n2,"/helix-db/CA-due.dat"); 
             strcpy(n3,"/helix-db/CA-tre.dat"); 
             strcpy(n4,"/helix-db/CA-quattro.dat"); 
             strcpy(n6,"/helix-db/CA-parm.dat"); 
             strcpy(n5,"/helix-db/CA-cin.dat"); 
             break;
    case CB: strcpy(n1,"/helix-db/CB-uno.dat"); 
             strcpy(n2,"/helix-db/CB-due.dat"); 
             strcpy(n3,"/helix-db/CB-tre.dat"); 
             strcpy(n4,"/helix-db/CB-quattro.dat"); 
             strcpy(n6,"/helix-db/CB-parm.dat"); 
             strcpy(n5,"/helix-db/CB-cin.dat"); 
             break;
    case CO: strcpy(n1,"/helix-db/CO-uno.dat"); 
             strcpy(n2,"/helix-db/CO-due.dat"); 
             strcpy(n3,"/helix-db/CO-tre.dat"); 
             strcpy(n4,"/helix-db/CO-quattro.dat"); 
             strcpy(n6,"/helix-db/CO-parm.dat");  
             strcpy(n5,"/helix-db/CO-cin.dat"); 
             break;
    case HN: strcpy(n1,"/helix-db/HN-uno.dat"); 
             strcpy(n2,"/helix-db/HN-due.dat"); 
             strcpy(n3,"/helix-db/HN-tre.dat"); 
             strcpy(n4,"/helix-db/HN-quattro.dat"); 
             strcpy(n6,"/helix-db/HN-parm.dat"); 
             strcpy(n5,"/helix-db/HN-cin.dat"); 
             break;
    case HA: strcpy(n1,"/helix-db/HA-uno.dat"); 
             strcpy(n2,"/helix-db/HA-due.dat"); 
             strcpy(n3,"/helix-db/HA-tre.dat"); 
             strcpy(n4,"/helix-db/HA-quattro.dat"); 
             strcpy(n6,"/helix-db/HA-parm.dat"); 
             strcpy(n5,"/helix-db/HA-cin.dat"); 
             break;
    case NH: strcpy(n1,"/helix-db/NH-uno.dat"); 
             strcpy(n2,"/helix-db/NH-due.dat"); 
             strcpy(n3,"/helix-db/NH-tre.dat"); 
             strcpy(n4,"/helix-db/NH-quattro.dat"); 
             strcpy(n6,"/helix-db/NH-parm.dat"); 
             strcpy(n5,"/helix-db/NH-cin.dat"); 
             break;
  }
  sprintf(str, "%s%s", dbpath, n1);
  uno = fopen(str,"r");
  if (uno==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  sprintf(str, "%s%s", dbpath, n2);
  due = fopen(str,"r");
  if (due==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  sprintf(str, "%s%s", dbpath, n3);
  tre = fopen(str,"r");
  if (tre==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  sprintf(str, "%s%s", dbpath, n4);
  qua = fopen(str,"r");
  if (qua==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  sprintf(str, "%s%s", dbpath, n5);
  cin = fopen(str,"r");
  if (cin==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  if(ph==0) {
    switch(ashift)
    {
      case CA: sprintf(str,"%s/helix-db/cs-ca-medi.dat",dbpath); break;
      case CB: sprintf(str,"%s/helix-db/cs-cb-medi.dat",dbpath); break;
      case CO: sprintf(str,"%s/helix-db/cs-co-medi.dat",dbpath); break;
      case HN: sprintf(str,"%s/helix-db/cs-hn-medi.dat",dbpath); break;
      case HA: sprintf(str,"%s/helix-db/cs-ha-medi.dat",dbpath); break;
      case NH: sprintf(str,"%s/helix-db/cs-n-medi.dat",dbpath); break;
    }
    m = fopen(str, "r");
    if (m==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}
  } else {
    switch(ashift)
    {
      case CA: sprintf(str,"%s/ph3-helix/cs-ca-medi.dat",dbpath); break;
      case CB: sprintf(str,"%s/ph3-helix/cs-cb-medi.dat",dbpath); break;
      case CO: sprintf(str,"%s/ph3-helix/cs-co-medi.dat",dbpath); break;
      case HN: sprintf(str,"%s/ph3-helix/cs-hn-medi.dat",dbpath); break;
      case HA: sprintf(str,"%s/ph3-helix/cs-ha-medi.dat",dbpath); break;
      case NH: sprintf(str,"%s/ph3-helix/cs-n-medi.dat",dbpath);  break;
    }
    m = fopen(str, "r");
    if (m==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}
  }

  sprintf(str, "%s%s", dbpath, n6);
  par = fopen(str,"r");
  if (par==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  err = fscanf(uno, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  err = fscanf(due, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  err = fscanf(tre, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  err = fscanf(qua, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  err = fscanf(cin, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
 
  format = "%*s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf";
  for(i=0;i<21;i++)
  {
    err = fscanf(m, "%*s %lf %*f %*f %*f", &mean[i]);

    err = fscanf(uno, format, &s1[i][0], &s1[i][1], &s1[i][2], &s1[i][3], &s1[i][4], &s1[i][5], &s1[i][6], &s1[i][7], &s1[i][8], &s1[i][9], &s1[i][10], &s1[i][11], &s1[i][12], &s1[i][13], &s1[i][14], &s1[i][15], &s1[i][16], &s1[i][17], &s1[i][18], &s1[i][19], &s1[i][20]);

    err = fscanf(due, format, &d1[i][0], &d1[i][1], &d1[i][2], &d1[i][3], &d1[i][4], &d1[i][5], &d1[i][6], &d1[i][7], &d1[i][8], &d1[i][9], &d1[i][10], &d1[i][11], &d1[i][12], &d1[i][13], &d1[i][14], &d1[i][15], &d1[i][16], &d1[i][17], &d1[i][18], &d1[i][19], &d1[i][20]);

    err = fscanf(tre, format, &d2[i][0], &d2[i][1], &d2[i][2], &d2[i][3], &d2[i][4], &d2[i][5], &d2[i][6], &d2[i][7], &d2[i][8], &d2[i][9], &d2[i][10], &d2[i][11], &d2[i][12], &d2[i][13], &d2[i][14], &d2[i][15], &d2[i][16], &d2[i][17], &d2[i][18], &d2[i][19], &d2[i][20]);

    err = fscanf(qua, format, &d3[i][0], &d3[i][1], &d3[i][2], &d3[i][3], &d3[i][4], &d3[i][5], &d3[i][6], &d3[i][7], &d3[i][8], &d3[i][9], &d3[i][10], &d3[i][11], &d3[i][12], &d3[i][13], &d3[i][14], &d3[i][15], &d3[i][16], &d3[i][17], &d3[i][18], &d3[i][19], &d3[i][20]);

    err = fscanf(cin, format, &d4[i][0], &d4[i][1], &d4[i][2], &d4[i][3], &d4[i][4], &d4[i][5], &d4[i][6], &d4[i][7], &d4[i][8], &d4[i][9], &d4[i][10], &d4[i][11], &d4[i][12], &d4[i][13], &d4[i][14], &d4[i][15], &d4[i][16], &d4[i][17], &d4[i][18], &d4[i][19], &d4[i][20]);
  }

  fclose(m);
  fclose(uno);
  fclose(due);
  fclose(tre);
  fclose(qua);
  fclose(cin);

  err = fscanf(par, "%lf", &as1);
  err = fscanf(par, "%lf", &ad1);
  err = fscanf(par, "%lf", &ad2);
  err = fscanf(par, "%lf", &ad3);
  err = fscanf(par, "%lf", &ad4);
  
  fclose(par);

  for(i=0;i<sql;i++)
  {
    if(nseq[i]==G&&(ashift==HA||ashift==CB)) continue;
    if(nseq[i]==P&&(ashift==HN||ashift==NH)) continue;
    if(nseq[i]==Z) continue;

    tcs[i] = mean[nseq[i]];

    if(i>0&&nseq[i-1]!=Z)                                                 tcs[i] += as1*s1[nseq[i-1]][nseq[i]]; 
    if(i<(sql-1)&&nseq[i+1]!=Z)                                           tcs[i] += ad1*d1[nseq[i+1]][nseq[i]]; 
    if(i<(sql-2)&&nseq[i+1]!=Z&&nseq[i+2]!=Z)                             tcs[i] += ad2*d2[nseq[i+2]][nseq[i]]; 
    if(i<(sql-3)&&nseq[i+1]!=Z&&nseq[i+2]!=Z&&nseq[i+3]!=Z)               tcs[i] += ad3*d3[nseq[i+3]][nseq[i]]; 
    if(i<(sql-4)&&nseq[i+1]!=Z&&nseq[i+2]!=Z&&nseq[i+3]!=Z&&nseq[i+4]!=Z) tcs[i] += ad4*d4[nseq[i+4]][nseq[i]];

    if(debug) {fprintf(stdout, "A %i_CS: %i %f \n", ashift, i, tcs[i]); fflush(stdout);}
  }

  return 0;
}

int do_beta(int ashift, double *tcs, int *nseq, int sql, int ph, char *dbpath)
{
  char   n1[25], n2[25], n3[25], n4[25], n5[25], *format, str[MAXLENGTH];
  double mean[21];
  double s1[21][21], s2[21][21], d1[21][21], d2[21][21];
  double as1, as2, ad1, ad2;
  int    i, err;
  FILE   *m=NULL, *uno, *due, *tre, *qua, *par;

  switch(ashift)
  {
    case CA: strcpy(n1,"beta-db/CA-uno.dat"); strcpy(n2,"beta-db/CA-due.dat"); strcpy(n3,"beta-db/CA-tre.dat"); strcpy(n4,"beta-db/CA-quattro.dat"); strcpy(n5,"beta-db/CA-parm.dat"); break;
    case CB: strcpy(n1,"beta-db/CB-uno.dat"); strcpy(n2,"beta-db/CB-due.dat"); strcpy(n3,"beta-db/CB-tre.dat"); strcpy(n4,"beta-db/CB-quattro.dat"); strcpy(n5,"beta-db/CB-parm.dat"); break;
    case CO: strcpy(n1,"beta-db/CO-uno.dat"); strcpy(n2,"beta-db/CO-due.dat"); strcpy(n3,"beta-db/CO-tre.dat"); strcpy(n4,"beta-db/CO-quattro.dat"); strcpy(n5,"beta-db/CO-parm.dat"); break;
    case HN: strcpy(n1,"beta-db/HN-uno.dat"); strcpy(n2,"beta-db/HN-due.dat"); strcpy(n3,"beta-db/HN-tre.dat"); strcpy(n4,"beta-db/HN-quattro.dat"); strcpy(n5,"beta-db/HN-parm.dat"); break;
    case HA: strcpy(n1,"beta-db/HA-uno.dat"); strcpy(n2,"beta-db/HA-due.dat"); strcpy(n3,"beta-db/HA-tre.dat"); strcpy(n4,"beta-db/HA-quattro.dat"); strcpy(n5,"beta-db/HA-parm.dat"); break;
    case NH: strcpy(n1,"beta-db/NH-uno.dat"); strcpy(n2,"beta-db/NH-due.dat"); strcpy(n3,"beta-db/NH-tre.dat"); strcpy(n4,"beta-db/NH-quattro.dat"); strcpy(n5,"beta-db/NH-parm.dat"); break;
  }
  sprintf(str, "%s/%s", dbpath, n1);
  uno = fopen(str,"r");
  if (uno==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  sprintf(str, "%s/%s", dbpath, n2);
  due = fopen(str,"r");
  if (due==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  sprintf(str, "%s/%s", dbpath, n3);
  tre = fopen(str,"r");
  if (tre==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  sprintf(str, "%s/%s", dbpath, n4);
  qua = fopen(str,"r");
  if (qua==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  if(ph==0) {
    switch(ashift)
    {
      case CA: sprintf(str,"%s/beta-db/cs-ca-medi.dat", dbpath); break;
      case CB: sprintf(str,"%s/beta-db/cs-cb-medi.dat", dbpath); break;
      case CO: sprintf(str,"%s/beta-db/cs-co-medi.dat", dbpath); break;
      case HN: sprintf(str,"%s/beta-db/cs-hn-medi.dat", dbpath); break;
      case HA: sprintf(str,"%s/beta-db/cs-ha-medi.dat", dbpath); break;
      case NH: sprintf(str,"%s/beta-db/cs-n-medi.dat", dbpath); break;
    }
    m = fopen(str, "r");
    if (m==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}
  } else {
    switch(ashift)
    {
      case CA: sprintf(str,"%s/ph3-beta/cs-ca-medi.dat", dbpath); break;
      case CB: sprintf(str,"%s/ph3-beta/cs-cb-medi.dat", dbpath); break;
      case CO: sprintf(str,"%s/ph3-beta/cs-co-medi.dat", dbpath); break;
      case HN: sprintf(str,"%s/ph3-beta/cs-hn-medi.dat", dbpath); break;
      case HA: sprintf(str,"%s/ph3-beta/cs-ha-medi.dat", dbpath); break;
      case NH: sprintf(str,"%s/ph3-beta/cs-n-medi.dat", dbpath);  break;
    }
    m = fopen(str, "r");
    if (m==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}
  }

  sprintf(str, "%s/%s", dbpath, n5);
  par = fopen(str,"r");
  if (par==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  err = fscanf(uno, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  err = fscanf(due, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  err = fscanf(tre, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  err = fscanf(qua, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
 
  format = "%*s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf";
  for(i=0;i<21;i++)
  {
    err = fscanf(m, "%*s %lf %*f", &mean[i]);

    err = fscanf(uno, format, &s1[i][0], &s1[i][1], &s1[i][2], &s1[i][3], &s1[i][4], &s1[i][5], &s1[i][6], &s1[i][7], &s1[i][8], &s1[i][9], &s1[i][10], &s1[i][11], &s1[i][12], &s1[i][13], &s1[i][14], &s1[i][15], &s1[i][16], &s1[i][17], &s1[i][18], &s1[i][19], &s1[i][20]);

    err = fscanf(due, format, &s2[i][0], &s2[i][1], &s2[i][2], &s2[i][3], &s2[i][4], &s2[i][5], &s2[i][6], &s2[i][7], &s2[i][8], &s2[i][9], &s2[i][10], &s2[i][11], &s2[i][12], &s2[i][13], &s2[i][14], &s2[i][15], &s2[i][16], &s2[i][17], &s2[i][18], &s2[i][19], &s2[i][20]);

    err = fscanf(tre, format, &d1[i][0], &d1[i][1], &d1[i][2], &d1[i][3], &d1[i][4], &d1[i][5], &d1[i][6], &d1[i][7], &d1[i][8], &d1[i][9], &d1[i][10], &d1[i][11], &d1[i][12], &d1[i][13], &d1[i][14], &d1[i][15], &d1[i][16], &d1[i][17], &d1[i][18], &d1[i][19], &d1[i][20]);

    err = fscanf(qua, format, &d2[i][0], &d2[i][1], &d2[i][2], &d2[i][3], &d2[i][4], &d2[i][5], &d2[i][6], &d2[i][7], &d2[i][8], &d2[i][9], &d2[i][10], &d2[i][11], &d2[i][12], &d2[i][13], &d2[i][14], &d2[i][15], &d2[i][16], &d2[i][17], &d2[i][18], &d2[i][19], &d2[i][20]);
  }

  fclose(m);
  fclose(uno);
  fclose(due);
  fclose(tre);
  fclose(qua);

  err = fscanf(par, "%lf", &as1);
  err = fscanf(par, "%lf", &as2);
  err = fscanf(par, "%lf", &ad1);
  err = fscanf(par, "%lf", &ad2);
  
  fclose(par);

  for(i=0;i<sql;i++)
  {
    if(nseq[i]==G&&(ashift==HA||ashift==CB)) continue;
    if(nseq[i]==P&&(ashift==HN||ashift==NH)) continue;
    if(nseq[i]==Z) continue;

    tcs[i] = mean[nseq[i]]; 

    if(i>1&&nseq[i-2]!=Z&&nseq[i-1]!=Z)       tcs[i] += as2*s2[nseq[i-2]][nseq[i]];
    if(i>0&&nseq[i-1]!=Z)                     tcs[i] += as1*s1[nseq[i-1]][nseq[i]];
    if(i<(sql-1)&&nseq[i+1]!=Z)               tcs[i] += ad1*d1[nseq[i+1]][nseq[i]];
    if(i<(sql-2)&&nseq[i+1]!=Z&&nseq[i+2]!=Z) tcs[i] += ad2*d2[nseq[i+2]][nseq[i]];

    if(debug) {fprintf(stdout, "B %i_CS: %i %f\n", ashift, i, tcs[i]); fflush(stdout);}
  }

  return 0;
}

int do_tamiola(int ashift, double *tcs, int *nseq, int sql, int ph, char *dbpath)
{
  char str[MAXLENGTH];
  double mean[21];
  double sx[21], dx[21];
  int  i, err;
  FILE *m=NULL;

  if(ph==0) {
    switch(ashift)
    {
      case CA: sprintf(str,"%s/tamiola-db/sCAdata.dat",dbpath); break;
      case CB: sprintf(str,"%s/tamiola-db/sCBdata.dat",dbpath); break;
      case CO: sprintf(str,"%s/tamiola-db/sCOdata.dat",dbpath); break;
      case HN: sprintf(str,"%s/tamiola-db/sHNdata.dat",dbpath); break;
      case HA: sprintf(str,"%s/tamiola-db/sHAdata.dat",dbpath); break;
      case NH: sprintf(str,"%s/tamiola-db/sNHdata.dat",dbpath); break;
    }
  }
  m = fopen(str,"r");
  if (m==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}
  for(i=0;i<21;i++)
    err = fscanf(m, "%*s %lf %lf %lf", &sx[i], &mean[i], &dx[i]);

  fclose(m);

  for(i=0;i<sql;i++)
  {
    if(nseq[i]==G&&(ashift==HA||ashift==CB)) continue;
    if(nseq[i]==P&&(ashift==HN||ashift==NH)) continue;
    if(nseq[i]==Z) continue;
    else tcs[i] = mean[nseq[i]]; 

    if(i>0) if(nseq[i-1]!=Z) { tcs[i] += sx[nseq[i-1]]; }
    if(i<sql-1) if(nseq[i+1]!=Z) {tcs[i] += dx[nseq[i+1]]; }

    if(debug) {fprintf(stdout, "C %i_CS: %i %f\n", ashift, i, tcs[i]); fflush(stdout);}
  }

  return 0;
}

int do_coil(int ashift, double *tcs, int *nseq, int sql, int ph, char *dbpath)
{
  char   n1[25], n2[25], n3[25], n4[25], n5[25], str[MAXLENGTH];
  char   *format;
  double mean[21];
  double s1[21][21], s2[21][21], d1[21][21], d2[21][21];
  double as1, as2, ad1, ad2;
  int    i, err;
  FILE   *m=NULL, *uno, *due, *tre, *qua, *par;

  switch(ashift)
  {
    case CA: strcpy(n1,"coil-db/CA-uno.dat"); strcpy(n2,"coil-db/CA-due.dat"); strcpy(n3,"coil-db/CA-tre.dat"); strcpy(n4,"coil-db/CA-quattro.dat"); strcpy(n5,"coil-db/CA-parm.dat"); break;
    case CB: strcpy(n1,"coil-db/CB-uno.dat"); strcpy(n2,"coil-db/CB-due.dat"); strcpy(n3,"coil-db/CB-tre.dat"); strcpy(n4,"coil-db/CB-quattro.dat"); strcpy(n5,"coil-db/CB-parm.dat"); break;
    case CO: strcpy(n1,"coil-db/CO-uno.dat"); strcpy(n2,"coil-db/CO-due.dat"); strcpy(n3,"coil-db/CO-tre.dat"); strcpy(n4,"coil-db/CO-quattro.dat"); strcpy(n5,"coil-db/CO-parm.dat"); break;
    case HN: strcpy(n1,"coil-db/HN-uno.dat"); strcpy(n2,"coil-db/HN-due.dat"); strcpy(n3,"coil-db/HN-tre.dat"); strcpy(n4,"coil-db/HN-quattro.dat"); strcpy(n5,"coil-db/HN-parm.dat"); break;
    case HA: strcpy(n1,"coil-db/HA-uno.dat"); strcpy(n2,"coil-db/HA-due.dat"); strcpy(n3,"coil-db/HA-tre.dat"); strcpy(n4,"coil-db/HA-quattro.dat"); strcpy(n5,"coil-db/HA-parm.dat"); break;
    case NH: strcpy(n1,"coil-db/NH-uno.dat"); strcpy(n2,"coil-db/NH-due.dat"); strcpy(n3,"coil-db/NH-tre.dat"); strcpy(n4,"coil-db/NH-quattro.dat"); strcpy(n5,"coil-db/NH-parm.dat"); break;
  }
  sprintf(str, "%s/%s", dbpath, n1);
  uno = fopen(str,"r");
  if (uno==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  sprintf(str, "%s/%s", dbpath, n2);
  due = fopen(str,"r");
  if (due==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  sprintf(str, "%s/%s", dbpath, n3);
  tre = fopen(str,"r");
  if (tre==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  sprintf(str, "%s/%s", dbpath, n4);
  qua = fopen(str,"r");
  if (qua==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  if(ph==0) {
    switch(ashift)
    {
      case CA: sprintf(str, "%s/coil-db/cs-ca-medi.dat", dbpath); break;
      case CB: sprintf(str, "%s/coil-db/cs-cb-medi.dat", dbpath); break;
      case CO: sprintf(str, "%s/coil-db/cs-co-medi.dat", dbpath); break;
      case HN: sprintf(str, "%s/coil-db/cs-hn-medi.dat", dbpath); break;
      case HA: sprintf(str, "%s/coil-db/cs-ha-medi.dat", dbpath); break;
      case NH: sprintf(str, "%s/coil-db/cs-n-medi.dat", dbpath); break;
    }
    m = fopen(str, "r");
    if (m==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}
  } else {
    switch(ashift)
    {
      case CA: sprintf(str, "%s/ph3-coil/cs-ca-medi.dat", dbpath); break;
      case CB: sprintf(str, "%s/ph3-coil/cs-cb-medi.dat", dbpath); break;
      case CO: sprintf(str, "%s/ph3-coil/cs-co-medi.dat", dbpath); break;
      case HN: sprintf(str, "%s/ph3-coil/cs-hn-medi.dat", dbpath); break;
      case HA: sprintf(str, "%s/ph3-coil/cs-ha-medi.dat", dbpath); break;
      case NH: sprintf(str, "%s/ph3-coil/cs-n-medi.dat", dbpath);  break;
    }
    m = fopen(str, "r");
    if (m==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}
  }

  sprintf(str, "%s/%s", dbpath, n5);
  par = fopen(str,"r");
  if (par==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  err = fscanf(uno, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  err = fscanf(due, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  err = fscanf(tre, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  err = fscanf(qua, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
 
  format = "%*s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf";
  for(i=0;i<21;i++)
  {
    err = fscanf(m, "%*s %lf %*f", &mean[i]);

    err = fscanf(uno, format, &s1[i][0], &s1[i][1], &s1[i][2], &s1[i][3], &s1[i][4], &s1[i][5], &s1[i][6], &s1[i][7], &s1[i][8], &s1[i][9], &s1[i][10], &s1[i][11], &s1[i][12], &s1[i][13], &s1[i][14], &s1[i][15], &s1[i][16], &s1[i][17], &s1[i][18], &s1[i][19], &s1[i][20]);

    err = fscanf(due, format, &s2[i][0], &s2[i][1], &s2[i][2], &s2[i][3], &s2[i][4], &s2[i][5], &s2[i][6], &s2[i][7], &s2[i][8], &s2[i][9], &s2[i][10], &s2[i][11], &s2[i][12], &s2[i][13], &s2[i][14], &s2[i][15], &s2[i][16], &s2[i][17], &s2[i][18], &s2[i][19], &s2[i][20]);

    err = fscanf(tre, format, &d1[i][0], &d1[i][1], &d1[i][2], &d1[i][3], &d1[i][4], &d1[i][5], &d1[i][6], &d1[i][7], &d1[i][8], &d1[i][9], &d1[i][10], &d1[i][11], &d1[i][12], &d1[i][13], &d1[i][14], &d1[i][15], &d1[i][16], &d1[i][17], &d1[i][18], &d1[i][19], &d1[i][20]);

    err = fscanf(qua, format, &d2[i][0], &d2[i][1], &d2[i][2], &d2[i][3], &d2[i][4], &d2[i][5], &d2[i][6], &d2[i][7], &d2[i][8], &d2[i][9], &d2[i][10], &d2[i][11], &d2[i][12], &d2[i][13], &d2[i][14], &d2[i][15], &d2[i][16], &d2[i][17], &d2[i][18], &d2[i][19], &d2[i][20]);
  }

  fclose(m);
  fclose(uno);
  fclose(due);
  fclose(tre);
  fclose(qua);

  err = fscanf(par, "%lf", &as1);
  err = fscanf(par, "%lf", &as2);
  err = fscanf(par, "%lf", &ad1);
  err = fscanf(par, "%lf", &ad2);
  
  fclose(par);

  for(i=0;i<sql;i++)
  {
    if(nseq[i]==G&&(ashift==HA||ashift==CB)) continue;
    if(nseq[i]==P&&(ashift==HN||ashift==NH)) continue;
    if(nseq[i]==Z) continue;

    tcs[i] = mean[nseq[i]]; 

    if(i>1&&nseq[i-2]!=Z&&nseq[i-1]!=Z)       tcs[i] += as2*s2[nseq[i-2]][nseq[i]];
    if(i>0&&nseq[i-1]!=Z)                     tcs[i] += as1*s1[nseq[i-1]][nseq[i]];
    if(i<(sql-1)&&nseq[i+1]!=Z)               tcs[i] += ad1*d1[nseq[i+1]][nseq[i]];
    if(i<(sql-2)&&nseq[i+1]!=Z&&nseq[i+2]!=Z) tcs[i] += ad2*d2[nseq[i+2]][nseq[i]];

    if(debug) {fprintf(stdout, "C %i_CS: %i %f\n", ashift, i, tcs[i]); fflush(stdout);}
  }

  return 0;
}

int do_ppii(int ashift, double *tcs, int *nseq, int sql, int ph, char *dbpath)
{
  char   n1[25], n2[25], n3[25], n4[25], n5[25], *format, str[MAXLENGTH];
  double mean[21];
  double s1[21][21], s2[21][21], d1[21][21], d2[21][21];
  double as1, as2, ad1, ad2;
  int    i, err;
  FILE   *m=NULL, *uno, *due, *tre, *qua, *par;

  switch(ashift)
  {
    case CA: strcpy(n1,"ppii-db/CA-uno.dat"); strcpy(n2,"ppii-db/CA-due.dat"); strcpy(n3,"ppii-db/CA-tre.dat"); strcpy(n4,"ppii-db/CA-quattro.dat"); strcpy(n5,"ppii-db/CA-parm.dat"); break;
    case CB: strcpy(n1,"ppii-db/CB-uno.dat"); strcpy(n2,"ppii-db/CB-due.dat"); strcpy(n3,"ppii-db/CB-tre.dat"); strcpy(n4,"ppii-db/CB-quattro.dat"); strcpy(n5,"ppii-db/CB-parm.dat"); break;
    case CO: strcpy(n1,"ppii-db/CO-uno.dat"); strcpy(n2,"ppii-db/CO-due.dat"); strcpy(n3,"ppii-db/CO-tre.dat"); strcpy(n4,"ppii-db/CO-quattro.dat"); strcpy(n5,"ppii-db/CO-parm.dat"); break;
    case HN: strcpy(n1,"ppii-db/HN-uno.dat"); strcpy(n2,"ppii-db/HN-due.dat"); strcpy(n3,"ppii-db/HN-tre.dat"); strcpy(n4,"ppii-db/HN-quattro.dat"); strcpy(n5,"ppii-db/HN-parm.dat"); break;
    case HA: strcpy(n1,"ppii-db/HA-uno.dat"); strcpy(n2,"ppii-db/HA-due.dat"); strcpy(n3,"ppii-db/HA-tre.dat"); strcpy(n4,"ppii-db/HA-quattro.dat"); strcpy(n5,"ppii-db/HA-parm.dat"); break;
    case NH: strcpy(n1,"ppii-db/NH-uno.dat"); strcpy(n2,"ppii-db/NH-due.dat"); strcpy(n3,"ppii-db/NH-tre.dat"); strcpy(n4,"ppii-db/NH-quattro.dat"); strcpy(n5,"ppii-db/NH-parm.dat"); break;
  }
  sprintf(str, "%s/%s", dbpath, n1);
  uno = fopen(str,"r");
  if (uno==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  sprintf(str, "%s/%s", dbpath, n2);
  due = fopen(str,"r");
  if (due==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  sprintf(str, "%s/%s", dbpath, n3);
  tre = fopen(str,"r");
  if (tre==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  sprintf(str, "%s/%s", dbpath, n4);
  qua = fopen(str,"r");
  if (qua==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  if(ph==0) {
    switch(ashift)
    {
      case CA: sprintf(str, "%s/ppii-db/cs-ca-medi.dat", dbpath); break;
      case CB: sprintf(str, "%s/ppii-db/cs-cb-medi.dat",dbpath); break;
      case CO: sprintf(str, "%s/ppii-db/cs-co-medi.dat",dbpath); break;
      case HN: sprintf(str, "%s/ppii-db/cs-hn-medi.dat",dbpath); break;
      case HA: sprintf(str, "%s/ppii-db/cs-ha-medi.dat",dbpath); break;
      case NH: sprintf(str, "%s/ppii-db/cs-n-medi.dat",dbpath); break;
    }
    m = fopen(str, "r");
    if (m==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}
  } else {
    switch(ashift)
    {
      case CA: sprintf(str, "%s/ph3-ppii/cs-ca-medi.dat", dbpath); break;
      case CB: sprintf(str, "%s/ph3-ppii/cs-cb-medi.dat", dbpath); break;
      case CO: sprintf(str, "%s/ph3-ppii/cs-co-medi.dat", dbpath); break;
      case HN: sprintf(str, "%s/ph3-ppii/cs-hn-medi.dat", dbpath); break;
      case HA: sprintf(str, "%s/ph3-ppii/cs-ha-medi.dat", dbpath); break;
      case NH: sprintf(str, "%s/ph3-ppii/cs-n-medi.dat", dbpath);  break;
    }
    m = fopen(str, "r");
    if (m==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}
  }

  sprintf(str, "%s/%s", dbpath, n5);
  par = fopen(str,"r");
  if (par==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  err = fscanf(uno, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  err = fscanf(due, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  err = fscanf(tre, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  err = fscanf(qua, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
 
  format = "%*s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf";
  for(i=0;i<21;i++)
  {
    err = fscanf(m, "%*s %lf %*f", &mean[i]);

    err = fscanf(uno, format, &s1[i][0], &s1[i][1], &s1[i][2], &s1[i][3], &s1[i][4], &s1[i][5], &s1[i][6], &s1[i][7], &s1[i][8], &s1[i][9], &s1[i][10], &s1[i][11], &s1[i][12], &s1[i][13], &s1[i][14], &s1[i][15], &s1[i][16], &s1[i][17], &s1[i][18], &s1[i][19], &s1[i][20]);

    err = fscanf(due, format, &s2[i][0], &s2[i][1], &s2[i][2], &s2[i][3], &s2[i][4], &s2[i][5], &s2[i][6], &s2[i][7], &s2[i][8], &s2[i][9], &s2[i][10], &s2[i][11], &s2[i][12], &s2[i][13], &s2[i][14], &s2[i][15], &s2[i][16], &s2[i][17], &s2[i][18], &s2[i][19], &s2[i][20]);

    err = fscanf(tre, format, &d1[i][0], &d1[i][1], &d1[i][2], &d1[i][3], &d1[i][4], &d1[i][5], &d1[i][6], &d1[i][7], &d1[i][8], &d1[i][9], &d1[i][10], &d1[i][11], &d1[i][12], &d1[i][13], &d1[i][14], &d1[i][15], &d1[i][16], &d1[i][17], &d1[i][18], &d1[i][19], &d1[i][20]);

    err = fscanf(qua, format, &d2[i][0], &d2[i][1], &d2[i][2], &d2[i][3], &d2[i][4], &d2[i][5], &d2[i][6], &d2[i][7], &d2[i][8], &d2[i][9], &d2[i][10], &d2[i][11], &d2[i][12], &d2[i][13], &d2[i][14], &d2[i][15], &d2[i][16], &d2[i][17], &d2[i][18], &d2[i][19], &d2[i][20]);
  }

  fclose(m);
  fclose(uno);
  fclose(due);
  fclose(tre);
  fclose(qua);

  err = fscanf(par, "%lf", &as1);
  err = fscanf(par, "%lf", &as2);
  err = fscanf(par, "%lf", &ad1);
  err = fscanf(par, "%lf", &ad2);
  
  fclose(par);

  for(i=0;i<sql;i++)
  {
    if(nseq[i]==G&&(ashift==HA||ashift==CB)) continue;
    if(nseq[i]==P&&(ashift==HN||ashift==NH)) continue;
    if(nseq[i]==Z) continue;

    tcs[i] = mean[nseq[i]]; 

    if(i>1&&nseq[i-2]!=Z&&nseq[i-1]!=Z)       tcs[i] += as2*s2[nseq[i-2]][nseq[i]];
    if(i>0&&nseq[i-1]!=Z)                     tcs[i] += as1*s1[nseq[i-1]][nseq[i]];
    if(i<(sql-1)&&nseq[i+1]!=Z)               tcs[i] += ad1*d1[nseq[i+1]][nseq[i]];
    if(i<(sql-2)&&nseq[i+1]!=Z&&nseq[i+2]!=Z) tcs[i] += ad2*d2[nseq[i+2]][nseq[i]];

    if(debug) {fprintf(stdout, "P %i_CS: %i %f\n", ashift, i, tcs[i]); fflush(stdout);}
  }

  return 0;
}

double do_predict(int ashift, char *bmrb, int ph, char *dbpath, char *out, int firstres, char *argv[], int argc)
{
  int    i, j, result, err, ti;
  int    start, end, by; 
  int    coeff;
  static int    sql=0;
  double sigma[4];
#ifndef LIB
  int    mtot;
  int    *star, *nseq;
  double wc[6][6];
  double **ecs, **intp;
  double *alphacs, *betacs, *coilcs, *ppiics;
  double **Aprob, **Bprob, **Cprob, **Pprob;
  double **pp, **spp;
  double mhelix, mbeta, mcoil, mppii;
  char   n7[256], backup[256], namecstable[256];
  FILE   *pred, *cstable=NULL;
#else
  static int    *star, *nseq;
  static double *alphacs, *betacs, *coilcs, *ppiics;
  static double **Aprob, **Bprob, **Cprob, **Pprob;
  static double **pp, **spp;
#endif
  double norm, prec, rerr, tmp;
  FILE   *ref, *tab=NULL;
  static char   seq[MAXLENGTH];
  char   *format=NULL, ss[MAXLENGTH], css[MAXLENGTH], str[MAXLENGTH];
  char   *line=NULL;
  size_t linecap = 0;
  ssize_t linelen;

  /* allocate everything */
#ifdef LIB
  if(!count) {
#endif
  nseq =    (int *) calloc (MAXLENGTH,sizeof(int));
  star =    (int *) calloc (MAXLENGTH,sizeof(int));
  alphacs = (double *) calloc(MAXLENGTH,sizeof(double));
  betacs =  (double *) calloc(MAXLENGTH,sizeof(double));
  coilcs =  (double *) calloc(MAXLENGTH,sizeof(double));
  ppiics =  (double *) calloc(MAXLENGTH,sizeof(double));
  ecs =  (double **) calloc(NSHIFT,sizeof(double*));
  intp = (double **) calloc(MAXLENGTH,sizeof(double*));
  for(i=0;i<NSHIFT;i++)    ecs[i] =  (double *) calloc(MAXLENGTH,sizeof(double));
  for(i=0;i<MAXLENGTH;i++) intp[i] = (double *) calloc(4,sizeof(double*));
  Aprob =   (double **) calloc(MAXLENGTH,sizeof(double*));
  Bprob =   (double **) calloc(MAXLENGTH,sizeof(double*));
  Cprob =   (double **) calloc(MAXLENGTH,sizeof(double*));
  Pprob =   (double **) calloc(MAXLENGTH,sizeof(double*));
  pp    =   (double **) calloc(MAXLENGTH,sizeof(double*));
  spp   =   (double **) calloc(MAXLENGTH,sizeof(double*));
  for(i=0;i<MAXLENGTH;i++) {
    Aprob[i] =   (double *) calloc(NSHIFT,sizeof(double));
    Bprob[i] =   (double *) calloc(NSHIFT,sizeof(double));
    Cprob[i] =   (double *) calloc(NSHIFT,sizeof(double));
    Pprob[i] =   (double *) calloc(NSHIFT,sizeof(double));
    pp[i]    =   (double *) calloc(4,sizeof(double));
    spp[i]   =   (double *) calloc(4,sizeof(double));
  }
#ifdef LIB
  }
#endif

#ifndef LIB
  if(debug) {fprintf(stdout, "In do predict: ashift is %i bmrb is %s out is %s ph is %i\n", ashift, bmrb, out, ph); fflush(stdout);}

  for(j=0;j<MAXLENGTH;j++) {seq[j]=0; css[j]=0; ss[j]=0; }

  if(ashift==NSHIFT) {start=0; end=NH; by=1; strcpy(n7,"SS-results.dat");}
  else {
    start = ashift; 
    end = ashift;
    by = 1;
    switch(ashift)
    {
      case CA: strcpy(n7,"SS-ca.dat"); break;
      case CB: strcpy(n7,"SS-cb.dat"); break;
      case CO: strcpy(n7,"SS-co.dat"); break;
      case HN: strcpy(n7,"SS-hn.dat"); break;
      case HA: strcpy(n7,"SS-ha.dat"); break;
      case NH: strcpy(n7,"SS-nh.dat"); break;
    }
  }

  if(out!=NULL) strcpy(n7,out);
  if(debug) {fprintf(stdout, "In do predict: out is %s n7 is %s \n", out, n7); fflush(stdout);}
  sprintf(backup, "%s.old", n7);
  rename(n7, backup);
  sprintf(str, "%s/other-db/weights.tab", dbpath);

  if(debug) {fprintf(stdout, "In do predict: I'm going to open %s ", str); fflush(stdout);}
  tab = fopen(str,"r");
  if (tab==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}
  err = fscanf(tab, "%*s %*s %*s %*s %*s %*s %*s");
  for(i=0;i<6;i++)
    err = fscanf(tab, "%*s %lf %lf %lf %lf %lf %lf", &wc[i][0], &wc[i][1], &wc[i][2], &wc[i][3], &wc[i][4], &wc[i][5]);
  fclose(tab);
  if(debug) {fprintf(stdout, "done\n"); fflush(stdout);}
#else 
  start=0; end=NH; by=1;
#endif

#ifdef LIB
  if(!count) {
    format = "%i %*i %c %*c %c %lf %lf %lf %lf %lf %lf";
#else
  if(debug) {fprintf(stdout, "In do_predict I am going to open %s ...", bmrb); fflush(stdout);}
  if(shifty) format = "%i %c %lf %lf %lf %lf %lf %lf"; 
  else format = "%i %*i %c %*c %c %lf %lf %lf %lf %lf %lf";
#endif 
  ref = fopen(bmrb,"r");
  if (ref==NULL) {fprintf(stderr, "%s file not found!\n", bmrb); return 1;}

  sql=0;
  while ((linelen = getline(&line, &linecap, ref)) != EOF) 
  {
    if(line[0]=='#'||linelen<=1) continue;
    if(!shifty) { 
      err = sscanf(line, format, &ti, &seq[sql], &css[sql], &ecs[0][sql], &ecs[1][sql], &ecs[2][sql], &ecs[3][sql], &ecs[4][sql], &ecs[5][sql]);
      if(!sql&&firstres==-999) firstres=ti;
      if(err!=9) {fprintf(stderr, "ERROR: WRONG FILE FORMAT AROUND LINE %i\n", sql+1); exit(1); }
    } else {
      err = sscanf(line, format, &ti, &seq[sql], &ecs[4][sql], &ecs[0][sql], &ecs[1][sql], &ecs[2][sql], &ecs[5][sql], &ecs[3][sql]); 
      if(err!=8) {fprintf(stderr, "ERROR: WRONG FILE FORMAT AROUND LINE %i\n", sql+1); exit(1); }
      css[sql]='U'; 
      if(!sql&&firstres==-999) firstres=ti;
    }
    sql++;
  }
  fclose(ref);
#ifdef LIB
  }
#endif

#ifndef LIB
  fprintf(stderr,"Protein length is %i. Last residue is %c\n", sql, seq[sql-1]); fflush(stderr);
#endif

  for(i=0;i<sql;i++) {
    switch(seq[i]) 
    {
      case 'A': nseq[i] = A; break;
      case 'C': nseq[i] = C; break;
      case 'D': nseq[i] = D; break;
      case 'E': nseq[i] = E; break;
      case 'F': nseq[i] = F; break;
      case 'G': nseq[i] = G; break;
      case 'H': nseq[i] = H; break;
      case 'I': nseq[i] = I; break;
      case 'K': nseq[i] = K; break;
      case 'L': nseq[i] = L; break;
      case 'M': nseq[i] = M; break;
      case 'N': nseq[i] = N; break;
      case 'O': nseq[i] = P; break; /* camcoil */
      case 'P': nseq[i] = P; break;
      case 'Q': nseq[i] = Q; break;
      case 'R': nseq[i] = R; break;
      case 'S': nseq[i] = S; break;
      case 'T': nseq[i] = T; break;
      case 'V': nseq[i] = V; break;
      case 'W': nseq[i] = W; break;
      case 'Y': nseq[i] = Y; break;
      case 'X': nseq[i] = X; break;
      case 'Z': nseq[i] = Z; break;
      default: fprintf(stderr, "ERROR: AminoAcids %c not recognized! line %d\n", seq[i], i); return 1; break;
    }
  }

  for(ashift=start;ashift<(end+1);ashift=ashift+by)
  {
    sprintf(str, "%s/other-db/tab.err", dbpath);
    if(debug) {fprintf(stdout, "In do predict: I'm going to open %s ", str); fflush(stdout);}
    tab = fopen(str,"r");
    if (tab==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

    err = fscanf(tab, "%*s %*s %*s %*s %*s");
    for(i=0;i<NSHIFT;i++)
    {
      if(i==ashift) err = fscanf(tab, "%*s %lf %lf %lf %lf", &sigma[0], &sigma[1], &sigma[2], &sigma[3]);
      else err = fscanf(tab, "%*s %*f %*f %*f %*f");
    }
    fclose(tab);
    if(debug) {fprintf(stdout, "done\n"); fflush(stdout);}

    for(i=0;i<MAXLENGTH;i++) { alphacs[i]=0.; betacs[i]=0.; coilcs[i]=0.; ppiics[i]=0.; }

    if(debug) {fprintf(stdout, "In do predict before calling do_alpha\n"); fflush(stdout);}
    result = do_alpha(ashift, alphacs, nseq, sql, ph, dbpath);
    if(result) return result;

#ifdef TAMIOLA
    if(debug) {fprintf(stdout, "In do predict before calling do_tamiola\n"); fflush(stdout);}
    result = do_tamiola(ashift, coilcs, nseq, sql, ph, dbpath);
    if(result) return result;
#else
    if(debug) {fprintf(stdout, "In do predict before calling do_coil\n"); fflush(stdout);}
    result = do_coil(ashift, coilcs, nseq, sql, ph, dbpath);
    if(result) return result;
#endif

    if(debug) {fprintf(stdout, "In do predict before calling do_ppii\n"); fflush(stdout);}
    result = do_ppii(ashift, ppiics, nseq, sql, ph, dbpath);
    if(result) return result;

    if(debug) {fprintf(stdout, "In do predict before calling do_beta\n"); fflush(stdout);}
    result = do_beta(ashift, betacs, nseq, sql, ph, dbpath);
    if(result) return result;

#ifndef LIB
    if(debug) {
      switch(ashift)
      {
        case CA: strcpy(namecstable,"pred-ca.dat"); break;
        case CB: strcpy(namecstable,"pred-cb.dat"); break;
        case CO: strcpy(namecstable,"pred-co.dat"); break;
        case HN: strcpy(namecstable,"pred-hn.dat"); break;
        case HA: strcpy(namecstable,"pred-ha.dat"); break;
        case NH: strcpy(namecstable,"pred-nh.dat"); break;
      }
      cstable = fopen(namecstable, "w");
      for(i=0;i<sql;i++) 
         fprintf(cstable, "%i %c %6.2f %6.2f %6.2f %6.2f\n", i, seq[i], alphacs[i], coilcs[i], ppiics[i], betacs[i]);
      fclose(cstable);
    }
#endif

    for(i=0;i<sql;i++) 
    {

      Aprob[i][ashift] = 999.999;
      Bprob[i][ashift] = 999.999;
      Cprob[i][ashift] = 999.999;
      Pprob[i][ashift] = 999.999;

      if(nseq[i]==G&&(ashift==HA||ashift==CB)) continue;
      if(nseq[i]==P&&(ashift==HN||ashift==NH)) continue;
      if(nseq[i]==Z) continue;

      if(ecs[ashift][i]<900.&&ecs[ashift][i]>0.) {
        if(((ecs[ashift][i]>=alphacs[i]-wc[5][5]*sigma[0])||
            (ecs[ashift][i]>=betacs[i]-wc[5][5]*sigma[1])||
            (ecs[ashift][i]>=coilcs[i]-(wc[5][5]+1.0)*sigma[2])||
            (ecs[ashift][i]>=ppiics[i]-wc[5][5]*sigma[3]))&&
           ((ecs[ashift][i]<=alphacs[i]+wc[5][5]*sigma[0])||
            (ecs[ashift][i]<=betacs[i]+wc[5][5]*sigma[1])||
            (ecs[ashift][i]<=coilcs[i]+(wc[5][5]+1.0)*sigma[2])||
            (ecs[ashift][i]<=ppiics[i]+wc[5][5]*sigma[3]))) {

            if(ashift==CA||ashift==CO) {
              /* are alpha and beta well predicted? */
              /* the correct order is Helix > Coil > PPII > beta */
              if(alphacs[i] < coilcs[i]) Aprob[i][ashift] = (1./M_PI)/(((wc[4][3]*wc[4][3]+1.0)*sigma[0]*0.95));
              if( coilcs[i] < betacs[i]) Bprob[i][ashift] = (1./M_PI)/(((wc[4][4]*wc[4][4]+1.0)*sigma[1]*0.95));

              /* if coil-beta is smaller than something, then ppii can be slightly smaller than beta */
              if(ppiics[i] < betacs[i] && (coilcs[i]-betacs[i])>wc[5][3]*sigma[2])  
                                         Pprob[i][ashift] = (1./M_PI)/(((wc[4][5]*wc[4][5]+1.0)*sigma[3]*0.95));  
              /* if ppii is slightly greater than coil can still be accepted */
              if((ppiics[i]-wc[5][4]*sigma[3]) > coilcs[i])  
                                         Pprob[i][ashift] = (1./M_PI)/(((wc[4][5]*wc[4][5]+1.0)*sigma[3]*0.95));

            } else if(ashift==HA||ashift==CB||ashift==NH||ashift==HN) {
              /* the correct order is Beta > PPII > Coil > Helix */
              if(coilcs[i] < alphacs[i]) Aprob[i][ashift] = (1./M_PI)/(((wc[4][3]*wc[4][3]+1.0)*sigma[0]*0.95));
              if(betacs[i] <  coilcs[i]) Bprob[i][ashift] = (1./M_PI)/(((wc[4][4]*wc[4][4]+1.0)*sigma[1]*0.95));
              if(ppiics[i] > betacs[i] && (betacs[i]-coilcs[i])>wc[5][3]*sigma[2]) 
                                         Pprob[i][ashift] = (1./M_PI)/(((wc[4][5]*wc[4][5]+1.0)*sigma[3]*0.95));  
              if((ppiics[i]+wc[5][4]*sigma[3]) < coilcs[i])
                                         Pprob[i][ashift] = (1./M_PI)/(((wc[4][5]*wc[4][5]+1.0)*sigma[3]*0.95));
            }

#ifdef LORENTZ
         /* experimental data probabilities */
         if(Aprob[i][ashift]>900) Aprob[i][ashift] = (0.95*sigma[0]/M_PI)/((ecs[ashift][i]-alphacs[i])*(ecs[ashift][i]-alphacs[i])+(sigma[0]*sigma[0]*0.95*0.95));
         if(Bprob[i][ashift]>900) Bprob[i][ashift] = (0.95*sigma[1]/M_PI)/((ecs[ashift][i]-betacs[i] )*(ecs[ashift][i]-betacs[i] )+(sigma[1]*sigma[1]*0.95*0.95));
         if(Cprob[i][ashift]>900) Cprob[i][ashift] = (0.95*sigma[2]/M_PI)/((ecs[ashift][i]-coilcs[i] )*(ecs[ashift][i]-coilcs[i] )+(sigma[2]*sigma[2]*0.95*0.95));
         if(Pprob[i][ashift]>900) Pprob[i][ashift] = (0.95*sigma[3]/M_PI)/((ecs[ashift][i]-ppiics[i] )*(ecs[ashift][i]-ppiics[i] )+(sigma[3]*sigma[3]*0.95*0.95));
#else
         if(Aprob[i][ashift]>900) Aprob[i][ashift] = 1./(sigma[0]*sqrt(2.*M_PI))*exp(-(ecs[ashift][i]-alphacs[i])*(ecs[ashift][i]-alphacs[i])/(2.*sigma[0]*sigma[0]));
         if(Bprob[i][ashift]>900) Bprob[i][ashift] = 1./(sigma[1]*sqrt(2.*M_PI))*exp(-(ecs[ashift][i]-betacs[i] )*(ecs[ashift][i]-betacs[i] )/(2.*sigma[1]*sigma[1]));
         if(Cprob[i][ashift]>900) Cprob[i][ashift] = 1./(sigma[2]*sqrt(2.*M_PI))*exp(-(ecs[ashift][i]-coilcs[i] )*(ecs[ashift][i]-coilcs[i] )/(2.*sigma[2]*sigma[2]));
         if(Pprob[i][ashift]>900) Pprob[i][ashift] = 1./(sigma[3]*sqrt(2.*M_PI))*exp(-(ecs[ashift][i]-ppiics[i] )*(ecs[ashift][i]-ppiics[i] )/(2.*sigma[3]*sigma[3]));
#endif

#ifndef LIB
        } else {
          fprintf(stderr, "WARNING: Chemical Shift %6.2f %2s for Res %3i %c is unexpected! ", ecs[ashift][i], shtostring(ashift), i+firstres, seq[i]);
          if(alphacs[i]>betacs[i]) {
            fprintf(stderr, "(expected range is: %6.2f - %6.2f)\n", betacs[i]-wc[5][5]*sigma[1], alphacs[i]+wc[5][5]*sigma[0]);
          } else {
            fprintf(stderr, "(expected range is: %6.2f - %6.2f)\n", alphacs[i]-wc[5][5]*sigma[0], betacs[i]+wc[5][5]*sigma[1]);
          }
#endif
        }
      }

      if(debug) { fprintf(stderr, "Scratch PROBs A_%i CS_%i H_%6.3f E_%6.3f C_%6.3f P_%6.3f\n", 
                          i, ashift, Aprob[i][ashift], Bprob[i][ashift], Cprob[i][ashift], Pprob[i][ashift]); 
                  fflush(stderr); }

    }

  }

  for(i=0;i<sql;i++) 
  {
    coeff=0;
    for(ashift=start;ashift<(end+1);ashift=ashift+by)
    {
      if(Aprob[i][ashift]<900) coeff++;
    }

    if(!force) star[i]=coeff;
    else star[i]=5;
    if(i==0||i==sql-1) star[i]=0;

  }

#ifndef LIB
  pred = fopen(n7, "w");
 
  if(debug) { fprintf(stdout, "Summing probabilities...\n"); fflush(stdout);}
#endif

  for(i=0;i<sql;i++)
  {
    pp[i][0] = pp[i][1] = pp[i][2] = pp[i][3] = 0.;
    for(j=0;j<NSHIFT;j++)
    {
       if(Aprob[i][j]<900)
       {
         pp[i][0] += wc[0][j]*log(Aprob[i][j]);
         pp[i][1] += wc[1][j]*log(Bprob[i][j]);
         pp[i][2] += wc[2][j]*log(Cprob[i][j]);
         pp[i][3] += wc[3][j]*log(Pprob[i][j]);
         if(debug) {fprintf(stderr, "PP SHIFTs A_%i CS_%i H_%8.3f E_%8.3f C_%8.3f P_%8.3f\n", i, j, pp[i][0], pp[i][1], pp[i][2], pp[i][3]); fflush(stderr); }
       }
    }
    /* these are the intrinsic probability per aminoacids */
    if(pp[i][0]!=0.||pp[i][1]!=0.||pp[i][2]!=0.) {
      if(ppii) { 
        norm = exp(pp[i][0])+exp(pp[i][1])+exp(pp[i][2])+exp(pp[i][3]);
        pp[i][3] = exp(pp[i][3])/norm;
      } else {
        norm = exp(pp[i][0])+exp(pp[i][1])+exp(pp[i][2]);
      }
      pp[i][0] = exp(pp[i][0])/norm;
      pp[i][1] = exp(pp[i][1])/norm;
      pp[i][2] = exp(pp[i][2])/norm;

      if(debug) {
        fprintf(stderr, "PP TOTALs A_%i CS_%i H_%8.3f E_%8.3f C_%8.3f P_%8.3f\n", i, j, pp[i][0], pp[i][1], pp[i][2], pp[i][3]); 
        fflush(stderr);
      }
    }
  }

  do_averaging(pp, spp, nseq, sql, wc, star, ss);

  /* helix first residue check */
  for(i=2;i<sql-1;i++)
  {
    /*
    EHH -> if EEHH we don't do anything
        -> else if CEHH then it is easily the result of the downshift 
           of the first residue of an helix so we can invert E with H
    CHH -> and if P(E) > P(H), again we can invert P(E) with P(H) in pp.
    */
    if(ss[i]=='H'&&ss[i+1]=='H'&&ss[i-1]=='E'&&ss[i-2]=='E') continue;
    else if(ss[i]=='H'&&ss[i+1]=='H'&&ss[i-1]=='E'&&ss[i-2]=='C'&&spp[i][0]>0.6) {
      tmp = pp[i-1][0];
      pp[i-1][0] = pp[i-1][1];
      pp[i-1][1] = tmp;
    } 
    else if(ss[i]=='H'&&ss[i+1]=='H'&&ss[i-1]=='C'&&pp[i-1][1]>pp[i-1][0]&&spp[i][0]>0.6) {
      tmp = pp[i-1][0];
      pp[i-1][0] = pp[i-1][1];
      pp[i-1][1] = tmp;
    } 
  }
  /* do_averaging(pp, spp, nseq, sql, wc, star, ss); */
  /* done */

  j=0; prec=0;
  for(i=0;i<sql;i++)
  {
    if(pp[i][0]!=0||pp[i][1]!=0||pp[i][2]!=0) {
      if(star[i]>2){
        if(css[i]!='T'&&css[i]!='C'&&
           css[i]!='P'&&css[i]!='H'&&
           css[i]!='I'&&css[i]!='G'&&
           css[i]!='B'&&css[i]!='E') continue;

        if(css[i]=='T') css[i]='C';
        if(css[i]!=ss[i]) prec++;

        if(ss[i]=='H'&&css[i]=='I') prec--;
        if(ss[i]=='H'&&css[i]=='G') prec--;

        if(ss[i]=='C'&&css[i]=='G') prec-=0.50;

        if(ss[i]=='E'&&css[i]=='B') prec-=0.50;
        if(ss[i]=='C'&&css[i]=='B') prec--;

        j++;
      }
    }
  }

#ifndef LIB
  fprintf(stderr,"SQ:%s\n", seq);
  fprintf(stderr,"SS:%s\n", ss); 
  fprintf(stderr,"SS:%s\n", css);
  rerr = prec/j*100.; 
  fprintf(stdout,"Err %7.3f\n", rerr); fflush(stdout);

/* Calculating total populations */
  mhelix = mbeta = mcoil = mppii = 0.;
  mtot = 0;
  for(i=0;i<sql;i++)
  {
    if(spp[i][0]!=0||spp[i][1]!=0||spp[i][2]!=0) {
      mhelix += spp[i][0];
      mbeta += spp[i][1];
      mcoil += spp[i][2];
      mtot++;
      if(ppii) mppii += spp[i][3];
    }
  }
  mhelix = mhelix/mtot*100.;
  mbeta = mbeta/mtot*100.;
  mcoil = mcoil/mtot*100.;
  mppii = mppii/mtot*100.;

/* Writing the output file */
  fprintf(pred,"#d2D - v. 1.4.0\n#(c) Carlo Camilloni\n#PLEASE CITE:\n");
  fprintf(pred,"#Camilloni C., De Simone A., Vranken W., and Vendruscolo M.\n");
  fprintf(pred,"#Determination of Secondary Structure Populations in Disordered States of Proteins using NMR Chemical Shifts\n");
  fprintf(pred,"#Biochemistry 2012, 51: 2224-2231\n\n");
  fprintf(pred,"#Executed as: ");
  for (i=0;i<argc;i++) fprintf(pred, "%s ", argv[i]);
  fprintf(pred, "\n\n");
  fprintf(pred,"#SQ:%s\n",seq);
  fprintf(pred,"#SS:%s\n\n",ss);
  fprintf(pred,"#Total Populations:\n");
  fprintf(pred,"#Helix(H): %4.1f%%\n", mhelix);
  fprintf(pred,"#Extended-Beta(E): %4.1f%%\n", mbeta);
  if(ppii) fprintf(pred,"#Polyproline II (PPII)(P): %4.1f%%\n", mppii);
  fprintf(pred,"#Coil(C): %4.1f%%\n\n", mcoil);
  fprintf(pred,"#Populations per residue (residues marked with a * are less reliable):\n");
  if(ppii) fprintf(pred,"#num \t res \t      Helix      Beta       Coil       PPII  SS \n");
  else     fprintf(pred,"#num \t res \t      Helix      Beta       Coil  SS\n");

  for(i=0;i<sql;i++)
  {
    if((pp[i][0]!=0||pp[i][1]!=0||pp[i][2]!=0)&&star[i]>2) {
      if(!ppii) 
      fprintf(pred, "%i \t %c \t %10.3f %10.3f %10.3f %c\n", i+firstres, seq[i], spp[i][0], spp[i][1], spp[i][2], ss[i]);
      else
      fprintf(pred, "%i \t %c \t %10.3f %10.3f %10.3f %10.3f %c\n", i+firstres, seq[i], spp[i][0], spp[i][1], spp[i][2], spp[i][3], ss[i]);
    } else if(star[i]<3&&ss[i]!=' ') {
      if(!ppii) 
      fprintf(pred, "%i \t %c \t %10.3f %10.3f %10.3f %c*\n", i+firstres, seq[i], spp[i][0], spp[i][1], spp[i][2], ss[i]);
      else
      fprintf(pred, "%i \t %c \t %10.3f %10.3f %10.3f %10.3f %c*\n", i+firstres, seq[i], spp[i][0], spp[i][1], spp[i][2], spp[i][3], ss[i]);
    } else {
      fprintf(pred, "#%i \t %c \t \n", i+firstres, seq[i]);
    }
  }
  fprintf(pred,"#DONE!\n");
  fprintf(pred,"#d2D - v. 1.4.0\n#(c) Carlo Camilloni\n#PLEASE CITE:\n");
  fprintf(pred,"#Camilloni C., De Simone A., Vranken W., and Vendruscolo M.\n");
  fprintf(pred,"#Determination of Secondary Structure Populations in Disordered States of Proteins using NMR Chemical Shifts\n");
  fprintf(pred,"#Biochemistry 2012, 51: 2224-2231\n");
  fclose(pred);
#endif

#ifdef LIB
  count++;
  if(!count) {
#endif
  free(nseq);
  free(star);
  for(i=0;i<NSHIFT;i++) free(ecs[i]);
  for(i=0;i<MAXLENGTH;i++) free(intp[i]);
  free(ecs);
  free(intp);
  free(alphacs);
  free(betacs);
  free(coilcs);
  free(ppiics);
  for(i=0;i<MAXLENGTH;i++) {
    free(Aprob[i]);
    free(Bprob[i]);
    free(Cprob[i]);
    free(Pprob[i]);
    free(pp[i]);
    free(spp[i]);
  }
  free(Aprob);
  free(Bprob);
  free(Cprob);
  free(Pprob);
  free(pp);
  free(spp); 
#ifdef LIB
  }
#endif

  rerr = (double)prec/j*100.;
  return rerr;
}

#ifndef LIB
void help()
{
  fprintf(stderr,"\nd2D - v. 1.4.0\n(c) Carlo Camilloni\n\nPLEASE CITE:\n");
  fprintf(stderr,"Camilloni C., De Simone A., Vranken W., and Vendruscolo M.\n");
  fprintf(stderr,"Determination of Secondary Structure Populations in Disordered States of Proteins\n");
  fprintf(stderr,"using NMR Chemical Shifts\n");
  fprintf(stderr,"Biochemistry 2012, 51: 2224-2231\n\n");
  fprintf(stderr,"\t-pH        \t([neutral], acid)\n");
  fprintf(stderr,"\t-file      \t(input file name)\n");
  fprintf(stderr,"\t-out       \t(output file name)\n");
  fprintf(stderr,"\t-fres      \t(rescale residue numbers from fres)\n");
  fprintf(stderr,"\t-noppii    \t(turn off the PolyProline II predictor)\n");
  fprintf(stderr,"\t-shifty    \t(input file in the shifty format (website format))\n");
  fprintf(stderr,"\t-force_few_res \t(accept less than 3 chemical shifts per residue as input, DANGEROUS!)\n");
  fprintf(stderr,"\t-debug     \t(more verbose)\n");
  fprintf(stderr,"\t-help      \t(here we are!)\n\n");
  fprintf(stderr,"EXAMPLE: (shifty format)\n");
  fprintf(stderr,"#NUM	AA	HA	CA	CB	CO	N	HN   \n");
  fprintf(stderr,"1	P	4.32	62.69	32.89	0	0	0    \n");
  fprintf(stderr,"2	N	5.03	52.39	38.99	174.40	0	0    \n");
  fprintf(stderr,"3	F	4.32	60.09	31.99	176.00	121.78	9.82 \n");
  fprintf(stderr,"4	S	4.30	60.39	64.19	174.30	112.38	8.49 \n");
  fprintf(stderr,"5	G	3.96	44.59	0	170.40	110.48	9.22 \n");
  fprintf(stderr,"6	N	5.52	52.69	39.89	174.00	118.48	8.15 \n");
  fprintf(stderr,"7	E	5.14	56.39	31.89	175.00	123.18	9.35 \n");
  fprintf(stderr,"8	K	5.28	53.79	35.79	174.90	123.08	10.17\n");
  fprintf(stderr,"9	I	3.92	55.19	33.29	175.80	125.48	9.05 \n\n");
  fprintf(stderr,"EXAMPLE: (standard format)\n");
  fprintf(stderr,"1 1	P X C	62.69	32.89	999.99	999.9 4.32 999.99\n");
  fprintf(stderr,"2 2	N X C	52.39	38.99	174.40	999.9 5.03 999.99\n");
  fprintf(stderr,"3 3	F X C	60.09	31.99	176.00	9.82  4.32 121.78\n");
  fprintf(stderr,"4 4	S X E	60.39	64.19	174.30	8.49  4.30 112.38\n");
  fprintf(stderr,"5 5	G X E	44.59	999.9	170.40	9.22  3.96 110.48\n");
  fprintf(stderr,"6 6	N X E	52.69	39.89	174.00	8.15  5.52 118.48\n");
  fprintf(stderr,"7 7	E X C	56.39	31.89	175.00	9.35  5.14 123.18\n");
  fprintf(stderr,"8 8	K X C	53.79	35.79	174.90	10.17 5.28 123.08\n");
  fprintf(stderr,"9 9	I X C	55.19	33.29	175.80	9.05  3.92 125.48\n\n");
}

int main(int argc, char *argv[])
{
  int ashift=NSHIFT;
  int status=0;
  int ph=0;
  int firstres=-999;
  int i, j=-1, o=-1;
  double dst=0.;
  char *dbpath;

  dbpath = getenv ("CAMDBV1");
  if (dbpath==NULL) {fprintf (stderr, "The CAMDBV1 is not set!\n"); return 1;}
  if(argc<3) {help(); return 0;}

  for(i=0;i<argc;i++)
  {
    if(!strcmp(argv[i],"-pH")) {
      if(!strcmp(argv[i+1],"acid")) ph = 1;
      else if(!strcmp(argv[i+1],"neutral")) ph = 0;
    }
    if(!strcmp(argv[i],"-file")) {
      if(argv[i+1]!=NULL) j=i+1;
    }
    if(!strcmp(argv[i],"-out")) {
      if(argv[i+1]!=NULL) o=i+1;
    }
    if(!strcmp(argv[i],"-debug")) debug=1;
    if(!strcmp(argv[i],"-noppii")) ppii=0;
    if(!strcmp(argv[i],"-fres")) firstres = atoi(argv[i+1]);
    if(!strcmp(argv[i],"-shifty")) shifty=1;
    if(!strcmp(argv[i],"-force_few_res")) force=1;
    if(!strcmp(argv[i],"-help")) {help(); return 0;}
  }

  if(debug) {fprintf(stdout, "In main: CAMDBV1 is %s\n", dbpath); fflush(stdout);}
  if(j==-1) {help(); return 0;}

    
  if(argv[j]!=NULL&&o!=-1) dst = do_predict(ashift, argv[j], ph, dbpath, argv[o], firstres, argv, argc);
  else if(argv[j]!=NULL&&o==-1) dst = do_predict(ashift, argv[j], ph, dbpath, NULL, firstres, argv, argc);
  else {fprintf(stderr, "Wrong file name: %s, expected one\n", argv[j]); return 1;}

  if(dst>=0.0) status = 0;
  return status;
}
#endif  
