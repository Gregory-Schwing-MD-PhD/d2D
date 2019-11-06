/* d2D - v. 2.1.0
 * (c) Carlo Camilloni
 * Camilloni C., De Simone A., Vranken W., and Vendruscolo M.
 * Determination of Secondary Structure Populations in Disordered States of Proteins using NMR Chemical Shifts
 * Biochemistry 2012, 51: 2224-2231 
 */

/* 
 * compile with
 * gcc d2D-predictor.c -lm -O3 -Wall -Wextra -ansi -pedantic -o d2D.x -I.  
 */
 
#include "d2D.h"

void do_averaging(double *pp[], double *spp[], int *nseq, int sql, double wc[5][6], int *star, char ss[MAXLENGTH])
{
  int    i, j, MAX;
  double dcount, maxi, norm;

  for(i=0;i<sql;i++)
  {
    if(pp[i][0]!=0.||pp[i][1]!=0.||pp[i][2]!=0.) 
    {
      if(star[i]<4) {
        if(((i>0)&&(i<sql-1))&&star[i-1]>2&&star[i+1]>3&&
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

        if(i>2) if(nseq[i-1]!=Z&&nseq[i-2]!=Z&&(pp[i-2][0]!=0.||pp[i-2][1]!=0.||pp[i-2][2]!=0.)&&star[i-2]>3) {
          spp[i][0] += wc[4][0]*pp[i-2][0];
          spp[i][1] += wc[4][0]*pp[i-2][1];
          spp[i][2] += wc[4][0]*pp[i-2][2];
          spp[i][3] += wc[4][0]*pp[i-2][3];
          dcount+=wc[4][0];
        }

        if(i>1) if(nseq[i-1]!=Z&&(pp[i-1][0]!=0.||pp[i-1][1]!=0.||pp[i-1][2]!=0.)&&star[i-1]>3) {
          spp[i][0] += wc[4][1]*pp[i-1][0];
          spp[i][1] += wc[4][1]*pp[i-1][1];
          spp[i][2] += wc[4][1]*pp[i-1][2];
          spp[i][3] += wc[4][1]*pp[i-1][3];
          dcount+=wc[4][1];
        }

        if(i<sql-2) if(nseq[i+1]!=Z&&(pp[i+1][0]!=0.||pp[i+1][1]!=0.||pp[i+1][2]!=0.)&&star[i+1]>3) {
          spp[i][0] += wc[4][1]*pp[i+1][0];
          spp[i][1] += wc[4][1]*pp[i+1][1];
          spp[i][2] += wc[4][1]*pp[i+1][2];
          spp[i][3] += wc[4][1]*pp[i+1][3];
          dcount+=wc[4][1];
        }

        if(i<sql-3) if(nseq[i+1]!=Z&&nseq[i+2]!=Z&&(pp[i+2][0]!=0.||pp[i+2][1]!=0.||pp[i+2][2]!=0.)&&star[i+2]>3) {
          spp[i][0] += wc[4][0]*pp[i+2][0];
          spp[i][1] += wc[4][0]*pp[i+2][1];
          spp[i][2] += wc[4][0]*pp[i+2][2];
          spp[i][3] += wc[4][0]*pp[i+2][3];
          dcount+=wc[4][0];
        }
      }
      spp[i][0] /= dcount;
      spp[i][1] /= dcount;
      spp[i][2] /= dcount;
      spp[i][3] /= dcount;

      norm = spp[i][0] + spp[i][1] + spp[i][2] + spp[i][3];
      spp[i][0] /= norm;
      spp[i][1] /= norm;
      spp[i][2] /= norm;
      spp[i][3] /= norm;

      maxi = spp[i][0]; MAX=0;
      for(j=1;j<4;j++) if( maxi < spp[i][j] ) {maxi = spp[i][j]; MAX=j;}

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
  char   n1[25], n2[25], n3[25], n4[25], n5[25], *format, str[MAXLENGTH];
  double mean[21];
  double s1[21][21], s2[21][21], d1[21][21], d2[21][21];
  double as1, as2, ad1, ad2;
  int    i;
  FILE   *m=NULL, *uno, *due, *tre, *qua, *par;

  switch(ashift)
  {
    case CA: strcpy(n1,"/helix-db/CA-uno.dat"); 
             strcpy(n2,"/helix-db/CA-due.dat"); 
             strcpy(n3,"/helix-db/CA-tre.dat"); 
             strcpy(n4,"/helix-db/CA-quattro.dat"); 
             strcpy(n5,"/helix-db/CA-parm.dat"); 
             break;

    case CB: strcpy(n1,"/helix-db/CB-uno.dat"); 
             strcpy(n2,"/helix-db/CB-due.dat"); 
             strcpy(n3,"/helix-db/CB-tre.dat"); 
             strcpy(n4,"/helix-db/CB-quattro.dat"); 
             strcpy(n5,"/helix-db/CB-parm.dat"); 
             break;

    case CO: strcpy(n1,"/helix-db/CO-uno.dat"); 
             strcpy(n2,"/helix-db/CO-due.dat"); 
             strcpy(n3,"/helix-db/CO-tre.dat"); 
             strcpy(n4,"/helix-db/CO-quattro.dat"); 
             strcpy(n5,"/helix-db/CO-parm.dat");  
             break;

    case HN: strcpy(n1,"/helix-db/HN-uno.dat"); 
             strcpy(n2,"/helix-db/HN-due.dat"); 
             strcpy(n3,"/helix-db/HN-tre.dat"); 
             strcpy(n4,"/helix-db/HN-quattro.dat"); 
             strcpy(n5,"/helix-db/HN-parm.dat"); 
             break;

    case HA: strcpy(n1,"/helix-db/HA-uno.dat"); 
             strcpy(n2,"/helix-db/HA-due.dat"); 
             strcpy(n3,"/helix-db/HA-tre.dat"); 
             strcpy(n4,"/helix-db/HA-quattro.dat"); 
             strcpy(n5,"/helix-db/HA-parm.dat"); 
             break;

    case NH: strcpy(n1,"/helix-db/NH-uno.dat"); 
             strcpy(n2,"/helix-db/NH-due.dat"); 
             strcpy(n3,"/helix-db/NH-tre.dat"); 
             strcpy(n4,"/helix-db/NH-quattro.dat"); 
             strcpy(n5,"/helix-db/NH-parm.dat"); 
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

  sprintf(str, "%s%s", dbpath, n5);
  par = fopen(str,"r");
  if (par==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}

  fscanf(uno, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  fscanf(due, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  fscanf(tre, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  fscanf(qua, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
 
  format = "%*s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf";
  for(i=0;i<21;i++)
  {
    fscanf(m, "%*s %lf %*f", &mean[i]);

    fscanf(uno, format, &s1[i][0],  &s1[i][1],  &s1[i][2],  &s1[i][3],  &s1[i][4],  &s1[i][5], 
                        &s1[i][6],  &s1[i][7],  &s1[i][8],  &s1[i][9],  &s1[i][10], &s1[i][11], 
                        &s1[i][12], &s1[i][13], &s1[i][14], &s1[i][15], &s1[i][16], &s1[i][17], 
                        &s1[i][18], &s1[i][19], &s1[i][20]);

    fscanf(due, format, &s2[i][0],  &s2[i][1],  &s2[i][2],  &s2[i][3],  &s2[i][4],  &s2[i][5], 
                        &s2[i][6],  &s2[i][7],  &s2[i][8],  &s2[i][9],  &s2[i][10], &s2[i][11], 
                        &s2[i][12], &s2[i][13], &s2[i][14], &s2[i][15], &s2[i][16], &s2[i][17], 
                        &s2[i][18], &s2[i][19], &s2[i][20]);

    fscanf(tre, format, &d1[i][0],  &d1[i][1],  &d1[i][2],  &d1[i][3],  &d1[i][4],  &d1[i][5], 
                        &d1[i][6],  &d1[i][7],  &d1[i][8],  &d1[i][9],  &d1[i][10], &d1[i][11], 
                        &d1[i][12], &d1[i][13], &d1[i][14], &d1[i][15], &d1[i][16], &d1[i][17], 
                        &d1[i][18], &d1[i][19], &d1[i][20]);

    fscanf(qua, format, &d2[i][0],  &d2[i][1],  &d2[i][2],  &d2[i][3],  &d2[i][4],  &d2[i][5], 
                        &d2[i][6],  &d2[i][7],  &d2[i][8],  &d2[i][9],  &d2[i][10], &d2[i][11], 
                        &d2[i][12], &d2[i][13], &d2[i][14], &d2[i][15], &d2[i][16], &d2[i][17], 
                        &d2[i][18], &d2[i][19], &d2[i][20]);

  }

  fclose(m);
  fclose(uno);
  fclose(due);
  fclose(tre);
  fclose(qua);

  fscanf(par, "%lf", &as1);
  fscanf(par, "%lf", &as2);
  fscanf(par, "%lf", &ad1);
  fscanf(par, "%lf", &ad2);
  
  fclose(par);

  for(i=0;i<sql;i++)
  {
    if(nseq[i]==G&&(ashift==HA||ashift==CB)) continue;
    if(nseq[i]==P&&(ashift==HN||ashift==NH)) continue;
    if(nseq[i]==Z) continue;

    tcs[i] = mean[nseq[i]];

    if(i>0&&nseq[i-1]!=Z) tcs[i] += as1*s1[nseq[i-1]][nseq[i]];/* i-1 */
    if(i<(sql-1)&&nseq[i+1]!=Z) tcs[i] += ad1*d1[nseq[i+1]][nseq[i]]; /* i+1 */
    if(i<(sql-3)&&nseq[i+3]!=Z&&nseq[i+2]!=Z&&nseq[i+1]!=Z) tcs[i] += as2*s2[nseq[i+3]][nseq[i]]; /* i+3 */
    if(i<(sql-4)&&nseq[i+4]!=Z&&nseq[i+3]!=Z&&nseq[i+2]!=Z&&nseq[i+1]!=Z) tcs[i] += ad2*d2[nseq[i+4]][nseq[i]]; /* i+4 */
  }

  return 0;
}

int do_beta(int ashift, double *tcs, int *nseq, int sql, int ph, char *dbpath)
{
  char   n1[25], n2[25], n3[25], n4[25], n5[25], *format, str[MAXLENGTH];
  double mean[21];
  double s1[21][21], s2[21][21], d1[21][21], d2[21][21];
  double as1, as2, ad1, ad2;
  int    i;
  FILE   *m=NULL, *uno, *due, *tre, *qua, *par;

  switch(ashift)
  {
    case CA: 
      strcpy(n1,"beta-db/CA-uno.dat"); 
      strcpy(n2,"beta-db/CA-due.dat"); 
      strcpy(n3,"beta-db/CA-tre.dat"); 
      strcpy(n4,"beta-db/CA-quattro.dat"); 
      strcpy(n5,"beta-db/CA-parm.dat"); 
      break;

    case CB: 
      strcpy(n1,"beta-db/CB-uno.dat"); 
      strcpy(n2,"beta-db/CB-due.dat"); 
      strcpy(n3,"beta-db/CB-tre.dat"); 
      strcpy(n4,"beta-db/CB-quattro.dat"); 
      strcpy(n5,"beta-db/CB-parm.dat"); 
      break;

    case CO: 
      strcpy(n1,"beta-db/CO-uno.dat"); 
      strcpy(n2,"beta-db/CO-due.dat"); 
      strcpy(n3,"beta-db/CO-tre.dat"); 
      strcpy(n4,"beta-db/CO-quattro.dat"); 
      strcpy(n5,"beta-db/CO-parm.dat"); 
      break;

    case HN: 
      strcpy(n1,"beta-db/HN-uno.dat"); 
      strcpy(n2,"beta-db/HN-due.dat"); 
      strcpy(n3,"beta-db/HN-tre.dat"); 
      strcpy(n4,"beta-db/HN-quattro.dat"); 
      strcpy(n5,"beta-db/HN-parm.dat"); 
      break;

    case HA: 
      strcpy(n1,"beta-db/HA-uno.dat"); 
      strcpy(n2,"beta-db/HA-due.dat"); 
      strcpy(n3,"beta-db/HA-tre.dat"); 
      strcpy(n4,"beta-db/HA-quattro.dat"); 
      strcpy(n5,"beta-db/HA-parm.dat"); 
      break;

    case NH: 
      strcpy(n1,"beta-db/NH-uno.dat"); 
      strcpy(n2,"beta-db/NH-due.dat"); 
      strcpy(n3,"beta-db/NH-tre.dat"); 
      strcpy(n4,"beta-db/NH-quattro.dat"); 
      strcpy(n5,"beta-db/NH-parm.dat"); 
      break;
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

  fscanf(uno, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  fscanf(due, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  fscanf(tre, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  fscanf(qua, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
 
  format = "%*s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf";
  for(i=0;i<21;i++)
  {
    fscanf(m, "%*s %lf %*f", &mean[i]);

    fscanf(uno, format, &s1[i][0],  &s1[i][1],  &s1[i][2],  &s1[i][3],  &s1[i][4],  &s1[i][5],  &s1[i][6],  &s1[i][7],  &s1[i][8],  &s1[i][9], 
                        &s1[i][10], &s1[i][11], &s1[i][12], &s1[i][13], &s1[i][14], &s1[i][15], &s1[i][16], &s1[i][17], &s1[i][18], &s1[i][19], 
                        &s1[i][20]);

    fscanf(due, format, &s2[i][0],  &s2[i][1],  &s2[i][2],  &s2[i][3],  &s2[i][4],  &s2[i][5],  &s2[i][6],  &s2[i][7],  &s2[i][8],  &s2[i][9], 
                        &s2[i][10], &s2[i][11], &s2[i][12], &s2[i][13], &s2[i][14], &s2[i][15], &s2[i][16], &s2[i][17], &s2[i][18], &s2[i][19], 
                        &s2[i][20]);

    fscanf(tre, format, &d1[i][0],  &d1[i][1],  &d1[i][2],  &d1[i][3],  &d1[i][4],  &d1[i][5],  &d1[i][6],  &d1[i][7],  &d1[i][8],  &d1[i][9], 
                        &d1[i][10], &d1[i][11], &d1[i][12], &d1[i][13], &d1[i][14], &d1[i][15], &d1[i][16], &d1[i][17], &d1[i][18], &d1[i][19], 
                        &d1[i][20]);

    fscanf(qua, format, &d2[i][0],  &d2[i][1],  &d2[i][2],  &d2[i][3],  &d2[i][4],  &d2[i][5],  &d2[i][6],  &d2[i][7],  &d2[i][8],  &d2[i][9], 
                        &d2[i][10], &d2[i][11], &d2[i][12], &d2[i][13], &d2[i][14], &d2[i][15], &d2[i][16], &d2[i][17], &d2[i][18], &d2[i][19], 
                        &d2[i][20]);
  }

  fclose(m);
  fclose(uno);
  fclose(due);
  fclose(tre);
  fclose(qua);

  fscanf(par, "%lf", &as1);
  fscanf(par, "%lf", &as2);
  fscanf(par, "%lf", &ad1);
  fscanf(par, "%lf", &ad2);
  
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
  int    i;
  FILE   *m=NULL, *uno, *due, *tre, *qua, *par;

  switch(ashift)
  {
    case CA: 
      strcpy(n1,"coil-db/CA-uno.dat"); 
      strcpy(n2,"coil-db/CA-due.dat"); 
      strcpy(n3,"coil-db/CA-tre.dat"); 
      strcpy(n4,"coil-db/CA-quattro.dat"); 
      strcpy(n5,"coil-db/CA-parm.dat"); 
      break;

    case CB: 
      strcpy(n1,"coil-db/CB-uno.dat"); 
      strcpy(n2,"coil-db/CB-due.dat"); 
      strcpy(n3,"coil-db/CB-tre.dat"); 
      strcpy(n4,"coil-db/CB-quattro.dat"); 
      strcpy(n5,"coil-db/CB-parm.dat"); 
      break;

    case CO: 
      strcpy(n1,"coil-db/CO-uno.dat"); 
      strcpy(n2,"coil-db/CO-due.dat"); 
      strcpy(n3,"coil-db/CO-tre.dat"); 
      strcpy(n4,"coil-db/CO-quattro.dat"); 
      strcpy(n5,"coil-db/CO-parm.dat"); 
      break;

    case HN: 
      strcpy(n1,"coil-db/HN-uno.dat"); 
      strcpy(n2,"coil-db/HN-due.dat"); 
      strcpy(n3,"coil-db/HN-tre.dat"); 
      strcpy(n4,"coil-db/HN-quattro.dat"); 
      strcpy(n5,"coil-db/HN-parm.dat"); 
      break;

    case HA: 
      strcpy(n1,"coil-db/HA-uno.dat"); 
      strcpy(n2,"coil-db/HA-due.dat"); 
      strcpy(n3,"coil-db/HA-tre.dat"); 
      strcpy(n4,"coil-db/HA-quattro.dat"); 
      strcpy(n5,"coil-db/HA-parm.dat"); 
      break;

    case NH: 
      strcpy(n1,"coil-db/NH-uno.dat"); 
      strcpy(n2,"coil-db/NH-due.dat"); 
      strcpy(n3,"coil-db/NH-tre.dat"); 
      strcpy(n4,"coil-db/NH-quattro.dat"); 
      strcpy(n5,"coil-db/NH-parm.dat"); 
      break;
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

  fscanf(uno, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  fscanf(due, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  fscanf(tre, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  fscanf(qua, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
 
  format = "%*s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf";
  for(i=0;i<21;i++)
  {
    fscanf(m, "%*s %lf %*f", &mean[i]);

    fscanf(uno, format, &s1[i][0],  &s1[i][1],  &s1[i][2],  &s1[i][3],  &s1[i][4],  &s1[i][5],  &s1[i][6],  &s1[i][7],  &s1[i][8],  &s1[i][9], 
                        &s1[i][10], &s1[i][11], &s1[i][12], &s1[i][13], &s1[i][14], &s1[i][15], &s1[i][16], &s1[i][17], &s1[i][18], &s1[i][19], 
                        &s1[i][20]);

    fscanf(due, format, &s2[i][0],  &s2[i][1],  &s2[i][2],  &s2[i][3],  &s2[i][4],  &s2[i][5],  &s2[i][6],  &s2[i][7],  &s2[i][8],  &s2[i][9], 
                        &s2[i][10], &s2[i][11], &s2[i][12], &s2[i][13], &s2[i][14], &s2[i][15], &s2[i][16], &s2[i][17], &s2[i][18], &s2[i][19], 
                        &s2[i][20]);

    fscanf(tre, format, &d1[i][0],  &d1[i][1],  &d1[i][2],  &d1[i][3],  &d1[i][4],  &d1[i][5],  &d1[i][6],  &d1[i][7],  &d1[i][8],  &d1[i][9], 
                        &d1[i][10], &d1[i][11], &d1[i][12], &d1[i][13], &d1[i][14], &d1[i][15], &d1[i][16], &d1[i][17], &d1[i][18], &d1[i][19], 
                        &d1[i][20]);

    fscanf(qua, format, &d2[i][0],  &d2[i][1],  &d2[i][2],  &d2[i][3],  &d2[i][4],  &d2[i][5],  &d2[i][6],  &d2[i][7],  &d2[i][8],  &d2[i][9], 
                        &d2[i][10], &d2[i][11], &d2[i][12], &d2[i][13], &d2[i][14], &d2[i][15], &d2[i][16], &d2[i][17], &d2[i][18], &d2[i][19], 
                        &d2[i][20]);
  }

  fclose(m);
  fclose(uno);
  fclose(due);
  fclose(tre);
  fclose(qua);

  fscanf(par, "%lf", &as1);
  fscanf(par, "%lf", &as2);
  fscanf(par, "%lf", &ad1);
  fscanf(par, "%lf", &ad2);
  
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
  }

  return 0;
}

int do_ppii(int ashift, double *tcs, int *nseq, int sql, int ph, char *dbpath)
{
  char   n1[25], n2[25], n3[25], n4[25], n5[25], *format, str[MAXLENGTH];
  double mean[21];
  double s1[21][21], s2[21][21], d1[21][21], d2[21][21];
  double as1, as2, ad1, ad2;
  int    i;
  FILE   *m=NULL, *uno, *due, *tre, *qua, *par;

  switch(ashift)
  {
    case CA: 
      strcpy(n1,"ppii-db/CA-uno.dat"); 
      strcpy(n2,"ppii-db/CA-due.dat"); 
      strcpy(n3,"ppii-db/CA-tre.dat"); 
      strcpy(n4,"ppii-db/CA-quattro.dat"); 
      strcpy(n5,"ppii-db/CA-parm.dat"); 
      break;

    case CB: 
      strcpy(n1,"ppii-db/CB-uno.dat"); 
      strcpy(n2,"ppii-db/CB-due.dat"); 
      strcpy(n3,"ppii-db/CB-tre.dat"); 
      strcpy(n4,"ppii-db/CB-quattro.dat"); 
      strcpy(n5,"ppii-db/CB-parm.dat"); 
      break;

    case CO: 
      strcpy(n1,"ppii-db/CO-uno.dat"); 
      strcpy(n2,"ppii-db/CO-due.dat"); 
      strcpy(n3,"ppii-db/CO-tre.dat"); 
      strcpy(n4,"ppii-db/CO-quattro.dat"); 
      strcpy(n5,"ppii-db/CO-parm.dat"); 
      break;

    case HN: 
      strcpy(n1,"ppii-db/HN-uno.dat"); 
      strcpy(n2,"ppii-db/HN-due.dat"); 
      strcpy(n3,"ppii-db/HN-tre.dat"); 
      strcpy(n4,"ppii-db/HN-quattro.dat"); 
      strcpy(n5,"ppii-db/HN-parm.dat"); 
      break;

    case HA: 
      strcpy(n1,"ppii-db/HA-uno.dat"); 
      strcpy(n2,"ppii-db/HA-due.dat"); 
      strcpy(n3,"ppii-db/HA-tre.dat"); 
      strcpy(n4,"ppii-db/HA-quattro.dat"); 
      strcpy(n5,"ppii-db/HA-parm.dat"); 
      break;

    case NH: 
      strcpy(n1,"ppii-db/NH-uno.dat"); 
      strcpy(n2,"ppii-db/NH-due.dat"); 
      strcpy(n3,"ppii-db/NH-tre.dat"); 
      strcpy(n4,"ppii-db/NH-quattro.dat"); 
      strcpy(n5,"ppii-db/NH-parm.dat"); 
      break;
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

  fscanf(uno, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  fscanf(due, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  fscanf(tre, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
  fscanf(qua, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
 
  format = "%*s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf";
  for(i=0;i<21;i++)
  {
    fscanf(m, "%*s %lf %*f", &mean[i]);

    fscanf(uno, format, &s1[i][0],  &s1[i][1],  &s1[i][2],  &s1[i][3],  &s1[i][4],  &s1[i][5],  &s1[i][6],  &s1[i][7],  &s1[i][8],  &s1[i][9], 
                        &s1[i][10], &s1[i][11], &s1[i][12], &s1[i][13], &s1[i][14], &s1[i][15], &s1[i][16], &s1[i][17], &s1[i][18], &s1[i][19], 
                        &s1[i][20]);

    fscanf(due, format, &s2[i][0],  &s2[i][1],  &s2[i][2],  &s2[i][3],  &s2[i][4],  &s2[i][5],  &s2[i][6],  &s2[i][7],  &s2[i][8],  &s2[i][9], 
                        &s2[i][10], &s2[i][11], &s2[i][12], &s2[i][13], &s2[i][14], &s2[i][15], &s2[i][16], &s2[i][17], &s2[i][18], &s2[i][19], 
                        &s2[i][20]);

    fscanf(tre, format, &d1[i][0],  &d1[i][1],  &d1[i][2],  &d1[i][3],  &d1[i][4],  &d1[i][5],  &d1[i][6],  &d1[i][7],  &d1[i][8],  &d1[i][9], 
                        &d1[i][10], &d1[i][11], &d1[i][12], &d1[i][13], &d1[i][14], &d1[i][15], &d1[i][16], &d1[i][17], &d1[i][18], &d1[i][19], 
                        &d1[i][20]);

    fscanf(qua, format, &d2[i][0],  &d2[i][1],  &d2[i][2],  &d2[i][3],  &d2[i][4],  &d2[i][5],  &d2[i][6],  &d2[i][7],  &d2[i][8],  &d2[i][9], 
                        &d2[i][10], &d2[i][11], &d2[i][12], &d2[i][13], &d2[i][14], &d2[i][15], &d2[i][16], &d2[i][17], &d2[i][18], &d2[i][19], 
                        &d2[i][20]);
  }

  fclose(m);
  fclose(uno);
  fclose(due);
  fclose(tre);
  fclose(qua);

  fscanf(par, "%lf", &as1);
  fscanf(par, "%lf", &as2);
  fscanf(par, "%lf", &ad1);
  fscanf(par, "%lf", &ad2);
  
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
  }

  return 0;
}

double do_predict(char *bmrb, int ph, char *dbpath, char *out, int firstres, char *argv[], int argc)
{
  int         ashift, bshift, a, b;
  int         i, j, result, err, ti;
  static int  sql=0;
  double      tmp_fix[4], sig_ta, sig_tb, prec, rerr=0., urange, lrange;
  FILE        *ref, *tab=NULL, *cstable=NULL;
  static char seq[MAXLENGTH];
  char        *format=NULL, ss[MAXLENGTH], css[MAXLENGTH], str[MAXLENGTH];
  char        *line=NULL, namecstable[256];
  size_t      linecap = 0;
  ssize_t     linelen;
  /* variables for multivariate gaussian, size is predefined for performance reasons */
  double      vec3h[3], vec4h[4], vec5h[5], vec6h[6], tmp3[3];
  double      vec3b[3], vec4b[4], vec5b[5], vec6b[6], tmp4[4];
  double      vec3c[3], vec4c[4], vec5c[5], vec6c[6], tmp5[5];
  double      vec3p[3], vec4p[4], vec5p[5], vec6p[6], tmp6[6];
  double      mat3[3][3], mat4[4][4], mat5[5][5], mat6[6][6];
  double      imat3[3][3], imat4[4][4], imat5[5][5], imat6[6][6];
  double      norm, vhMvh=0., vbMvb=0., vcMvc=0., vpMvp=0.;
#ifndef LIB
  int         mtot, *star, *nseq;
  double      sig_a[21][6], sig_b[21][6], sig_c[21][6], sig_p[21][6];
  double      wc[5][6], cor[6][6];
  double      **ecs, **intp, **pp, **spp;
  double      *alphacs[6], *betacs[6], *coilcs[6], *ppiics[6];
  double      *Aprob, *Bprob, *Cprob, *Pprob;
  double      mhelix, mbeta, mcoil, mppii;
  char        n7[256], backup[256];
  FILE        *pred;
#else
  static int    *star, *nseq;
  static double *alphacs[6], *betacs[6], *coilcs[6], *ppiics[6];
  static double *Aprob, *Bprob, *Cprob, *Pprob;
  static double **pp, **spp;
#endif

  /* allocate everything */
#ifdef LIB
  if(!count) {
#endif
    nseq = (int *) calloc (MAXLENGTH,sizeof(int));
    star = (int *) calloc (MAXLENGTH,sizeof(int));
    ecs =  (double **) calloc(NSHIFT,sizeof(double*));
    for(i=0;i<NSHIFT;i++) {
      alphacs[i]= (double *) calloc(MAXLENGTH,sizeof(double));
      betacs[i] = (double *) calloc(MAXLENGTH,sizeof(double));
      coilcs[i] = (double *) calloc(MAXLENGTH,sizeof(double));
      ppiics[i] = (double *) calloc(MAXLENGTH,sizeof(double));
      ecs[i]    = (double *) calloc(MAXLENGTH,sizeof(double));
    }
    intp = (double **) calloc(MAXLENGTH,sizeof(double*));
    Aprob = (double *) calloc(MAXLENGTH,sizeof(double));
    Bprob = (double *) calloc(MAXLENGTH,sizeof(double));
    Cprob = (double *) calloc(MAXLENGTH,sizeof(double));
    Pprob = (double *) calloc(MAXLENGTH,sizeof(double));
    pp    = (double **) calloc(MAXLENGTH,sizeof(double*));
    spp   = (double **) calloc(MAXLENGTH,sizeof(double*));
    for(i=0;i<MAXLENGTH;i++) {
      intp[i] = (double *) calloc(4,sizeof(double*));
      pp[i]   = (double *) calloc(4,sizeof(double));
      spp[i]  = (double *) calloc(4,sizeof(double));
    }
#ifdef LIB
  }
#endif

#ifndef LIB
  for(j=0;j<MAXLENGTH;j++) {seq[j]=0; css[j]=0; ss[j]=0; star[j]=0;}

  /* output file */
  strcpy(n7,"SS-results.dat");
  if(out!=NULL) strcpy(n7,out);
  sprintf(backup, "%s.old", n7);
  rename(n7, backup);

  /* read weights file */
  sprintf(str, "%s/other-db/weights.tab", dbpath);
  tab = fopen(str,"r");
  if (tab==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}
  err = fscanf(tab, "%*s %*s %*s %*s %*s %*s %*s");
  for(i=0;i<5;i++)
    err = fscanf(tab, "%*s %lf %lf %lf %lf %lf %lf", &wc[i][0], &wc[i][1], &wc[i][2], &wc[i][3], &wc[i][4], &wc[i][5]);
  fclose(tab);

  /* read correlations file */ 
  sprintf(str, "%s/other-db/correlations.tab", dbpath);
  tab = fopen(str,"r");
  if (tab==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}
  for(i=0;i<6;i++)
    err = fscanf(tab, "%lf %lf %lf %lf %lf %lf", &cor[i][0], &cor[i][1], &cor[i][2], &cor[i][3], &cor[i][4], &cor[i][5]);
  fclose(tab);

  /* read sigma helices */
  sprintf(str, "%s/other-db/helix.err", dbpath);
  tab = fopen(str,"r");
  if (tab==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}
  for(i=0;i<21;i++)
  {
    err = fscanf(tab, "%*s %lf %lf %lf %lf %lf %lf", &sig_a[i][0], &sig_a[i][1], &sig_a[i][2], &sig_a[i][3], &sig_a[i][4], &sig_a[i][5]);
  }
  fclose(tab);

  /* read sigma betas */
  sprintf(str, "%s/other-db/beta.err", dbpath);
  tab = fopen(str,"r");
  if (tab==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}
  for(i=0;i<21;i++)
  {
    err = fscanf(tab, "%*s %lf %lf %lf %lf %lf %lf", &sig_b[i][0], &sig_b[i][1], &sig_b[i][2], &sig_b[i][3], &sig_b[i][4], &sig_b[i][5]);
  }
  fclose(tab);

  /* read sigma coils */
  sprintf(str, "%s/other-db/coil.err", dbpath);
  tab = fopen(str,"r");
  if (tab==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}
  for(i=0;i<21;i++)
  {
    err = fscanf(tab, "%*s %lf %lf %lf %lf %lf %lf", &sig_c[i][0], &sig_c[i][1], &sig_c[i][2], &sig_c[i][3], &sig_c[i][4], &sig_c[i][5]);
  }
  fclose(tab);

  /* read sigma ppiis */
  sprintf(str, "%s/other-db/ppii.err", dbpath);
  tab = fopen(str,"r");
  if (tab==NULL) {fprintf(stderr, "%s file not found!\n", str); return 1;}
  for(i=0;i<21;i++)
  {
    err = fscanf(tab, "%*s %lf %lf %lf %lf %lf %lf", &sig_p[i][0], &sig_p[i][1], &sig_p[i][2], &sig_p[i][3], &sig_p[i][4], &sig_p[i][5]);
  }
  fclose(tab);
#endif

#ifdef LIB
  if(!count) {
#endif
    if(shifty) format = "%i %c %lf %lf %lf %lf %lf %lf"; 
    else       format = "%i %*i %c %*c %c %lf %lf %lf %lf %lf %lf";
    /* read the input file */ 
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
#else
  fprintf(stderr,"Protein length is %i. Last residue is %c\n", sql, seq[sql-1]); fflush(stderr);
#endif

#ifdef LIB
  if(!count){
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
        default: fprintf(stderr, "ERROR: Aminoacid %c not recognized! line %d\n", seq[i], i); return 1; break;
      }
    }

    /* predict all the chemical shifts */
    for(i=0;i<MAXLENGTH;i++) {
      for(ashift=0;ashift<NSHIFT;ashift++) { 
        alphacs[ashift][i]=0.; 
        betacs[ashift][i]=0.; 
        coilcs[ashift][i]=0.; 
        ppiics[ashift][i]=0.; 
      }
    }
    for(ashift=0;ashift<NSHIFT;ashift++) {
      result = do_alpha(ashift, alphacs[ashift], nseq, sql, ph, dbpath);
      result = do_coil(ashift, coilcs[ashift], nseq, sql, ph, dbpath);
      result = do_ppii(ashift, ppiics[ashift], nseq, sql, ph, dbpath);
      result = do_beta(ashift, betacs[ashift], nseq, sql, ph, dbpath);

      if(debug) { /* write the tables with the predicted chemical shifts */
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
          fprintf(cstable, "%i %c %6.2f %6.2f %6.2f %6.2f\n", i, seq[i], alphacs[ashift][i], coilcs[ashift][i], 
                                                                         ppiics[ashift][i], betacs[ashift][i]);
        fclose(cstable);
      }

      /* fix the predictions sorting them in the expected order */
      for(i=0;i<sql;i++) {
        if(ashift==CA||ashift==CO) {
          /* the correct order is Helix > Coil > PPII > beta */
          tmp_fix[3]=alphacs[ashift][i]; tmp_fix[2]=coilcs[ashift][i]; tmp_fix[1]=ppiics[ashift][i]; tmp_fix[0]=betacs[ashift][i];
          qsort(tmp_fix, sizeof(tmp_fix)/sizeof(tmp_fix[0]), sizeof(tmp_fix[0]), cmp);
          alphacs[ashift][i] = tmp_fix[3]; coilcs[ashift][i]=tmp_fix[2]; ppiics[ashift][i]=tmp_fix[1]; betacs[ashift][i]=tmp_fix[0];
        } else if(ashift==HA||ashift==CB||ashift==NH||ashift==HN) {
          /* the correct order is Beta > PPII > Coil > Helix */
          tmp_fix[0]=alphacs[ashift][i]; tmp_fix[1]=coilcs[ashift][i]; tmp_fix[2]=ppiics[ashift][i]; tmp_fix[3]=betacs[ashift][i];
          qsort(tmp_fix, sizeof(tmp_fix)/sizeof(tmp_fix[0]), sizeof(tmp_fix[0]), cmp);
          alphacs[ashift][i] = tmp_fix[0]; coilcs[ashift][i]=tmp_fix[1]; ppiics[ashift][i]=tmp_fix[2]; betacs[ashift][i]=tmp_fix[3];
        }
      }
    }
#ifdef LIB
  }
#endif

  /* Here we have all the calculated chemical shifts for the four secondary structures for each residue
   * it is now possible to cycle over the residues and the chemical shifts to calculate the populations */
  for(i=0;i<sql;i++) {
    Aprob[i] = 999.999;
    Bprob[i] = 999.999;
    Cprob[i] = 999.999;
    Pprob[i] = 999.999;
    star[i]  = 0;
    if(nseq[i]==Z||(i>0&&nseq[i-1]==Z)||(i+1<sql&&nseq[i+1]==Z)) continue;

    /* define the dimensionality for the multivariate function (from 3 to 6) */
    for(ashift=0;ashift<NSHIFT;ashift++) {
      if(nseq[i]==G&&(ashift==HA||ashift==CB)) continue;
      if(nseq[i]==P&&(ashift==HN||ashift==NH)) continue;
      /* if we have an experimental chemical shift for a nucleus-residue */
      if(ecs[ashift][i]<900.&&ecs[ashift][i]>0.) {
        /* and the chemical shift is within a reasonable range of values */
        switch(ashift) {
          /* here helix > beta */
          case CA:
          case CO:
            urange = alphacs[ashift][i]+NSIGMA*sig_a[nseq[i]][ashift];
            lrange = betacs[ashift][i] -NSIGMA*sig_b[nseq[i]][ashift];
            if(ecs[ashift][i]<urange && ecs[ashift][i]>lrange) star[i]++;
            else {
#ifndef LIB
              fprintf(stderr, "WARNING: %4i%c %2s (%6.2f) is out of range (%6.2f - %6.2f)\n", 
                      i+firstres, seq[i], shtostring(ashift), ecs[ashift][i], lrange, urange);
#endif
              /* this is very important otherwise the chemical shifts will be used in place of a good one */
              ecs[ashift][i]=999.999;
            }
            break;
          /* here helix < beta */
          case CB:
          case HA:
          case HN:
          case NH:
            urange = betacs[ashift][i] +NSIGMA*sig_b[nseq[i]][ashift];
            lrange = alphacs[ashift][i]-NSIGMA*sig_a[nseq[i]][ashift];
            if(ecs[ashift][i]<urange && ecs[ashift][i]>lrange) star[i]++;
            else {
#ifndef LIB
              fprintf(stderr, "WARNING: %4i%c %2s (%6.2f) is out of range (%6.2f - %6.2f)\n", 
                      i+firstres, seq[i], shtostring(ashift), ecs[ashift][i], lrange, urange);
#endif
              /* this is very important otherwise the chemical shifts will be used in place of a good one */
              ecs[ashift][i]=999.999;
            }
            break;
        }
      }
    }
    if(debug) fprintf(stdout, "res %d type %d star %d\n", i, nseq[i], star[i]); 
   
    /* generate the vector and matrix needed */
    a=0;
    for(ashift=0;ashift<NSHIFT;ashift++) {
      if(nseq[i]==G&&(ashift==HA||ashift==CB)) continue;
      if(nseq[i]==P&&(ashift==HN||ashift==NH)) continue;
      if(ecs[ashift][i]<900.&&ecs[ashift][i]>0.) {
      sig_ta = sqrt(0.25*(sig_a[nseq[i]][ashift]*sig_a[nseq[i]][ashift]+
                          sig_b[nseq[i]][ashift]*sig_b[nseq[i]][ashift]+
                          sig_c[nseq[i]][ashift]*sig_c[nseq[i]][ashift]+
                          sig_p[nseq[i]][ashift]*sig_p[nseq[i]][ashift]));

        b=0;
        for(bshift=0;bshift<NSHIFT;bshift++) {
          if(nseq[i]==G&&(bshift==HA||bshift==CB)) continue;
          if(nseq[i]==P&&(bshift==HN||bshift==NH)) continue;
          if(ecs[bshift][i]<900.&&ecs[bshift][i]>0.) {
            sig_tb = sqrt(0.25*(sig_a[nseq[i]][bshift]*sig_a[nseq[i]][bshift]+
                                sig_b[nseq[i]][bshift]*sig_b[nseq[i]][bshift]+
                                sig_c[nseq[i]][bshift]*sig_c[nseq[i]][bshift]+
                                sig_p[nseq[i]][bshift]*sig_p[nseq[i]][bshift]));
            switch(star[i]) {
              case 3: 
                mat3[a][b] = cor[ashift][bshift]*sig_ta*sig_tb; 
                break;
              case 4: 
                mat4[a][b] = cor[ashift][bshift]*sig_ta*sig_tb; 
                break;
              case 5: 
                mat5[a][b] = cor[ashift][bshift]*sig_ta*sig_tb; 
                break;
              case 6: 
                mat6[a][b] = cor[ashift][bshift]*sig_ta*sig_tb; 
                break;
              default: 
                break;
            }
            b++;
          }
        }

        switch(star[i]) {
          case 3: 
            vec3h[a] = (ecs[ashift][i]-alphacs[ashift][i]); 
            vec3b[a] = (ecs[ashift][i]-betacs[ashift][i]); 
            vec3c[a] = (ecs[ashift][i]-coilcs[ashift][i]); 
            vec3p[a] = (ecs[ashift][i]-ppiics[ashift][i]);
            break;
          case 4: 
            vec4h[a] = (ecs[ashift][i]-alphacs[ashift][i]); 
            vec4b[a] = (ecs[ashift][i]-betacs[ashift][i]); 
            vec4c[a] = (ecs[ashift][i]-coilcs[ashift][i]); 
            vec4p[a] = (ecs[ashift][i]-ppiics[ashift][i]); 
            break;
          case 5: 
            vec5h[a] = (ecs[ashift][i]-alphacs[ashift][i]); 
            vec5b[a] = (ecs[ashift][i]-betacs[ashift][i]); 
            vec5c[a] = (ecs[ashift][i]-coilcs[ashift][i]); 
            vec5p[a] = (ecs[ashift][i]-ppiics[ashift][i]); 
            break;
          case 6: 
            vec6h[a] = (ecs[ashift][i]-alphacs[ashift][i]); 
            vec6b[a] = (ecs[ashift][i]-betacs[ashift][i]); 
            vec6c[a] = (ecs[ashift][i]-coilcs[ashift][i]); 
            vec6p[a] = (ecs[ashift][i]-ppiics[ashift][i]); 
            break;
          default: 
            break;
         }
         a++;  
       }
    }

    if(debug) {
      fprintf(stdout, "Covariance Matrix:\n");
      switch(star[i]) {
        case 3: 
          for(ashift=0;ashift<star[i];ashift++) {
            for(bshift=0;bshift<star[i];bshift++) {
              fprintf(stdout, "%9.6f ", mat3[ashift][bshift]);
            }
            fprintf(stdout, "\n");
          }
          break;
        case 4: 
          for(ashift=0;ashift<star[i];ashift++) {
            for(bshift=0;bshift<star[i];bshift++) {
              fprintf(stdout, "%9.6f ", mat4[ashift][bshift]);
            }
            fprintf(stdout, "\n");
          }
          break;
        case 5: 
          for(ashift=0;ashift<star[i];ashift++) {
            for(bshift=0;bshift<star[i];bshift++) {
              fprintf(stdout, "%9.6f ", mat5[ashift][bshift]);
            }
            fprintf(stdout, "\n");
          }
          break;
        case 6: 
          for(ashift=0;ashift<star[i];ashift++) {
            for(bshift=0;bshift<star[i];bshift++) {
              fprintf(stdout, "%9.6f ", mat6[ashift][bshift]);
            }
            fprintf(stdout, "\n");
          }
          break;
        default:
          break;
      }
    } 
    
    switch(star[i]) {
      case 3:
        m_inv3(mat3,imat3);               /* calculate the inverse of the matrix */
        m_v_m3(imat3, vec3h, tmp3);       /* calculate the mat*vec product */
        vhMvh = v_dot_p3(vec3h, tmp3);    /* calculate the dot product*/
        m_v_m3(imat3, vec3b, tmp3);
        vbMvb = v_dot_p3(vec3b, tmp3);
        m_v_m3(imat3, vec3c, tmp3);
        vcMvc = v_dot_p3(vec3c, tmp3);
        m_v_m3(imat3, vec3p, tmp3);
        vpMvp = v_dot_p3(vec3p, tmp3);
        if(debug) {
          fprintf(stdout, "Inverted Matrix:\n");
          for(ashift=0;ashift<3;ashift++) {
            for(bshift=0;bshift<3;bshift++) {
              fprintf(stdout, "%9.6f ", imat3[ashift][bshift]);
            }
            fprintf(stdout, "\n");
          }
          fprintf(stdout, "vh: %f\nvb: %f\nvc: %f\nvp: %f\n", vhMvh, vbMvb, vcMvc, vpMvp);
        }
        break;
      case 4:
        m_inv4(mat4,imat4);               /* calculate the inverse of the matrix */
        m_v_m4(imat4, vec4h, tmp4);       /* calculate the mat*vec product */
        vhMvh = v_dot_p4(vec4h, tmp4);    /* calculate the dot product*/
        m_v_m4(imat4, vec4b, tmp4);
        vbMvb = v_dot_p4(vec4b, tmp4);
        m_v_m4(imat4, vec4c, tmp4);
        vcMvc = v_dot_p4(vec4c, tmp4);
        m_v_m4(imat4, vec4p, tmp4);
        vpMvp = v_dot_p4(vec4p, tmp4);
        if(debug) {
          fprintf(stdout, "Inverted Matrix:\n");
          for(ashift=0;ashift<4;ashift++) {
            for(bshift=0;bshift<4;bshift++) {
              fprintf(stdout, "%9.6f ", imat4[ashift][bshift]);
            }
            fprintf(stdout, "\n");
          }
          fprintf(stdout, "vh: %f\nvb: %f\nvc: %f\nvp: %f\n", vhMvh, vbMvb, vcMvc, vpMvp);
        }
        break;
      case 5:
        m_inv5(mat5,imat5);               /* calculate the inverse of the matrix */
        m_v_m5(imat5, vec5h, tmp5);       /* calculate the mat*vec product */
        vhMvh = v_dot_p5(vec5h, tmp5);    /* calculate the dot product*/
        m_v_m5(imat5, vec5b, tmp5);
        vbMvb = v_dot_p5(vec5b, tmp5);
        m_v_m5(imat5, vec5c, tmp5);
        vcMvc = v_dot_p5(vec5c, tmp5);
        m_v_m5(imat5, vec5p, tmp5);
        vpMvp = v_dot_p5(vec5p, tmp5);
        if(debug) {
          fprintf(stdout, "Inverted Matrix:\n");
          for(ashift=0;ashift<5;ashift++) {
            for(bshift=0;bshift<5;bshift++) {
              fprintf(stdout, "%9.6f ", imat5[ashift][bshift]);
            }
            fprintf(stdout, "\n");
          }
          fprintf(stdout, "vh: %f\nvb: %f\nvc: %f\nvp: %f\n", vhMvh, vbMvb, vcMvc, vpMvp);
        }
        break;
      case 6:
        m_inv6(mat6,imat6);               /* calculate the inverse of the matrix */
        m_v_m6(imat6, vec6h, tmp6);       /* calculate the mat*vec product */
        vhMvh = v_dot_p6(vec6h, tmp6);    /* calculate the dot product*/
        m_v_m6(imat6, vec6b, tmp6);
        vbMvb = v_dot_p6(vec6b, tmp6);
        m_v_m6(imat6, vec6c, tmp6);
        vcMvc = v_dot_p6(vec6c, tmp6);
        m_v_m6(imat6, vec6p, tmp6);
        vpMvp = v_dot_p6(vec6p, tmp6);
        if(debug) {
          fprintf(stdout, "Inverted Matrix:\n");
          for(ashift=0;ashift<6;ashift++) {
            for(bshift=0;bshift<6;bshift++) {
              fprintf(stdout, "%9.6f ", imat6[ashift][bshift]);
            }
            fprintf(stdout, "\n");
          }
          fprintf(stdout, "vh: %f\nvb: %f\nvc: %f\nvp: %f\n", vhMvh, vbMvb, vcMvc, vpMvp);
        }
        break;
      default: 
        break;
    }

    if(star[i]>3) {
      /* calculate the gaussians normalisation (no det))*/ 
      norm = 1./sqrt(pow(2.*M_PI,star[i])); 
      /* the following numbers say what is the most probable sec struc, in the 
       * following these numbers are combined to give the relative populations */
      Aprob[i] = norm*exp(-0.5*vhMvh);
      Bprob[i] = norm*exp(-0.5*vbMvb);
      Cprob[i] = norm*exp(-0.5*vcMvc);
      Pprob[i] = norm*exp(-0.5*vpMvp);
      /* the wc relative weigths are optimised for the number of available chemical shifts */
      /* this is done so to make more reliable the predictions with less than 6 nuclei */
      if(nseq[i]!=P&&nseq[i]!=G) {
        pp[i][0] = pow(Aprob[i],wc[0][6-star[i]]);  
        pp[i][1] = pow(Bprob[i],wc[1][6-star[i]]);  
        pp[i][2] = pow(Cprob[i],wc[2][6-star[i]]);  
        pp[i][3] = pow(Pprob[i],wc[3][6-star[i]]);
      } else {
        pp[i][0] = pow(Aprob[i],wc[0][4-star[i]]);  
        pp[i][1] = pow(Bprob[i],wc[1][4-star[i]]);  
        pp[i][2] = pow(Cprob[i],wc[2][4-star[i]]);  
        pp[i][3] = pow(Pprob[i],wc[3][4-star[i]]);
      }
      norm = pp[i][0] + pp[i][1] + pp[i][2] + pp[i][3];  
      pp[i][0] /= norm;  
      pp[i][1] /= norm;  
      pp[i][2] /= norm;  
      pp[i][3] /= norm;  
    }
  }

  /* smoothing window */
  star[0]=star[sql-1]=0;
  do_averaging(pp, spp, nseq, sql, wc, star, ss);

  /* quality check */
  if(!shifty) {
    j=0; prec=0;
    for(i=0;i<sql;i++)
    {
      if(star[i]>3){
        if(css[i]!='T'&&css[i]!='C'&&css[i]!='S'&&
           css[i]!='P'&&css[i]!='H'&&css[i]!='I'&&
           css[i]!='G'&&css[i]!='B'&&css[i]!='E') continue;

        if(css[i]=='T'||css[i]=='S') css[i]='C';
        if(css[i]!=ss[i]) prec++;
   
        if(ss[i]=='C'&&css[i]=='P'&&spp[i][0]<spp[i][3]&&spp[i][1]<spp[i][3]) prec-=0.75;
        if(ss[i]=='H'&&css[i]=='I') prec--;
        if(ss[i]=='H'&&css[i]=='G') prec--;
        if(ss[i]=='C'&&css[i]=='G') prec-=0.50;
        if(ss[i]=='E'&&css[i]=='B') prec-=0.50;
        if(ss[i]=='C'&&css[i]=='B') prec--;

        j++;
      }
    }

    rerr = prec/(double)j*100.; 
  }

#ifndef LIB
  fprintf(stderr,"SQ:%s\n", seq);
  fprintf(stderr,"SS:%s\n", ss);
  if(!shifty) {  
    fprintf(stderr,"SS:%s\n", css);
    fprintf(stdout,"Err %7.3f\n", rerr); 
  }

  /* Calculating total populations */
  mhelix = mbeta = mcoil = mppii = 0.;
  mtot = 0;
  for(i=0;i<sql;i++)
  {
    if(spp[i][0]!=0||spp[i][1]!=0||spp[i][2]!=0) {
      mhelix += spp[i][0];
      mbeta  += spp[i][1];
      mcoil  += spp[i][2];
      mppii  += spp[i][3];
      mtot++;
    }
  }
  mhelix = mhelix/mtot*100.;
  mbeta  = mbeta/mtot*100.;
  mcoil  = mcoil/mtot*100.;
  mppii  = mppii/mtot*100.;

  /* Writing the output file */
  pred = fopen(n7, "w");
  fprintf(pred,"#d2D - v. 2.1.0\n#(c) Carlo Camilloni\n#PLEASE CITE:\n");
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
  fprintf(pred,"#Polyproline II (PPII)(P): %4.1f%%\n", mppii);
  fprintf(pred,"#Coil(C): %4.1f%%\n\n", mcoil);
  fprintf(pred,"#Populations per residue (residues marked with a * are less reliable):\n");
  fprintf(pred,"#num \t res \t      Helix      Beta       Coil       PPII  SS \n");
  for(i=0;i<sql;i++)
  {
    if((pp[i][0]!=0||pp[i][1]!=0||pp[i][2]!=0)&&star[i]>3) {
      fprintf(pred, "%i \t %c \t %10.6f %10.6f %10.6f %10.6f %c\n", i+firstres, seq[i], spp[i][0], spp[i][1], spp[i][2], spp[i][3], ss[i]);
    } else if(star[i]<4&&ss[i]!=' ') {
      fprintf(pred, "%i \t %c \t %10.3f %10.3f %10.3f %10.3f %c*\n", i+firstres, seq[i], spp[i][0], spp[i][1], spp[i][2], spp[i][3], ss[i]);
    } else {
      if(!dbformat) fprintf(pred, "#%i \t %c \t \n", i+firstres, seq[i]);
      else fprintf(pred, "@%i \t %c \t \n", i+firstres, seq[i]);
    }
  }
  fprintf(pred,"#DONE!\n");
  fprintf(pred,"#d2D - v. 2.1.0\n#(c) Carlo Camilloni\n#PLEASE CITE:\n");
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
    free(Aprob);
    free(Bprob);
    free(Cprob);
    free(Pprob);
    for(i=0;i<NSHIFT;i++) {
      free(ecs[i]);
      free(alphacs[i]);
      free(betacs[i]);
      free(coilcs[i]);
      free(ppiics[i]);
    }
    free(ecs);
    free(intp);
    for(i=0;i<MAXLENGTH;i++) {
      free(intp[i]);
      free(pp[i]);
      free(spp[i]);
    }
    free(pp);
    free(spp); 
#ifdef LIB
  }
#endif

  return rerr;
}

#ifndef LIB
void help()
{
  fprintf(stderr,"\nd2D - v. 2.1.0\n(c) Carlo Camilloni\n\nPLEASE CITE:\n");
  fprintf(stderr,"Camilloni C., De Simone A., Vranken W., and Vendruscolo M.\n");
  fprintf(stderr,"Determination of Secondary Structure Populations in Disordered States of Proteins\n");
  fprintf(stderr,"using NMR Chemical Shifts\n");
  fprintf(stderr,"Biochemistry 2012, 51: 2224-2231\n\n");
  fprintf(stderr,"\t-pH        \t([neutral], acid)\n");
  fprintf(stderr,"\t-file      \t(input file name)\n");
  fprintf(stderr,"\t-out       \t(output file name)\n");
  fprintf(stderr,"\t-fres      \t(rescale residue numbers from fres)\n");
  fprintf(stderr,"\t-shifty    \t(input file in the shifty format (website format))\n");
  fprintf(stderr,"\t-debug     \t(more verbose)\n");
  fprintf(stderr,"\t-help      \t(here we are!)\n\n");
  fprintf(stderr,"Use single letter aminoacid format, C for reduced cysteine X for oxydaised\n\n");
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
  fprintf(stderr,"EXAMPLE: (internal format)\n");
  fprintf(stderr,"1  1	P X C	62.69	32.89	999.99	999.9 4.32 999.99\n");
  fprintf(stderr,"2  2	N X C	52.39	38.99	174.40	999.9 5.03 999.99\n");
  fprintf(stderr,"3  3	F X C	60.09	31.99	176.00	9.82  4.32 121.78\n");
  fprintf(stderr,"4  4	S X E	60.39	64.19	174.30	8.49  4.30 112.38\n");
  fprintf(stderr,"5  5	G X E	44.59	999.9	170.40	9.22  3.96 110.48\n");
  fprintf(stderr,"6  6	N X E	52.69	39.89	174.00	8.15  5.52 118.48\n");
  fprintf(stderr,"7  7	E X C	56.39	31.89	175.00	9.35  5.14 123.18\n");
  fprintf(stderr,"8  8	K X C	53.79	35.79	174.90	10.17 5.28 123.08\n");
  fprintf(stderr,"9  9	I X C	55.19	33.29	175.80	9.05  3.92 125.48\n\n");
}

int main(int argc, char *argv[])
{
  int    status=0, ph=0, firstres=-999;
  int    i, j=-1, o=-1;
  double dst=0.;
  char   *dbpath;

  dbpath = getenv("CAMDBV3");
  if(dbpath==NULL) {fprintf (stderr, "The CAMDBPATH is not set!\n"); return 1;}
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
    if(!strcmp(argv[i],"-fres")) firstres = atoi(argv[i+1]);
    if(!strcmp(argv[i],"-shifty")) shifty=1;
    if(!strcmp(argv[i],"-dbformat")) dbformat=1;
    if(!strcmp(argv[i],"-help")) {help(); return 0;}
  }

  if(debug) {fprintf(stdout, "In main: CAMDBPATH is %s\n", dbpath); fflush(stdout);}
  if(j==-1) {help(); return 0;}

    
  if(argv[j]!=NULL&&o!=-1) dst = do_predict(argv[j], ph, dbpath, argv[o], firstres, argv, argc);
  else if(argv[j]!=NULL&&o==-1) dst = do_predict(argv[j], ph, dbpath, NULL, firstres, argv, argc);
  else {fprintf(stderr, "Wrong file name: %s, expected one\n", argv[j]); return 1;}

  if(dst>=0.0) status = 0;
  return status;
}
#endif  
