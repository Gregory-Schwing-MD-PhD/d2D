/* 
 * compile with
 *  gcc  ppii_pred.c -lm -O9 -Wall -Wextra -ansi -pedantic -o ppii.x
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAXLENGTH 70000

enum {
  A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y, X, Z, NR
};
enum {
  CA, CB, CO, HN, HA, NH, NSHIFT
};

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

int do_test(char *inputfile)
{
  unsigned int  sql, i, count;
  int err, ashift;
  int    nseq[MAXLENGTH];
  char n1[21], n2[21], n3[21], n4[25], n5[25], n6[25];
  char   *format, seq[MAXLENGTH], str[MAXLENGTH], *dbpath;
  double mean[21], mean2[21];
  double ecs[MAXLENGTH], tcs[MAXLENGTH];
  double tcs2[MAXLENGTH], ecs2[MAXLENGTH];
  double s1[21][21], s2[21][21], d1[21][21], d2[21][21];
  double as1, as2, ad1, ad2;
  double chi;
  double cor=0, mx=0, my=0, mxx=0, myy=0, t1, t2;
  FILE *m, *ref, *uno, *due, *tre, *qua, *par, *pred;

  dbpath = getenv ("CAMDBPATH");
  for(ashift=0;ashift<NSHIFT;ashift++) {
    switch(ashift)
    {
      case CA: strcpy(n1,"/ppii-db/CA-uno.dat"); strcpy(n2,"/ppii-db/CA-due.dat"); strcpy(n3,"/ppii-db/CA-tre.dat"); 
               strcpy(n4,"/ppii-db/CA-quattro.dat"); strcpy(n5,"/ppii-db/CA-parm.dat"); break;
      case CB: strcpy(n1,"/ppii-db/CB-uno.dat"); strcpy(n2,"/ppii-db/CB-due.dat"); strcpy(n3,"/ppii-db/CB-tre.dat"); 
               strcpy(n4,"/ppii-db/CB-quattro.dat"); strcpy(n5,"/ppii-db/CB-parm.dat"); break;
      case CO: strcpy(n1,"/ppii-db/CO-uno.dat"); strcpy(n2,"/ppii-db/CO-due.dat"); strcpy(n3,"/ppii-db/CO-tre.dat"); 
               strcpy(n4,"/ppii-db/CO-quattro.dat"); strcpy(n5,"/ppii-db/CO-parm.dat"); break;
      case HN: strcpy(n1,"/ppii-db/HN-uno.dat"); strcpy(n2,"/ppii-db/HN-due.dat"); strcpy(n3,"/ppii-db/HN-tre.dat"); 
               strcpy(n4,"/ppii-db/HN-quattro.dat"); strcpy(n5,"/ppii-db/HN-parm.dat"); break;
      case HA: strcpy(n1,"/ppii-db/HA-uno.dat"); strcpy(n2,"/ppii-db/HA-due.dat"); strcpy(n3,"/ppii-db/HA-tre.dat"); 
               strcpy(n4,"/ppii-db/HA-quattro.dat"); strcpy(n5,"/ppii-db/HA-parm.dat"); break;
      case NH: strcpy(n1,"/ppii-db/NH-uno.dat"); strcpy(n2,"/ppii-db/NH-due.dat"); strcpy(n3,"/ppii-db/NH-tre.dat"); 
               strcpy(n4,"/ppii-db/NH-quattro.dat"); strcpy(n5,"/ppii-db/NH-parm.dat"); break;
    }
    sprintf(str, "%s%s", dbpath, n1);
    uno = fopen(str,"r");
    sprintf(str, "%s%s", dbpath, n2);
    due = fopen(str,"r");
    sprintf(str, "%s%s", dbpath, n3);
    tre = fopen(str,"r");
    sprintf(str, "%s%s", dbpath, n4);
    qua = fopen(str,"r");
    sprintf(str, "%s%s", dbpath, n5);
    par = fopen(str,"r");

    switch(ashift)
    {
      case CA: format =  "%*i %*i %c %*c %*c %lf %*lf %*lf %*lf %*lf %*lf"; strcpy(n6,"pred-ca.dat"); break;
      case CB: format =  "%*i %*i %c %*c %*c %*lf %lf %*lf %*lf %*lf %*lf"; strcpy(n6,"pred-cb.dat"); break;
      case CO: format =  "%*i %*i %c %*c %*c %*lf %*lf %lf %*lf %*lf %*lf"; strcpy(n6,"pred-co.dat"); break;
      case HN: format =  "%*i %*i %c %*c %*c %*lf %*lf %*lf %lf %*lf %*lf"; strcpy(n6,"pred-hn.dat"); break;
      case HA: format =  "%*i %*i %c %*c %*c %*lf %*lf %*lf %*lf %lf %*lf"; strcpy(n6,"pred-ha.dat"); break;
      case NH: format =  "%*i %*i %c %*c %*c %*lf %*lf %*lf %*lf %*lf %lf"; strcpy(n6,"pred-nh.dat"); break;
    }

    ref = fopen(inputfile,"r");

    sql=0;
    while(1)
    {
      err = fscanf(ref, format, &seq[sql], &ecs[sql]);
      if(err==EOF) break;
      sql++;
    }
    fclose(ref);
    fprintf(stderr,"Reference database length %i. Last residue is %c\n", sql, seq[sql-1]);

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
        case 'O': nseq[i] = P; break;
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
        default: fprintf(stderr, "ERROR: AminoAcids %c not recognized! line %d\n",seq[i], i); return 1; break;
      }
    }
    err = fscanf(uno, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
    err = fscanf(due, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
    err = fscanf(tre, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
    err = fscanf(qua, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
 
    format = "%*s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf";
    switch(ashift)
    {
      case CA: sprintf(str,"%s/ppii-db/cs-ca-medi.dat",dbpath); break;
      case CB: sprintf(str,"%s/ppii-db/cs-cb-medi.dat",dbpath); break;
      case CO: sprintf(str,"%s/ppii-db/cs-co-medi.dat",dbpath); break;
      case HN: sprintf(str,"%s/ppii-db/cs-hn-medi.dat",dbpath); break;
      case HA: sprintf(str,"%s/ppii-db/cs-ha-medi.dat",dbpath); break;
      case NH: sprintf(str,"%s/ppii-db/cs-n-medi.dat",dbpath); break;
    }

    m = fopen(str,"r");

    for(i=0;i<21;i++)
    {
      err = fscanf(m, "%*s %lf %lf", &mean[i], &mean2[i]);

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

    count=0;
    for(i=0;i<sql;i++)
    {
      if(nseq[i]==G&&(ashift==HA||ashift==CB)) {ecs[i]=999.999; continue;}
      if(nseq[i]==P&&(ashift==HN||ashift==NH)) {ecs[i]=999.999; continue;}

      if(nseq[i]==Z) continue;
      else tcs[i] = mean[nseq[i]]; 

      if(i>1&&nseq[i-2]!=Z&&nseq[i-1]!=Z) tcs[i] += as2*s2[nseq[i-2]][nseq[i]];
      if(i>0&&nseq[i-1]!=Z) tcs[i] += as1*s1[nseq[i-1]][nseq[i]];
      if(i<(strlen(seq)-1)&&nseq[i+1]!=Z) tcs[i] += ad1*d1[nseq[i+1]][nseq[i]];
      if(i<(strlen(seq)-2)&&nseq[i+2]!=Z&&nseq[i+1]!=Z) tcs[i] += ad2*d2[nseq[i+2]][nseq[i]];
      if(ecs[i]<999&&fabs(ecs[i]-mean[nseq[i]])<3.*mean2[nseq[i]]) {ecs2[count]=ecs[i]; tcs2[count]=tcs[i]; count++;}
    }

    mx=0; my=0; cor=0; mxx=0; myy=0;
    for(i=0;i<count;i++) {
      mx += tcs2[i];	
      my += ecs2[i];
    }
    mx /= count;
    my /= count;

    for(i=0;i<count;i++) { 
      t1 = (tcs2[i]-mx); 
      t2 = (ecs2[i]-my); 
      mxx += t1*t1;
      myy += t2*t2;
      cor += t1*t2;
    }

    cor =  cor/sqrt(mxx*myy);

    chi = 0.;
    count = 0;
    for(i=1;i<strlen(seq)-1;i++) 
    {
      if(ecs[i]<900) {
        chi += ((ecs[i]-tcs[i])*(ecs[i]-tcs[i]));
        count++;
      }
    }
    chi = sqrt(chi/count);
    fprintf(stderr,"CS %s Err %f Cor %f\n", shtostring(ashift), chi, cor);

    pred = fopen(n6, "w"); 
    for(i=0;i<sql;i++) 
    {
      if(ecs[i]<999) 
        fprintf(pred, "%i \t %c \t %10.3f %10.3f \n", i+1, seq[i], tcs[i], ecs[i]);
    }
    fclose(pred);

  }
  return 0;
}

int main(int argc, char *argv[])
{
  int status;

  status = 0;
  fprintf(stderr,"Arguments: %d %s %s\n", argc, argv[0], argv[1]);

  status = do_test(argv[1]);

  return status ;
}  
