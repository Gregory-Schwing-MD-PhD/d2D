/* 
 * compile with
 *  gcc-4.2  tamiola.c -lm -O9 -Wall -Wextra -ansi -pedantic -o tamiola-coil.x
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
  int    nseq[MAXLENGTH];
  int  err, ashift;
  unsigned int i, sql, count;
  char n6[25];
  char   *format=NULL, seq[MAXLENGTH], str[MAXLENGTH], *dbpath;
  double mean[21], sx[21], dx[21];
  double ecs[MAXLENGTH], tcs[MAXLENGTH];
  double tcs2[MAXLENGTH], ecs2[MAXLENGTH];
  double cor=0, mx=0, my=0, mxx=0, myy=0, t1, t2;
  double chi;
  FILE *m, *ref, *pred;

  dbpath = getenv ("CAMDBPATH");
  for(ashift=0;ashift<NSHIFT;ashift++) {
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

    switch(ashift)
    {
      case CA: sprintf(str,"%s/tamiola-db/sCAdata.dat",dbpath); break;
      case CB: sprintf(str,"%s/tamiola-db/sCBdata.dat",dbpath); break;
      case CO: sprintf(str,"%s/tamiola-db/sCOdata.dat",dbpath); break;
      case HN: sprintf(str,"%s/tamiola-db/sHNdata.dat",dbpath); break;
      case HA: sprintf(str,"%s/tamiola-db/sHAdata.dat",dbpath); break;
      case NH: sprintf(str,"%s/tamiola-db/sNHdata.dat",dbpath); break;
    }

    m = fopen(str,"r");

    for(i=0;i<21;i++)
      err = fscanf(m, "%*s %lf %lf %lf", &sx[i], &mean[i], &dx[i]);

    fclose(m);

    count=0;
    for(i=0;i<sql;i++)
    {
      if(nseq[i]==G&&(ashift==HA||ashift==CB)) {ecs[i]=999.999; continue;}
      if(nseq[i]==P&&(ashift==HN||ashift==NH)) {ecs[i]=999.999; continue;}

      if(nseq[i]==Z) continue;
      else tcs[i] = mean[nseq[i]]; 

      if(i>0&&nseq[i-1]!=Z) tcs[i] += sx[nseq[i-1]];
      if(i<(strlen(seq)-1)&&nseq[i+1]!=Z) tcs[i] += dx[nseq[i+1]];

      if(ecs[i]<999) {ecs2[count]=ecs[i]; tcs2[count]=tcs[i]; count++;}
    }

    cor=0; mx=0; my=0; mxx=0; myy=0;

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
