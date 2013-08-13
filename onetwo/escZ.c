/*********  Autonomous Routines for Real Space Equilibrium ********/
#ifdef Doc_ESI
\begin{verbatim}
/*------------------------ FORTRAN USE ----------------------*/
      program esitest
      implicit none
      integer i,j,k
      real*8 da,dgq

      integer na1,np1,n
      parameter (na1=21,np1=65,n=na1*np1)
      integer sw(n)
      real*8 q(na1)
      real*8 r(n),ra(n),rq(n),z(n),za(n),zq(n),B(n),Ba(n),Bq(n),gh(n)
     &     ,gha(n),ghq(n),ne(n),dne(n)
      real*8 a(n),gq(n),F(n),Fa(n),gFa(n),gFaa(n),gYa(n),gYaa(n),T(n)
     &     ,Ta(n),P(n),Pa(n)

      integer esiread,esigetprofiles,esiget2dfunctions
     &     ,esigetrzb,gcmotion
      external esiread,esilink2c,esigetprofiles,esiget2dfunctions
     &     ,esigetrzb,gcmotion

c1.     Making link with C-ESI. Provides addresses of arrays to be used by 
ccc     esiget2dfunctions() to direct its output. Should be called only 
ccc     once (if there is no intentional redirection of the output).
      call esilink2c4ne(
     &ne,dne)		! density and its radial derivative

      call esilink2c(
     &	   F,Fa		! F(a)=rBtor function and its radial derivative
     &     ,gFa,gFaa	! first and second derivatives of the toroidal flux/2pi
     &     ,gYa,gYaa	! first and second derivatives of the poloidal flux/2pi
     &     ,T,Ta	! Shafranovs's FF'(Psi) and its first derivative
     &     ,P,Pa	! Shafranovs's p'(Psi) and its first derivative
     &     ,r,ra,rq	! cylindrical r(a,q), and its first derivatives
     &     ,z,za,zq	! z(a,q), and its first derivatives
     &     ,B,Ba,Bq	! module B(a,q), and its first derivative
     &     ,gh,gha,ghq	! correction for PEST coordinates and its derivatives
     &     ,sw) 	! integer label of the particle. Make it sw(i)=0
c All radial derivatives are taken with respect to a normalized 0<= a <=1

c2.    Reading data file (s) File names: <A-Z,a-z>,<0-9>,<_>,<.>
      k=esiread('esiA.00 ') 	! one read per equiulibrium

      if(k.ne.0) then
         write(*,*)'bad file for reading'
         stop
      endif

ccc   Example of getting plasma profiles as a function of minor radius

      write(*,'(a2,7a10)')'i','a','F','gFa','gYa','q','T','P'
      da=1./(na1-1)
      do i =1,na1
         a(i)=da*(i-1)
         k=esigetprofiles(F(i),Fa(i), gFa(i),gFaa(i), gYa(i),gYaa(i)
     &        ,T(i),Ta(i), P(i),Pa(i), a(i))
         if(i.eq.1) then
            q(i)=-gYaa(i)/gFaa(i)
         else
            q(i)=-gYa(i)/gFa(i)
         endif
         write(*,'(i2,1p7e10.3)')i,a(i),F(i),gFa(i),gYa(i),q(i),T(i),P(i
     &        )
      enddo
ccc

ccc   Example of setting coordinates of particles

      da=0.8/(na1-1)
      dgq=8.*datan(1.d0)/(np1-1)
      k=0
      do i=1,na1
         do j=1,np1
            k=k+1
            a(k)=da*(i-1)+0.1
            gq(k)=dgq*(j-1)+0.2
         enddo
      enddo
ccc

c3.   Example of getting 2-D functions. Output goes to the arrays given by
ccc   esilink2c(). Designed for multiple calls inside the time step loops.

      k=esiget2dfunctions(a,gq,n)

      if(k.ne.0) then
         write(*,*)'wrong a or gq for particle #',k
         call esifree()
         stop
      endif

      write(*,'(1p6a10)')'r','z','ra','za','rq','zq'
      do i=4,6
         j=np1*i
         write(*,'(a,i2)')"i=",i
         write(*,'(1p6e10.3)')(r(k),z(k),ra(k),za(k),rq(k),zq(k)
     &        ,k=j-np1+1,j)
         j=j+np1
      enddo
c end of 3.

c4.   Example of getting r,z,B in a number of points

      k=esigetrzb(r,z,B,a,gq,n)
      if(k.ne.0) then
         write(*,*)'wrong a or gq for particle #',k
         call esifree()
         stop
      endif
c end of 4.

c5.   Example of getting time derivatives of the guiding center 
c     coordinates gr_parallel,a,theta,phi of n particles with magnetic 
c     momentum mu

      k=gcmotion(dgr,da,dgq,dgf,gr,a,gq,gm,n)
c end of 5

c6    Example of setting output for density value at the particle position
c     next routine sends adresses of arrays ne,dne for density and it 
c     derivative to ESI
c     Should be called just once


c     Specifying tempterature at the point (nT is known only by ESI)

      Te=3.
      Ti=3.
      call settemperatures(Te,Ti)
c     After this density will be calculated during calls of 
c     esiget2dfunctions()

      call esifree()
c end of 2.   Freeing ESI    
      end
\end{verbatim}

#endif
/*===================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifndef single
typedef double REAL;
#else
typedef float REAL;
#endif

double *ESsr,*ESsz,*ESaB;

static double cr2,cr3,cr4,cr6,cr12,cgp,c2gp,cgp_4,cr2gp,cgm0,crgm0;
static double c2rqgp;

static int Na,Na1=0,Np,Np1=0,MemFl=0;
static REAL *sa,*gq,*aF,*daF,*dgF,*d2gF,*dgY,*d2gY,*aT,*daT,*aP,*daP;
static REAL *gF,*gY;
static REAL *sne,*dsne;
static REAL *sr,*sra,*srq,*sraq;
static REAL *sz,*sza,*szq,*szaq;
static REAL *aB,*aBa,*aBq,*aBaq;
static REAL *gH,*gHa,*gHq,*gHaq,*gHqq,*gHaqq;

static int i0=0,i1=1,j00,j10,j01,j11; 
static REAL A,ha,rha,hq,rhq,cgq0,crgq0;/* period and inverse period in gq */
static REAL xF,XF,xD,XD,dX,dxD,dXD,yF,YF,yD,YD,dY,dyD,dYD;
static REAL y,Y,yy,YY;

static int nF=128;
static double F0[128],F1[128],F2[128],F3[128],tF1=5.,hF,rhF;

int InitLocalConstants()
{
  cr2	=0.5;
  cr3	=1./3.;
  cr4	=0.25;
  cr6	=1./6.;
  cr12	=1./12.;
  cgp_4	=atan((double)1.);
  cgp	=4.*cgp_4;
  c2gp	=8.*cgp_4;
  cr2gp	=1./c2gp;
  cgm0	=0.4*cgp;
  crgm0	=1./cgm0;

  c2rqgp=2./sqrt(cgp);
  return(0);
}

int esifree_()
{
  if((MemFl&0x20)) free(gF);
  if((MemFl&0x10)) free(gH);
  if((MemFl&0x08)) free(aB);
  if((MemFl&0x04)) free(sz);
  if((MemFl&0x02)) free(sr);
  if((MemFl&0x01)) free(gq);
  MemFl=0;
  return(0);
}

int SetFunctionErf()
{
  int i;
  double x,t,dt,s,r6h,r2h,F,dF;
  double S,dY0,dY1,d2Y0,d2Y1,f;
  
  hF	=tF1/nF;
  rhF	=1./hF;
  r2h	=0.5*hF;
  r6h	=hF/6.;

  S	=0.;
  dY1	=0.;
  d2Y1	=0.;

  f	=2.*cr3*c2rqgp;
  F	=0.;
  dF	=0.;
  for(i=0; i < nF; i++){
    F0[i]	=F;
    F1[i]	=dF;
    t	=hF*(i+1);
    dY0	=dY1;
    d2Y0=d2Y1;
    x	=t*t;
    s	=exp(-x);
    dY1	=x*x*x*s;
    d2Y1=t*x*x*s*(6.-2.*x);
    S	+=(dY0+dY1+r6h*(d2Y0-d2Y1))*r2h;
    F	=f*(s*(1.+0.4*x)+0.8*S/(t*x));
    dF	=F-t*f*(1.2*t*s+2.4*S/(x*x));
    F	*=t;
    F2[i]=3.*(F-F0[i])-hF*(2.*F1[i]+dF);
    F3[i]=-2.*(F-F0[i])+hF*(F1[i]+dF);
  }
  return(0);
}

int ESISetMem(double *gq0, int na1,int np1)
{
  int i;

  Np1	=np1;
  Na1	=na1;
  Na	=Na1-1;
  Np	=Np1-1;

  InitLocalConstants();
  SetFunctionErf();

  gq	=gq0;
  sa	=gq	+Np1;
  aF	=sa	+Na1;
  daF	=aF	+Na1;
  dgF	=daF	+Na1;
  d2gF	=dgF	+Na1;
  dgY	=d2gF	+Na1;
  d2gY	=dgY	+Na1;
  aT	=d2gY	+Na1;
  daT	=aT	+Na1;
  aP	=daT	+Na1;
  daP	=aP	+Na1;

  i	=Na1*Np1;
  sr	=daP	+Na1;
  sra	=sr	+i;
  srq	=sra	+i;
  sraq	=srq	+i;

  sz	=sraq	+i;
  sza	=sz	+i;
  szq	=sza	+i;
  szaq	=szq	+i;

  aB	=szaq	+i;
  aBa	=aB	+i;
  aBq	=aBa	+i;
  aBaq	=aBq	+i;

  gH	=aBaq	+i;
  gHa	=gH	+i;
  gHq	=gHa	+i;
  gHaq	=gHq	+i;
  gHqq	=gHaq	+i;
  gHaqq	=gHqq	+i;

  ESsr	=sr;
  ESsz	=sz;
  ESaB	=aB;
  if((MemFl&0x20)){
    free(gF);
  }
  gF	=(REAL*)malloc(2*Na1*sizeof(REAL));
  gY	=gF	+Na1;
  if(gF != NULL){
    MemFl	|=0x20;
  }
  return(0);
}

int ESIMemAlloc(int na1,int np1)
{
  int i,j;

  Np1	=np1;
  Na1	=na1;
  Na	=Na1-1;
  Np	=Np1-1;

  InitLocalConstants();
  SetFunctionErf();
  gq	=(REAL*)malloc((Np1+13*Na1)*sizeof(REAL));
  if(gq != NULL){
    MemFl	+=0x01;
  }
  i	=Na1*Np1;
  j	=4*i;
  sr	=(REAL *)malloc(j*sizeof(REAL));
  if(sr != NULL){
    MemFl	|=0x02;
  }
  sz	=(REAL *)malloc(j*sizeof(REAL));
  if(sz != NULL){
    MemFl	|=0x04;
  }
  aB	=(REAL *)malloc(j*sizeof(REAL));
  if(aB != NULL){
    MemFl	|=0x08;
  }
  gH	=(REAL *)malloc((j+2*i)*sizeof(REAL));
  if(gH != NULL){
    MemFl	|=0x10;
  }
  if(MemFl < 0x1F){
    printf("No memory for BiCubic ESI MemFl=0x%x%c\n",7,MemFl);
    esifree_();
    return(1);
  }

  sa	=gq	+Np1;
  aF	=sa	+Na1;
  daF	=aF	+Na1;
  dgF	=daF	+Na1;
  d2gF=dgF	+Na1;
  dgY	=d2gF	+Na1;
  d2gY=dgY	+Na1;
  aT	=d2gY	+Na1;
  daT	=aT	+Na1;
  aP	=daT	+Na1;
  daP	=aP	+Na1;
  sne	=daP	+Na1;
  dsne	=sne	+Na1;
  
  sra	=sr	+i;
  srq	=sra	+i;
  sraq=srq	+i;
  sza	=sz	+i;
  szq	=sza	+i;
  szaq=szq	+i;
  aBa	=aB	+i;
  aBq	=aBa	+i;
  aBaq=aBq	+i;
  gHa	=gH	+i;
  gHq	=gHa	+i;
  gHaq=gHq	+i;
  gHqq=gHaq	+i;
  gHaqq=gHqq	+i;

  return(0);
}

int ESICopy(double *esr,double *esz,double *eaB,double *egH,int na1,int np1)
{
  int i,n;
  double *esra,*esrt,*esrat,*esza,*eszt,*eszat,*eaBa,*eaBt,*eaBat
    ,*egHa,*egHt,*egHat;

  n	=np1*na1;
  esra	=esr	+n;
  esrt	=esra	+n;
  esrat	=esrt	+n;
  esza	=esz	+n;
  eszt	=esza	+n;
  eszat	=eszt	+n;
  eaBa	=eaB	+n;
  eaBt	=eaBa	+n;
  eaBat	=eaBt	+n;
  egHa	=egH	+n;
  egHt	=egHa	+n;
  egHat	=egHt	+n;
  for(i=0; i < n; i++){
    esr[i]	=sr[i];
    esra[i]	=sra[i];
    esrt[i]	=-srq[i];
    esrat[i]	=-sraq[i];
    esz[i]	=sz[i];
    esza[i]	=sza[i];
    eszt[i]	=-szq[i];
    eszat[i]	=-szaq[i];
    eaB[i]	=aB[i];
    eaBa[i]	=aBa[i];
    eaBt[i]	=-aBq[i];
    eaBat[i]	=-aBaq[i];
    egH[i]	=-gH[i];
    egHa[i]	=-gHa[i];
    egHt[i]	=gHq[i];
    egHat[i]	=gHaq[i];

    egH[i]	=gHq[i];
    egHa[i]	=gHaq[i];
    egHt[i]	=-gHqq[i];
    egHat[i]	=-gHaqq[i];
  }
  return(0);
}

/* Binary Double */
int ESIReadBD(char* FNm) 
{
  int i,j,k;
  FILE *lf;
  REAL s,ss;

  lf	=fopen(FNm,"r");
  if(lf == NULL){
    printf("%s !!! cannot be open%c\n",FNm,7);
    return(1);
  }
 
  fread(&j,sizeof(int),1,lf);
  fread(&i,sizeof(int),1,lf);
  if(MemFl == 0 || Np1 != j || Na1 != i){
    if(MemFl){
      esifree_();
    }
    if(ESIMemAlloc(i,j)){
      return(1);
    }
  } 
  fread(gq,sizeof(REAL),Np1,lf);
  fread(sa,sizeof(REAL),Na1,lf);
  fread(aF,sizeof(REAL),Na1,lf);
  fread(daF,sizeof(REAL),Na1,lf);
  fread(dgF,sizeof(REAL),Na1,lf);
  fread(d2gF,sizeof(REAL),Na1,lf);
  fread(dgY,sizeof(REAL),Na1,lf);
  fread(d2gY,sizeof(REAL),Na1,lf);
  fread(aT,sizeof(REAL),Na1,lf);
  fread(daT,sizeof(REAL),Na1,lf);
  fread(aP,sizeof(REAL),Na1,lf);
  fread(daP,sizeof(REAL),Na1,lf);
  j	=0;
  for(i=0; i < Na1; i++){
    fread(sr+j,sizeof(REAL),Np1,lf);
    fread(sra+j,sizeof(REAL),Np1,lf);
    fread(srq+j,sizeof(REAL),Np1,lf);
    fread(sraq+j,sizeof(REAL),Np1,lf);
    fread(sz+j,sizeof(REAL),Np1,lf);
    fread(sza+j,sizeof(REAL),Np1,lf);
    fread(szq+j,sizeof(REAL),Np1,lf);
    fread(szaq+j,sizeof(REAL),Np1,lf);
    fread(aB+j,sizeof(REAL),Np1,lf);
    fread(aBa+j,sizeof(REAL),Np1,lf);
    fread(aBq+j,sizeof(REAL),Np1,lf);
    fread(aBaq+j,sizeof(REAL),Np1,lf);
    fread(gHq+j,sizeof(REAL),Np1,lf);
    fread(gHaq+j,sizeof(REAL),Np1,lf);
    fread(gHqq+j,sizeof(REAL),Np1,lf);
    if(fread(gHaqq+j,sizeof(REAL),Np1,lf) != Np1){
      fclose(lf);
      printf("%s Bad binary file\n",FNm);
      return(1);
    }
    j	+=Np1;
  }
  fclose(lf);
  return(0);
}

int ESIReadAD(char* FNm) 
{
  int i,j,k;
  double s;
  FILE *lf;
  char ln[256],*lc;

  lf	=fopen(FNm,"r");
  if(lf == NULL){
    printf("File=<%s> !!! cannot be open%c\n",FNm,7);
    return(1);
  }
  lc	=ln;
  while(feof(lf) == 0 && (*lc=fgetc(lf)) != '\n'){
    if(lc-ln < 25){
      lc++;
    }
  }
  *lc	='\0';
  i	=lc-ln;
  lc	="!!! Do not edit this file";
  if(i < strlen(lc) || strncmp(ln,lc,strlen(lc))){
#ifdef H
    printf("%s first line is not !!! Do not edit this file%c\n",FNm,7);
#endif
    fclose(lf);
    return(2);
  }
  if(fscanf(lf,"%d x%d",&j,&i) != 2){
    printf("%s - wrong data for Np1, Na1%c\n",FNm,7);
    return(1);
  }
  lc	=ln;
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  if(MemFl == 0 || Np1 != j || Na1 != i){
    if(MemFl){
      esifree_();
    }
    if(ESIMemAlloc(i,j)){
      return(1);
    }
  } 
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  for(j=0; j < Np1; j++){
    if(fscanf(lf,"%lg",gq+j) != 1){
      printf("%s - wrong data for gq[%d]%c\n",FNm,j,7);
      return(1);
    }
  }
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }

  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  for(i=0; i < Na1; i++){
    if(fscanf(lf,"%lg%lg%lg",sa+i,aF+i,daF+i) != 3){
      printf("%s - wrong data for a,F,F'[%d]%c\n",FNm,i,7);
      return(1);
    }
  }
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }

  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  for(i=0; i < Na1; i++){
    if(fscanf(lf,"%lg%lg%lg%lg",dgF+i,d2gF+i,dgY+i,d2gY+i) != 4){
      printf("%s - wrong data for gF',gF'',gY',gY'',[%d]%c\n",FNm,i,7);
      return(1);
    }
  }
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }

  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }

  for(i=0; i < Na1; i++){
    if(fscanf(lf,"%lg%lg%lg%lg",aT+i,daT+i,aP+i,daP+i) != 4){
      printf("%s - wrong data for T,T',P,P'[%d]%c\n",FNm,i,7);
      return(1);
    }
  }
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }

  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  k	=0;
  for(i=0; i < Na1; i++){
    for(j=0; j < Np1; j++){
      if(fscanf(lf,"%lg%lg%lg%lg",sr+k,sra+k,srq+k,sraq+k) != 4){
	printf("%s - wrong data [j=%d,i=%d] for r,r'_a,r'_gq,r''_{a,gq}%c\n"
	       ,FNm,j,i,7);
	return(	1);
      }
      k++;
    }
  }
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }

  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  k	=0;
  for(i=0; i < Na1; i++){
    for(j=0; j < Np1; j++){
      if(fscanf(lf,"%lg%lg%lg%lg",sz+k,sza+k,szq+k,szaq+k) != 4){
	printf("%s - wrong data [j=%d,i=%d] for z,z'_a,z'_gq,z''_{a,gq}%c\n"
	       ,FNm,j,i,7);
	return(	1);
      }
      k++;
    }
  }
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }

  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  k	=0;
  for(i=0; i < Na1; i++){
    for(j=0; j < Np1; j++){
      if(fscanf(lf,"%lg%lg%lg%lg",aB+k,aBa+k,aBq+k,aBaq+k) != 4){
	printf("%s - wrong data [j=%d,i=%d] for B,B'_a,B'_gq,B''_{a,gq}%c\n"
	       ,FNm,j,i,7);
	return(	1);
      }
      k++;
    }
  }
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }

  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  k	=0;
  for(i=0; i < Na1; i++){
    for(j=0; j < Np1; j++){
      if(fscanf(lf,"%lg%lg%lg%lg",gHq+k,gHaq+k,gHqq+k,gHaqq+k) != 4){
	printf("%s -wrong data [j=%d,i=%d] for gh,gh'_a,gh'_gq,gh''_{a,gq}%c\n"
	       ,FNm,j,i,7);
	return(	1);
      }
      k++;
    }
  }
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }

  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  for(i=0; i < Na1; i++){
    if(fscanf(lf,"%lg%lg%lg",&s,sne+i,dsne+i) != 3){
      if(i == 0){
	break;
      }
      printf("%s -wrong data [i=%d] for ne,dne%c\n",FNm,i,7);
      return(1);
    }
  }
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  fclose(lf);
  return(0);
}

#ifdef single
/* Binary Float */
int ESIReadBF(char* FNm) 
{
  int i,j,k;
  FILE *lf;
  REAL s,ss,*ld;
  double *d;

  lf	=fopen(FNm,"r");
  if(lf == NULL){
    printf("%s !!! cannot be open%c\n",FNm,7);
    return(1);
  }
  fread(&j,sizeof(int),1,lf);
  fread(&i,sizeof(int),1,lf);
  if(MemFl == 0 || Np1 != j || Na1 != i){
    if(MemFl){
      esifree_();
    }
    if(ESIMemAlloc(i,j)){
      return(1);
    }
  } 
  d	=(double*)malloc(Np1*sizeof(double));
  fread(d,sizeof(double),Np1,lf);
  ld	=gq;
  for(k=0; k < Np1; k++){
    *ld++	=d[k];
  }
  fread(d,sizeof(double),Na1,lf);
  ld	=sa;
  for(k=0; k < Na1; k++){
    *ld++	=d[k];
  }
  fread(d,sizeof(double),Na1,lf);
  ld	=aF;
  for(k=0; k < Na1; k++){
    *ld++	=d[k];
  }
  fread(d,sizeof(double),Na1,lf);
  ld	=daF;
  for(k=0; k < Na1; k++){
    *ld++	=d[k];
  }
  fread(d,sizeof(double),Na1,lf);
  ld	=dgF;
  for(k=0; k < Na1; k++){
    *ld++	=d[k];
  }
  fread(d,sizeof(double),Na1,lf);
  ld	=d2gF;
  for(k=0; k < Na1; k++){
    *ld++	=d[k];
  }
  fread(d,sizeof(double),Na1,lf);
  ld	=dgY;
  for(k=0; k < Na1; k++){
    *ld++	=d[k];
  }
  fread(d,sizeof(double),Na1,lf);
  ld	=d2gY;
  for(k=0; k < Na1; k++){
    *ld++	=d[k];
  }
  fread(d,sizeof(double),Na1,lf);
  ld	=aT;
  for(k=0; k < Na1; k++){
    *ld++	=d[k];
  }
  fread(d,sizeof(double),Na1,lf);
  ld	=daT;
  for(k=0; k < Na1; k++){
    *ld++	=d[k];
  }
  fread(d,sizeof(double),Na1,lf);
  ld	=aP;
  for(k=0; k < Na1; k++){
    *ld++	=d[k];
  }
  fread(d,sizeof(double),Na1,lf);
  ld	=daP;
  for(k=0; k < Na1; k++){
    *ld++	=d[k];
  }
  j	=0;
  for(i=0; i < Na1; i++){
    fread(d,sizeof(double),Np1,lf);
    ld	=sr+j;
    for(k=0; k < Na1; k++){
      *ld++	=d[k];
    }
    fread(d,sizeof(double),Np1,lf);
    ld	=sra+j;
    for(k=0; k < Na1; k++){
      *ld++	=d[k];
    }
    fread(d,sizeof(double),Np1,lf);
    ld	=srq+j;
    for(k=0; k < Na1; k++){
      *ld++	=d[k];
    }
    fread(d,sizeof(double),Np1,lf);
    ld	=sraq+j;
    for(k=0; k < Na1; k++){
      *ld++	=d[k];
    }
    fread(d,sizeof(double),Np1,lf);
    ld	=sz+j;
    for(k=0; k < Na1; k++){
      *ld++	=d[k];
    }
    fread(d,sizeof(double),Np1,lf);
    ld	=sza+j;
    for(k=0; k < Na1; k++){
      *ld++	=d[k];
    }
    fread(d,sizeof(double),Np1,lf);
    ld	=szq+j;
    for(k=0; k < Na1; k++){
      *ld++	=d[k];
    }
    fread(d,sizeof(double),Np1,lf);
    ld	=szaq+j;
    for(k=0; k < Na1; k++){
      *ld++	=d[k];
    }
    fread(d,sizeof(double),Np1,lf);
    ld	=aB+j;
    for(k=0; k < Na1; k++){
      *ld++	=d[k];
    }
    fread(d,sizeof(double),Np1,lf);
    ld	=aBa+j;
    for(k=0; k < Na1; k++){
      *ld++	=d[k];
    }
    fread(d,sizeof(double),Np1,lf);
    ld	=aBq+j;
    for(k=0; k < Na1; k++){
      *ld++	=d[k];
    }
    fread(d,sizeof(double),Np1,lf);
    ld	=aBaq+j;
    for(k=0; k < Na1; k++){
      *ld++	=d[k];
    }
    fread(d,sizeof(double),Np1,lf);
    ld	=gHq+j;
    for(k=0; k < Na1; k++){
      *ld++	=d[k];
    }
    fread(d,sizeof(double),Np1,lf);
    ld	=gHaq+j;
    for(k=0; k < Na1; k++){
      *ld++	=d[k];
    }
    fread(d,sizeof(double),Np1,lf);
    ld	=gHqq+j;
    for(k=0; k < Na1; k++){
      *ld++	=d[k];
    }
    fread(d,sizeof(double),Np1,lf);
    ld	=gHaqq+j;
    for(k=0; k < Na1; k++){
      *ld++	=d[k];
    }
    j	+=Np1;
  }
  free(d);
  fclose(lf);
  return(0);
}

int ESIReadAF(char* FNm) 
{
  int i,j,k;
  FILE *lf;
  char ln[256],*lc;
  double d[4];

  lf	=fopen(FNm,"r");
  if(lf == NULL){
    printf("%s !!! cannot be open%c\n",FNm,7);
    return(1);
  }
  lc	=ln;
  while(feof(lf) == 0 && (*lc=fgetc(lf)) != '\n'){
    if(lc-ln < 25){
      lc++;
    }
  }
  *lc	='\0';
  i	=lc-ln;
  lc	="!!! Do not edit this file";
  if(i < strlen(lc) || strncmp(ln,lc,strlen(lc))){
    printf("%s first line is not !!! Do not edit this file%c\n",FNm,7);
    fclose(lf);
    return(2);
  }
  if(fscanf(lf,"%d x%d",&j,&i) != 2){
    printf("%s - wrong data on Np, Na%c\n",FNm,7);
    return(1);
  }
  *lc	=ln;
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  if(MemFl == 0 || Np1 != j || Na1 != i){
    if(MemFl){
      esifree_();
    }
    if(ESIMemAlloc(i,j)){
      return(1);
    }
  } 
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  for(j=0; j < Np1; j++){
    if(fscanf(lf,"%lg",d) != 1){
      printf("%s - wrong data for gq[%d]%c\n",FNm,j,7);
      return(1);
    }
    gq[j]	=d[0];
  }
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }

  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  for(i=0; i < Na1; i++){
    if(fscanf(lf,"%lg%lg%lg",d,d+1,d+2) != 3){
      printf("%s - wrong data for a,F,F'[%d]%c\n",FNm,i,7);
      return(1);
    }
    sa[i]	=d[0];
    aF[i]	=d[1];
    daF[i]	=d[2];
  }
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }

  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  for(i=0; i < Na1; i++){
    if(fscanf(lf,"%lg%lg%lg%lg",d,d+1,d+2,d+3) != 4){
      printf("%s - wrong data for gF',gF'',gY',gY''[%d]%c\n",FNm,i,7);
      return(1);
    }
    dgF[i]	=d[0];
    d2gF[i]	=d[1];
    dgY[i]	=d[2];
    d2gY[i]	=d[3];
  }
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }

  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  for(i=0; i < Na1; i++){
    if(fscanf(lf,"%lg%lg%lg%lg",d,d+1,d+2,d+3) != 4){
      printf("%s - wrong data for T,T',P,P'[%d]%c\n",FNm,i,7);
      return(1);
    }
    aT[i]	=d[0];
    daT[i]	=d[1];
    aP[i]	=d[2];
    daP[i]	=d[3];
  }	
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }

  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  k	=0;
  for(i=0; i < Na1; i++){
    for(j=0; j < Np1; j++){
      if(fscanf(lf,"%lg%lg%lg%lg",d,d+1,d+2,d+3) != 4){
	printf("%s - wrong data [j=%d,i=%d] for r,r'_a,r'_gq,r''_{a,gq}%c\n"
	       ,FNm,j,i,7);
	return(	1);
      }
      sr[k]	=d[0];
      sra[k]	=d[1];
      srq[k]	=d[2];
      sraq[k]	=d[3];
      k++;
    }
  }
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }

  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  k	=0;
  for(i=0; i < Na1; i++){
    for(j=0; j < Np1; j++){
      if(fscanf(lf,"%lg%lg%lg%lg",d,d+1,d+2,d+3) != 4){
	printf("%s - wrong data [j=%d,i=%d] for z,z'_a,z'_gq,z''_{a,gq}%c\n"
	       ,FNm,j,i,7);
	return(	1);
      }
      sz[k]	=d[0];
      sza[k]	=d[1];
      szq[k]	=d[2];
      szaq[k]	=d[3];
      k++;
    }
  }
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }

  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  k	=0;
  for(i=0; i < Na1; i++){
    for(j=0; j < Np1; j++){
      if(fscanf(lf,"%lg%lg%lg%lg",d,d+1,d+2,d+3) != 4){
	printf("%s - wrong data [j=%d,i=%d] for B,B'_a,B'_gq,B''_{a,gq}%c\n"
	       ,FNm,j,i,7);
	return(	1);
      }
      aB[k]	=d[0];
      aBa[k]	=d[1];
      aBq[k]	=d[2];
      aBaq[k]	=d[3];
      k++;
    }
  }
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }

  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  k	=0;
  for(i=0; i < Na1; i++){
    for(j=0; j < Np1; j++){
      if(fscanf(lf,"%lg%lg%lg%lg",d,d+1,d+2,d+3) != 4){
	printf("%s -wrong data [j=%d,i=%d] for gh,gh'_a,gh'_gq,gh''_{a,gq}%c\n"
	       ,FNm,j,i,7);
	return(	1);
      }
      gHq[k]	=d[0];
      gHaq[k]	=d[1];
      gHqq[k]	=d[2];
      gHaqq[k]	=d[3];
      k++;
    }
  }
  while(feof(lf) == 0 && (*lc=fgetc(lf)) !='\n'){
    ;
  }
  fclose(lf);
  return(0);
}
#endif

int ESIInit()
{
  int i,i1,j,ji,ji1;
  REAL s,ss;

  cgq0	=gq[Np];
  crgq0	=1./cgq0;
  hq	=gq[1]-gq[0];
  rhq	=1./hq;
  ha	=sa[1]-sa[0];
  rha	=1./ha;
  ji	=0;
  for(i=0; i < Na1; i++){
    s	=0;
    ss	=0;
    gH[ji]	=0.;
    gHa[ji]	=0.;
    for(j=1; j < Np1; j++){
      ji1	=ji;
      ji++;
      gH[ji]	=gH[ji1]+((gHq[ji1]+gHq[ji])+cr6*hq*(gHqq[ji1]-gHqq[ji]))
	*hq*0.5;
      gHa[ji]	=gHa[ji1]+((gHaq[ji1]+gHaq[ji])
			   +cr6*hq*(gHaqq[ji1]-gHaqq[ji]))*hq*0.5;
      s	+=gH[ji];
      ss+=gHa[ji];
    }
    ji++;
    if(i){
      i1	=i-1;
      gF[i]	=gF[i1]
	+(0.5*sa[i1]*(dgF[i1]+dgF[i]+cr6*ha*(d2gF[i1]-d2gF[i]))
	  +0.05*ha*(3.*dgF[i1]+7.*dgF[i]+cr3*ha*(2.*d2gF[i1]-3.*d2gF[i])))*ha;
      gY[i]	=gY[i1]
	+(0.5*(dgY[i1]+dgY[i]+cr6*ha*(d2gY[i1]-d2gY[i]))
	  +0.05*ha*(3.*dgY[i1]+7.*dgY[i]+cr3*ha*(2.*d2gY[i1]-3.*d2gY[i])))*ha;
    }
    else{
      gF[i]	=0.;
      gY[i]	=0.;
    }
#ifdef H
    s	/=Np;
    ss	/=Np;
    ji	-=Np;
    ji++;
    for(j=1; j < Np1; j++){
      gH[ji]	-=s;
      gHa[ji]	-=ss;
      ji++;
    }
#endif
  }
  s	=gY[Na];
  for(i=0; i < Na1; i++){
    gY[i]	-=s;
  }
  return(0);
}

int esiread_(char* FName)
{
  int i;
  char FNm[128],*lc,*ls;
  int *ll0,*ll1,*ll2;

  ls	=FName;
  while(isspace(*ls)){
    ls++;
  }
  lc	=FNm;
  while(*ls != '\0' && lc-FNm < 127 && !isspace(*ls)){
    *lc++	=*ls++;
  }
#ifdef H
  while(!isspace(*ls) && lc-FNm < 127){
    *lc++	=*ls++;
  }
#endif
  *lc	='\0';
#ifndef single
  if((i=ESIReadAD(FNm)) == 2){ 
    i	=ESIReadBD(FNm);
  }
#else
  if((i=ESIReadAF(FNm)) == 2){ 
    i	=ESIReadBF(FNm);
  }
#endif
  if(i){	
    return(i);
  }
  ESIInit();
  return(0);
}

/***************** FORTRAN-C Link Routine ****************************/

static int *sw;
static REAL *r=NULL,*ra=NULL,*rq=NULL,*z=NULL,*za=NULL,*zq=NULL
,*B=NULL,*Ba=NULL,*Bq=NULL,*gh=NULL,*gha=NULL,*ghq=NULL,*ap=NULL,*qp=NULL;
static REAL *F=NULL,*Fa=NULL,*gFa=NULL,*gFaa=NULL,*gYa=NULL,*gYaa=NULL
,*T=NULL,*Ta=NULL,*P=NULL,*Pa=NULL;

static REAL *Ne=NULL,*dNe=NULL;
static REAL aTe=15.,aTi=15.,raT;

void esilink2c_(REAL *XaF,REAL *XaFa	 
		,REAL *XgFa,REAL *XgFaa
		,REAL *XgYa,REAL *XgYaa
		,REAL *XT,REAL *XTa
		,REAL *XP,REAL *XPa
		,REAL *Xr,REAL *Xra,REAL *Xrq
		,REAL *Xz,REAL *Xza,REAL *Xzq
		,REAL *XB,REAL *XBa,REAL *XBq
		,REAL *Xgh,REAL *Xgha,REAL *Xghq
		,int *Xsw)
{
  F	=XaF;
  Fa	=XaFa;
  gFa	=XgFa;
  gFaa	=XgFaa;
  gYa	=XgYa;
  gYaa	=XgYaa;
  T	=XT;
  Ta	=XTa;
  P	=XP;
  Pa	=XPa;
  r	=Xr;
  ra	=Xra;
  rq	=Xrq;
  z	=Xz;
  za	=Xza;
  zq	=Xzq;
  B	=XB;
  Ba	=XBa;
  Bq	=XBq;
  gh	=Xgh;
  gha	=Xgha;
  ghq	=Xghq;
  sw	=Xsw;

  return;
}

void esilink2c4ne_(double *Xne,double *Xdne)
{
  Ne	=Xne;
  dNe	=Xdne;
  return;
}

void settemperatures_(double *Te,double *Ti)
{
  aTe	=*Te;
  aTi	=*Ti;
  raT	=(aTe+aTi) != 0. ? 30./(aTe+aTi) : 1.;
  return;
}

/***************** Main Reconstruction Routines **********************/
int esigetprofiles_(REAL *sF,REAL *sFa,	 
		       /* \baF,\R{d\baF}{da}*/
		       REAL *sgFa,REAL *sgFaa,
		       /* \R{d\bgF}{ada},(\R{d\bgF}{ada})'_a */
		       REAL *sgYa,REAL *sgYaa,
		       /* \R{d\bgY}{ada},(\R{d\bgY}{ada})'_a */
		       REAL *sT,REAL *sTa,
		       /* T=\baF\R{d\baF}{d\bgY},\R{dT}{da} */
		       REAL *sP,REAL *sPa
		       /* P=\R{d\bsp}{d\bgY},\R{dP}{da} */
		       ,REAL *a)
{
  REAL x,X,xx,XX;

  A	=*a;
  if(A < 0.){
    A	=-A;
  }
  i0	=(A-sa[0])*rha;
  if(i0 >= Na){
    i0	=Na-1;
  }
  i1	=i0+1;

  x	=(A-sa[i0])*rha;
  X	=1.-x;

  XX	=X*X;
  xx	=x*x;
  XF	=XX*(3.-2.*X);
  XD	=XX*x*ha;
  xF	=xx*(3.-2.*x);
  xD	=-xx*X*ha;

  dX	=3.*x*X;
  dXD	=X-dX;
  dxD	=x-dX;
  dX	*=2.*rha;

  *sF	=XF*aF[i0]+xF*aF[i1]+XD*daF[i0]+xD*daF[i1];
  *sFa	=dX*(aF[i1]-aF[i0])+dXD*daF[i0]+dxD*daF[i1];

  *sgFa	=XF*dgF[i0]+xF*dgF[i1]+XD*d2gF[i0]+xD*d2gF[i1];
  *sgFaa=dX*(dgF[i1]-dgF[i0])+dXD*d2gF[i0]+dxD*d2gF[i1];
  *sgFaa=A*(*gFaa)+(*gFa);
  *sgYa	=XF*dgY[i0]+xF*dgY[i1]+XD*d2gY[i0]+xD*d2gY[i1];
  *sgYaa=dX*(dgY[i1]-dgY[i0])+dXD*d2gY[i0]+dxD*d2gY[i1];
  *sgYaa=A*(*gYaa)+(*gYa);

  if(A != 0.){
    *sgFa	*=A;
    *sgYa	*=A;
  }
  *sT	=XF*aT[i0]+xF*aT[i1]+XD*daT[i0]+xD*daT[i1];
  *sTa	=dX*(aT[i1]-aT[i0])+dXD*daT[i0]+dxD*daT[i1];

  *sP	=XF*aP[i0]+xF*aP[i1]+XD*daP[i0]+xD*daP[i1];
  *sPa	=dX*(aP[i1]-aP[i0])+dXD*daP[i0]+dxD*daP[i1];

  if(Ne != NULL){
    *Ne	=(XF*sne[i0]+xF*sne[i1]+XD*dsne[i0]+xD*dsne[i1])*raT;
    *dNe=(dX*(sne[i1]-sne[i0])+dXD*dsne[i0]+dxD*dsne[i1])*raT;
  }
  return(0);
}

int esiget2dfunctions0(REAL a,REAL q,int k)
{
  REAL x,X,XX;
  REAL f0,f1,fq0,fq1,r0;
 
  j01	=j00+1;
  j10	=j00+Np1;
  j11	=j10+1;
  
  x	=a*rha;
  X	=1.-x;
  
  XX	=X*X;
  XF	=XX*(3.-2.*X);
  XD	=XX;
  xF	=x*(3.-2.*x)*rha;
  xD	=-x*X;
  
  dX	=3.*x*X;
  dXD	=X-dX;
  dxD	=x-dX;
  dX	*=2.*rha;
  
  r0	=sr[0];
  f0	=xF*sr[j10]+XD*sra[j00]+xD*sra[j10];
  f1	=xF*sr[j11]+XD*sra[j01]+xD*sra[j11];
  fq0	=xF*srq[j10]+XD*sraq[j00]+xD*sraq[j10];
  fq1	=xF*srq[j11]+XD*sraq[j01]+xD*sraq[j11];
  r[k]=(YF+yF)*XF*sr[0]+a*(YF*f0+yF*f1+YD*fq0+yD*fq1);
  rq[k]=dY*(f1-f0)+dYD*fq0+dyD*fq1;
  f0	=dX*(sr [j10]-r0)+dXD*sra [j00]+dxD*sra [j10];
  f1	=dX*(sr [j11]-r0)+dXD*sra [j01]+dxD*sra [j11];
  fq0	=dX*(srq[j10]-srq[j00])+dXD*sraq[j00]+dxD*sraq[j10];
  fq1	=dX*(srq[j11]-srq[j01])+dXD*sraq[j01]+dxD*sraq[j11];
  ra[k]=YF*f0+yF*f1+YD*fq0+yD*fq1;
  
  f0	=xF*sz [j10]+XD*sza [j00]+xD*sza [j10];
  f1	=xF*sz [j11]+XD*sza [j01]+xD*sza [j11];
  fq0	=xF*szq[j10]+XD*szaq[j00]+xD*szaq[j10];
  fq1	=xF*szq[j11]+XD*szaq[j01]+xD*szaq[j11];
  zq[k]=dY*(f1-f0)+dYD*fq0+dyD*fq1;
  z[k]=(YF+yF)*XF*sz[0]+a*(YF*f0+yF*f1+YD*fq0+yD*fq1);
  f0	=dX*(sz [j10]-sz [j00])+dXD*sza [j00]+dxD*sza [j10];
  f1	=dX*(sz [j11]-sz [j01])+dXD*sza [j01]+dxD*sza [j11];
  fq0	=dX*(szq[j10]-szq[j00])+dXD*szaq[j00]+dxD*szaq[j10];
  fq1	=dX*(szq[j11]-szq[j01])+dXD*szaq[j01]+dxD*szaq[j11];
  za[k]=YF*f0+yF*f1+YD*fq0+yD*fq1;
  
  f0	=xF*aB [j10]+XD*aBa [j00]+xD*aBa [j10];
  f1	=xF*aB [j11]+XD*aBa [j01]+xD*aBa [j11];
  fq0	=xF*aBq[j10]+XD*aBaq[j00]+xD*aBaq[j10];
  fq1	=xF*aBq[j11]+XD*aBaq[j01]+xD*aBaq[j11];
  B[k]=(YF+yF)*XF*aB[0]+a*(YF*f0+yF*f1+YD*fq0+yD*fq1);
  Bq[k]=dY*(f1-f0)+dYD*fq0+dyD*fq1;
  f0	=dX*(aB [j10]-aB [j00])+dXD*aBa [j00]+dxD*aBa [j10];
  f1	=dX*(aB [j11]-aB [j01])+dXD*aBa [j01]+dxD*aBa [j11];
  fq0	=dX*(aBq[j10]-aBq[j00])+dXD*aBaq[j00]+dxD*aBaq[j10];
  fq1	=dX*(aBq[j11]-aBq[j01])+dXD*aBaq[j01]+dxD*aBaq[j11];
  Ba[k]=YF*f0+yF*f1+YD*fq0+yD*fq1;
  
  XD	*=a;
  xF	*=a;
  xD	*=a;
  
  F[k]	=XF*aF[i0]+xF*aF[i1]+XD*daF[i0]+xD*daF[i1];
  Fa[k]	=dX*(aF[i1]-aF[i0])+dXD*daF[i0]+dxD*daF[i1];
  gFa[k]	=XF*dgF[i0]+xF*dgF[i1]+XD*d2gF[i0]+xD*d2gF[i1];
  gFaa[k]	=dX*(dgF[i1]-dgF[i0])+dXD*d2gF[i0]+dxD*d2gF[i1];
  gFaa[k]	=a*gFaa[k]+gFa[k];

  gYa[k]	=XF*dgY[i0]+xF*dgY[i1]+XD*d2gY[i0]+xD*d2gY[i1];
  gYaa[k]	=dX*(dgY[i1]-dgY[i0])+dXD*d2gY[i0]+dxD*d2gY[i1];
  gYaa[k]	=a*gYaa[k]+gYa[k];

  T[k]	=XF*aT[i0]+xF*aT[i1]+XD*daT[i0]+xD*daT[i1];
  Ta[k]	=dX*(aT[i1]-aT[i0])+dXD*daT[i0]+dxD*daT[i1];
  P[k]	=XF*aP[i0]+xF*aP[i1]+XD*daP[i0]+xD*daP[i1];
  Pa[k]	=dX*(aP[i1]-aP[i0])+dXD*daP[i0]+dxD*daP[i1];
  if(Ne != NULL){
    Ne[k]=(XF*sne[i0]+xF*sne[i1]+XD*dsne[i0]+xD*dsne[i1])*raT;
    dNe[k]=(dX*(sne[i1]-sne[i0])+dXD*dsne[i0]+dxD*dsne[i1])*raT;
  }
  return(0);
}

int esiget2d_(int m, REAL *F, REAL *Fa, REAL *Fq, REAL *a0,REAL *gq0,int *n)
{
  int j0,j1,k;
  REAL x,X,xx,XX;
  REAL f0,f1,fq0,fq1;
  REAL a,q;
  double *sZ,*sZa,*sZq,*sZaq;
  
  switch(m){
  case 0:
    sZ	=sr;
    sZa	=sra;
    sZq	=srq;
    sZaq=sraq;
    break;
  case 1:
    sZ	=sz;
    sZa	=sza;
    sZq	=szq;
    sZaq=szaq;
    break;
  case 2:
    sZ	=aB;
    sZa	=aBa;
    sZq	=aBq;
    sZaq=aBaq;
    break;
  default:
    break;
  }

  for(k =0; k < *n; k++){
    a	=a0[k];
    q	=gq0[k];
    j1	=q*crgq0;
    if(q < 0.){
      j1--;
    }
    q	-=cgq0*j1;
    j0	=q*rhq;
    if(j0 >= Np){
      j0	-=Np;
      q	-=cgq0;
    }
    if(j0 < 0){
      j0	+=Np;
      q	+=cgq0;
    }
    j1	=j0+1;
    y	=(q-gq[j0])*rhq;
    Y	=1.-y;
   
    YY	=Y*Y;
    yy	=y*y;
    YF	=YY*(3.-2.*Y);
    YD	=YY*y*hq;
    yF	=yy*(3.-2.*y);
    yD	=-yy*Y*hq;
    
    dY	=3.*y*Y;
    dYD	=Y-dY;
    dyD	=y-dY;
    dY	*=2.*rhq;
    if(a < 0.){
      a	=-a;
    }
    i0	=(a-sa[0])*rha;
    if(a != 0.){
      if(i0 >= Na){
	i0	=Na-1;
      }
      i1	=i0+1;
      
      x	=(a-sa[i0])*rha;
      X	=1.-x;

      XX	=X*X;
      xx	=x*x;
      XF	=XX*(3.-2.*X);
      XD	=XX*x*ha;
      xF	=xx*(3.-2.*x);
      xD	=-xx*X*ha;
      
      dX	=3.*x*X;
      dXD	=X-dX;
      dxD	=x-dX;
      dX	*=2.*rha;
      
      j00	=Np1*i0+j0;
      j01	=j00+1;
      j10	=j00+Np1;
      j11	=j10+1;

      f0	=XF*sZ [j00]+xF*sZ [j10]+XD*sZa [j00]+xD*sZa [j10];
      f1	=XF*sZ [j01]+xF*sZ [j11]+XD*sZa [j01]+xD*sZa [j11];
      fq0	=XF*sZq[j00]+xF*sZq[j10]+XD*sZaq[j00]+xD*sZaq[j10];
      fq1	=XF*sZq[j01]+xF*sZq[j11]+XD*sZaq[j01]+xD*sZaq[j11];
      F[k]=YF*f0+yF*f1+YD*fq0+yD*fq1;
      Fq[k]=dY*(f1-f0)+dYD*fq0+dyD*fq1;
      f0	=dX*(sZ [j10]-sZ [j00])+dXD*sZa [j00]+dxD*sZa [j10];
      f1	=dX*(sZ [j11]-sZ [j01])+dXD*sZa [j01]+dxD*sZa [j11];
      fq0	=dX*(sZq[j10]-sZq[j00])+dXD*sZaq[j00]+dxD*sZaq[j10];
      fq1	=dX*(sZq[j11]-sZq[j01])+dXD*sZaq[j01]+dxD*sZaq[j11];
      Fa[k]=YF*f0+yF*f1+YD*fq0+yD*fq1;
    }
    else{	/* a == 0. */
      F[k]	=sZ[0];
      Fa[k]	=YF*sZa[j0]+yF*sZa[j1]+YD*sZaq[j0]+yD*sZaq[j1];
      Fq[k]	=dY*(sZa[j1]-sZa[j0])+dYD*sZaq[j0]+dyD*sZaq[j1];
    }
  }
  return(0);
}

int esiget2dfunctions_(REAL *a0,REAL *gq0,int *n)
{
  /* Exception: at a=0 the routine returns:
     dr_q/a, dz_q/a, dB_q/a, dgF/a, dgY/a
     rather than
     dr_q=0, dz_q=0, dB_q=0, dgF=0, dgY=0 */
  int j0,j1,k;
  REAL x,X,xx,XX;
  REAL f0,f1,fq0,fq1;
  REAL a,q;

  for(k =0; k < *n; k++){
    a	=a0[k];
    q	=gq0[k];
    j1	=q*crgq0;
    if(q < 0.){
      j1--;
    }
    q	-=cgq0*j1;
    j0	=q*rhq;
    if(j0 >= Np){
      j0	-=Np;
      q	-=cgq0;
    }
    if(j0 < 0){
      j0	+=Np;
      q	+=cgq0;
    }
    j1	=j0+1;
    y	=(q-gq[j0])*rhq;
    Y	=1.-y;
   
    YY	=Y*Y;
    yy	=y*y;
    YF	=YY*(3.-2.*Y);
    YD	=YY*y*hq;
    yF	=yy*(3.-2.*y);
    yD	=-yy*Y*hq;
    
    dY	=3.*y*Y;
    dYD	=Y-dY;
    dyD	=y-dY;
    dY	*=2.*rhq;
    if(a < 0.){
      a	=-a;
    }
    i0	=(a-sa[0])*rha;
    if(a != 0.){
      if(i0 >= Na){
	i0	=Na-1;
      }
      i1	=i0+1;
      
      x	=(a-sa[i0])*rha;
      X	=1.-x;

      XX	=X*X;
      xx	=x*x;
      XF	=XX*(3.-2.*X);
      XD	=XX*x*ha;
      xF	=xx*(3.-2.*x);
      xD	=-xx*X*ha;
      
      dX	=3.*x*X;
      dXD	=X-dX;
      dxD	=x-dX;
      dX	*=2.*rha;
      
      j00	=Np1*i0+j0;
      j01	=j00+1;
      j10	=j00+Np1;
      j11	=j10+1;
      
      /* r,r'_a,r'_q */
      f0	=XF*sr [j00]+xF*sr [j10]+XD*sra [j00]+xD*sra [j10];
      f1	=XF*sr [j01]+xF*sr [j11]+XD*sra [j01]+xD*sra [j11];
      fq0	=XF*srq[j00]+xF*srq[j10]+XD*sraq[j00]+xD*sraq[j10];
      fq1	=XF*srq[j01]+xF*srq[j11]+XD*sraq[j01]+xD*sraq[j11];
      r[k]	=YF*f0+yF*f1+YD*fq0+yD*fq1;
      rq[k]	=dY*(f1-f0)+dYD*fq0+dyD*fq1;
      f0	=dX*(sr [j10]-sr [j00])+dXD*sra [j00]+dxD*sra [j10];
      f1	=dX*(sr [j11]-sr [j01])+dXD*sra [j01]+dxD*sra [j11];
      fq0	=dX*(srq[j10]-srq[j00])+dXD*sraq[j00]+dxD*sraq[j10];
      fq1	=dX*(srq[j11]-srq[j01])+dXD*sraq[j01]+dxD*sraq[j11];
      ra[k]	=YF*f0+yF*f1+YD*fq0+yD*fq1;
      
      /* z,z'_a,z'_q */
      f0	=XF*sz [j00]+xF*sz [j10]+XD*sza [j00]+xD*sza [j10];
      f1	=XF*sz [j01]+xF*sz [j11]+XD*sza [j01]+xD*sza [j11];
      fq0	=XF*szq[j00]+xF*szq[j10]+XD*szaq[j00]+xD*szaq[j10];
      fq1	=XF*szq[j01]+xF*szq[j11]+XD*szaq[j01]+xD*szaq[j11];
      z[k]=YF*f0+yF*f1+YD*fq0+yD*fq1;
      zq[k]=dY*(f1-f0)+dYD*fq0+dyD*fq1;
      f0	=dX*(sz [j10]-sz [j00])+dXD*sza [j00]+dxD*sza [j10];
      f1	=dX*(sz [j11]-sz [j01])+dXD*sza [j01]+dxD*sza [j11];
      fq0	=dX*(szq[j10]-szq[j00])+dXD*szaq[j00]+dxD*szaq[j10];
      fq1	=dX*(szq[j11]-szq[j01])+dXD*szaq[j01]+dxD*szaq[j11];
      za[k]=YF*f0+yF*f1+YD*fq0+yD*fq1;
      
      /* B,B'_a,B'_q */
      f0	=XF*aB [j00]+xF*aB [j10]+XD*aBa [j00]+xD*aBa [j10];
      f1	=XF*aB [j01]+xF*aB [j11]+XD*aBa [j01]+xD*aBa [j11];
      fq0	=XF*aBq[j00]+xF*aBq[j10]+XD*aBaq[j00]+xD*aBaq[j10];
      fq1	=XF*aBq[j01]+xF*aBq[j11]+XD*aBaq[j01]+xD*aBaq[j11];
      B[k]=YF*f0+yF*f1+YD*fq0+yD*fq1;
      Bq[k]=dY*(f1-f0)+dYD*fq0+dyD*fq1;
      f0	=dX*(aB [j10]-aB [j00])+dXD*aBa [j00]+dxD*aBa [j10];
      f1	=dX*(aB [j11]-aB [j01])+dXD*aBa [j01]+dxD*aBa [j11];
      fq0	=dX*(aBq[j10]-aBq[j00])+dXD*aBaq[j00]+dxD*aBaq[j10];
      fq1	=dX*(aBq[j11]-aBq[j01])+dXD*aBaq[j01]+dxD*aBaq[j11];
      Ba[k]=YF*f0+yF*f1+YD*fq0+yD*fq1;

      /* gh,gh'_a,gh'_q */
      f0	=XF*gHq [j00]+xF*gHq [j10]+XD*gHaq [j00]+xD*gHaq [j10];
      f1	=XF*gHq [j01]+xF*gHq [j11]+XD*gHaq [j01]+xD*gHaq [j11];
      fq0	=XF*gHqq[j00]+xF*gHqq[j10]+XD*gHaqq[j00]+xD*gHaqq[j10];
      fq1	=XF*gHqq[j01]+xF*gHqq[j11]+XD*gHaqq[j01]+xD*gHaqq[j11];
      ghq[k]	=YF*f0+yF*f1+YD*fq0+yD*fq1;
      gh[k]	=XF*gH[j00]+xF*gH[j10]+XD*gHa[j00]+xD*gHa[j10]+
	hq*0.5*((1.+YY*(YY-2.*Y))*f0+yy*(2.*y-yy)*f1
		+cr6*hq*((1.-YY*(YY+4.*Y*y))*fq0-yy*(4.*y*Y+yy)*fq1));
      f0	=dX*(gHq [j10]-gHq [j00])+dXD*gHaq [j00]+dxD*gHaq [j10];
      f1	=dX*(gHq [j11]-gHq [j01])+dXD*gHaq [j01]+dxD*gHaq [j11];
      fq0	=dX*(gHqq[j10]-gHqq[j00])+dXD*gHaqq[j00]+dxD*gHaqq[j10];
      fq1	=dX*(gHqq[j11]-gHqq[j01])+dXD*gHaqq[j01]+dxD*gHaqq[j11];
      gha[k]	=dX*(gH[j10]-gH[j00])+dXD*gHa[j00]+dxD*gHa[j10]+
	hq*0.5*((1.+YY*(YY-2.*Y))*f0+yy*(2.*y-yy)*f1
		+cr6*hq*((1.-YY*(YY+4.*Y*y))*fq0-yy*(4.*y*Y+yy)*fq1));

      F[k]	=XF*aF[i0]+xF*aF[i1]+XD*daF[i0]+xD*daF[i1];
      Fa[k]	=dX*(aF[i1]-aF[i0])+dXD*daF[i0]+dxD*daF[i1];
      gFa[k]	=XF*dgF[i0]+xF*dgF[i1]+XD*d2gF[i0]+xD*d2gF[i1];
      gFaa[k]	=dX*(dgF[i1]-dgF[i0])+dXD*d2gF[i0]+dxD*d2gF[i1];
      gFaa[k]	=a*gFaa[k]+gFa[k];
      gYa[k]	=XF*dgY[i0]+xF*dgY[i1]+XD*d2gY[i0]+xD*d2gY[i1];
      gYaa[k]	=dX*(dgY[i1]-dgY[i0])+dXD*d2gY[i0]+dxD*d2gY[i1];
      gYaa[k]	=a*gYaa[k]+gYa[k];
      T[k]	=XF*aT[i0]+xF*aT[i1]+XD*daT[i0]+xD*daT[i1];
      Ta[k]	=dX*(aT[i1]-aT[i0])+dXD*daT[i0]+dxD*daT[i1];
      P[k]	=XF*aP[i0]+xF*aP[i1]+XD*daP[i0]+xD*daP[i1];
      Pa[k]	=dX*(aP[i1]-aP[i0])+dXD*daP[i0]+dxD*daP[i1];
      if(Ne != NULL){
	Ne[k]	=(XF*sne[i0]+xF*sne[i1]+XD*dsne[i0]+xD*dsne[i1])*raT;
	dNe[k]	=(dX*(sne[i1]-sne[i0])+dXD*dsne[i0]+dxD*dsne[i1])*raT;
      }
      gFa[k]	*=a;
      gYa[k]	*=a;
    }
    else{	/* a == 0. */
      r[k]	=sr[0];
      ra[k]	=YF*sra[j0]+yF*sra[j1]+YD*sraq[j0]+yD*sraq[j1];
      rq[k]	=dY*(sra[j1]-sra[j0])+dYD*sraq[j0]+dyD*sraq[j1];
      z[k]	=sz[0];
      za[k]	=YF*sza[j0]+yF*sza[j1]+YD*szaq[j0]+yD*szaq[j1];
      zq[k]	=dY*(sza[j1]-sza[j0])+dYD*szaq[j0]+dyD*szaq[j1];
      B[k]	=aB[0];
      Ba[k]	=YF*aBa[j0]+yF*aBa[j1]+YD*aBaq[j0]+yD*aBaq[j1];
      Bq[k]	=dY*(aBa[j1]-aBa[j0])+dYD*aBaq[j0]+dyD*aBaq[j1];
      ghq[k]	=YF*gHaq[j0]+yF*gHaq[j1]+YD*gHaqq[j0]+yD*gHaqq[j1];
      gh[k]	=0.;
      gha[k]	=gHa[j0]+hq*0.5*((1.+YY*(YY-2.*Y))*gHaq[j0]
				 +yy*(2.*y-yy)*gHaq[j1]
				 +cr6*hq*((1.-YY*(YY+4.*Y*y))*gHaqq[j0]
					   -yy*(4.*y*Y+yy)*gHaqq[j1]));
      F[k]	=aF[0];
      Fa[k]	=0.;
      gFa[k]	=dgF[0];
      gFaa[k]	=dgF[0];
      gYa[k]	=dgY[0];
      gYaa[k]	=dgY[0];
      T[k]	=aT[0];
      Ta[k]	=0.;
      P[k]	=aP[0];
      Pa[k]	=0.;
      if(Ne != NULL){
	Ne[k]=(XF*sne[i0]+xF*sne[i1]+XD*dsne[i0]+xD*dsne[i1])*raT;
	dNe[k]=(dX*(sne[i1]-sne[i0])+dXD*dsne[i0]+dxD*dsne[i1])*raT;
      }
    }
  }
  return(0);
}

int esigetmagfluxes_(REAL *gf,REAL *gy,REAL *a0)
{
  REAL x,a;
  REAL f0,f1,df0,df1,d2f,d3f;

  a	=*a0;
  if(a < 0.){
    a	=-a;
  }
  i0	=(a-sa[0])*rha;
  if(a != 0.){
    if(i0 >= Na){
      i0	=Na-1;
    }
    i1	=i0+1;
    x	=(a-sa[i0])*rha;
    f0	=dgF[i0];
    f1	=dgF[i1];
    df0	=ha*d2gF[i0];
    df1	=ha*d2gF[i1];
    d2f	=3.*(f1-f0)-2.*df0-df1;
    d3f	=df0+df1-2.*(f1-f0);
    *gf	=gF[i0]+ha*x*(sa[i0]*(f0+(0.5*df0+(cr3*d2f+0.25*d3f*x)*x)*x)
		      +ha*x*(0.5*f0+(cr3*df0+(0.25*d2f+0.2*d3f*x)*x)*x));
    f0	=dgY[i0];
    f1	=dgY[i1];
    df0	=ha*d2gY[i0];
    df1	=ha*d2gY[i1];
    d2f	=3.*(f1-f0)-2.*df0-df1;
    d3f	=df0+df1-2.*(f1-f0);
    *gy	=gY[i0]+ha*x*(sa[i0]*(f0+(0.5*df0+(cr3*d2f+0.25*d3f*x)*x)*x)
		      +ha*x*(0.5*f0+(cr3*df0+(0.25*d2f+0.2*d3f*x)*x)*x));
  }
  else{/* a == 0. */
    *gf	=0.;
    *gy	=gY[0];
  }
  *gf	*=c2gp;
  *gy	*=c2gp;
  return(0);
}

int esigetrzb_(REAL *R, REAL *Z, REAL *Bm,REAL *a0,REAL *q0,int *n)
{
  int j0,j1,k;
  REAL x,X,xx,XX;
  REAL y,Y,yy,YY;
  REAL f0,f1,fq0,fq1;
  REAL a,q;

  for(k =0; k < *n; k++){
    a	=a0[k];
    q	=q0[k];
    
    if(a < 0.){
      a	=-a;
    }
    i0	=(a-sa[0])*rha;
    if(i0 >= Na){
      i0	=Na-1;
    }
    i1	=i0+1;
    
    x	=(a-sa[i0])*rha;
    X	=1.-x;
    
    XX	=X*X;
    xx	=x*x;
    XF	=XX*(3.-2.*X);
    XD	=XX*x*ha;
    xF	=xx*(3.-2.*x);
    xD	=-xx*X*ha;
    
    dX	=3.*x*X;
    dXD	=X-dX;
    dxD	=x-dX;
    dX	*=2.*rha;
    
    j1	=q*crgq0;
    if(q < 0.){
      j1--;
    }
    q	-=cgq0*j1;
    j0	=q*rhq;
    if(j0 >= Np){
      j0	-=Np;
      q	-=cgq0;
    }
    if(j0 < 0){
      j0	+=Np;
      q	+=cgq0;
    }
    j1	=j0+1;
    
    Y	=(gq[j1]-q)*rhq;
    y	=1.-Y;
    
    YY	=Y*Y;
    yy	=y*y;
    YF	=YY*(3.-2.*Y);
    YD	=YY*y*hq;
    yF	=yy*(3.-2.*y);
    yD	=-yy*Y*hq;
    
    dY	=3.*y*Y;
    dYD	=Y-dY;
    dyD	=y-dY;
    dY	*=2.*rhq;
    
    j00	=Np1*i0+j0;
    j01	=j00+1;
    j10	=j00+Np1;
    j11	=j10+1;
    
    f0	=XF*sr [j00]+xF*sr [j10]+XD*sra [j00]+xD*sra [j10];
    f1	=XF*sr [j01]+xF*sr [j11]+XD*sra [j01]+xD*sra [j11];
    fq0	=XF*srq[j00]+xF*srq[j10]+XD*sraq[j00]+xD*sraq[j10];
    fq1	=XF*srq[j01]+xF*srq[j11]+XD*sraq[j01]+xD*sraq[j11];
    R[k]=YF*f0+yF*f1+YD*fq0+yD*fq1;
    
    f0	=XF*sz [j00]+xF*sz [j10]+XD*sza [j00]+xD*sza [j10];
    f1	=XF*sz [j01]+xF*sz [j11]+XD*sza [j01]+xD*sza [j11];
    fq0	=XF*szq[j00]+xF*szq[j10]+XD*szaq[j00]+xD*szaq[j10];
    fq1	=XF*szq[j01]+xF*szq[j11]+XD*szaq[j01]+xD*szaq[j11];
    Z[k]=YF*f0+yF*f1+YD*fq0+yD*fq1;
    
    f0	=XF*aB [j00]+xF*aB [j10]+XD*aBa [j00]+xD*aBa [j10];
    f1	=XF*aB [j01]+xF*aB [j11]+XD*aBa [j01]+xD*aBa [j11];
    fq0	=XF*aBq[j00]+xF*aBq[j10]+XD*aBaq[j00]+xD*aBaq[j10];
    fq1	=XF*aBq[j01]+xF*aBq[j11]+XD*aBaq[j01]+xD*aBaq[j11];
    Bm[k]=YF*f0+yF*f1+YD*fq0+yD*fq1;
  }
  return(0);
}

int gcmotion_(REAL *dXgr,REAL *dXa,REAL *dXgq,REAL *dXgf
	      ,REAL *Xgr,REAL *Xa,REAL *Xgq,REAL *Xgm,int*n)
{
  int j0,j1,k;
  REAL r,ra,rq,za,zq,B,Ba,Bq,F,gYa,T,P;
  REAL a,q,a11,a12,a21,a22;

  REAL x,X,xx,XX;
  REAL f0,f1,fq0,fq1;
  REAL y,Y,yy,YY;
  REAL D;
  for(k=0; k < *n; k++){
    a	=Xa[k];
    if(a < 1e-100){
      a	=0.;
    }
    q	=Xgq[k];

    j1	=q*crgq0;
    if(q < 0.){
      j1--;
    }
    q	-=cgq0*j1;
    j0	=q*rhq;
    if(j0 >= Np){
      j0	-=Np;
      q	-=cgq0;
    }
    if(j0 < 0){
      j0	+=Np;
      q	+=cgq0;
    }
    j1	=j0+1;
    y	=(q-gq[j0])*rhq;

    Y	=1.-y;
    YY	=Y*Y;
    yy	=y*y;
    YF	=YY*(3.-2.*Y);
    YD	=YY*y*hq;
    yF	=yy*(3.-2.*y);
    yD	=-yy*Y*hq;
    dY	=3.*y*Y;
    dYD	=Y-dY;
    dyD	=y-dY;
    dY	*=2.*rhq;

    if(a < 0.){
      a	=-a;
    }
    i0	=(a-sa[0])*rha;
    if(a != 0.){
      if(i0 >= Na){
	i0	=Na-1;
      }
      i1	=i0+1;
      x		=(a-sa[i0])*rha;
      X		=1.-x;
      XX	=X*X;
      xx	=x*x;
      XF	=XX*(3.-2.*X);
      XD	=XX*x*ha;
      xF	=xx*(3.-2.*x);
      xD	=-xx*X*ha;
      
      dX	=3.*x*X;
      dXD	=X-dX;
      dxD	=x-dX;
      dX	*=2.*rha;
      
      j00	=Np1*i0+j0;
      j01	=j00+1;
      j10	=j00+Np1;
      j11	=j10+1;
      
      /* r,r'_a,r'_q */
      f0	=XF*sr [j00]+xF*sr [j10]+XD*sra [j00]+xD*sra [j10];
      f1	=XF*sr [j01]+xF*sr [j11]+XD*sra [j01]+xD*sra [j11];
      fq0	=XF*srq[j00]+xF*srq[j10]+XD*sraq[j00]+xD*sraq[j10];
      fq1	=XF*srq[j01]+xF*srq[j11]+XD*sraq[j01]+xD*sraq[j11];
      r		=YF*f0+yF*f1+YD*fq0+yD*fq1;
      rq	=dY*(f1-f0)+dYD*fq0+dyD*fq1;
      f0	=dX*(sr [j10]-sr [j00])+dXD*sra [j00]+dxD*sra [j10];
      f1	=dX*(sr [j11]-sr [j01])+dXD*sra [j01]+dxD*sra [j11];
      fq0	=dX*(srq[j10]-srq[j00])+dXD*sraq[j00]+dxD*sraq[j10];
      fq1	=dX*(srq[j11]-srq[j01])+dXD*sraq[j01]+dxD*sraq[j11];
      ra	=YF*f0+yF*f1+YD*fq0+yD*fq1;
      
      /* z,z'_a,z'_q */
      f0	=XF*sz [j00]+xF*sz [j10]+XD*sza [j00]+xD*sza [j10];
      f1	=XF*sz [j01]+xF*sz [j11]+XD*sza [j01]+xD*sza [j11];
      fq0	=XF*szq[j00]+xF*szq[j10]+XD*szaq[j00]+xD*szaq[j10];
      fq1	=XF*szq[j01]+xF*szq[j11]+XD*szaq[j01]+xD*szaq[j11];
      zq	=dY*(f1-f0)+dYD*fq0+dyD*fq1;
      f0	=dX*(sz [j10]-sz [j00])+dXD*sza [j00]+dxD*sza [j10];
      f1	=dX*(sz [j11]-sz [j01])+dXD*sza [j01]+dxD*sza [j11];
      fq0	=dX*(szq[j10]-szq[j00])+dXD*szaq[j00]+dxD*szaq[j10];
      fq1	=dX*(szq[j11]-szq[j01])+dXD*szaq[j01]+dxD*szaq[j11];
      za	=YF*f0+yF*f1+YD*fq0+yD*fq1;
      
      /* B,B'_a,B'_q */
      f0	=XF*aB [j00]+xF*aB [j10]+XD*aBa [j00]+xD*aBa [j10];
      f1	=XF*aB [j01]+xF*aB [j11]+XD*aBa [j01]+xD*aBa [j11];
      fq0	=XF*aBq[j00]+xF*aBq[j10]+XD*aBaq[j00]+xD*aBaq[j10];
      fq1	=XF*aBq[j01]+xF*aBq[j11]+XD*aBaq[j01]+xD*aBaq[j11];
      B		=YF*f0+yF*f1+YD*fq0+yD*fq1;
      Bq	=dY*(f1-f0)+dYD*fq0+dyD*fq1;
      f0	=dX*(aB [j10]-aB [j00])+dXD*aBa [j00]+dxD*aBa [j10];
      f1	=dX*(aB [j11]-aB [j01])+dXD*aBa [j01]+dxD*aBa [j11];
      fq0	=dX*(aBq[j10]-aBq[j00])+dXD*aBaq[j00]+dxD*aBaq[j10];
      fq1	=dX*(aBq[j11]-aBq[j01])+dXD*aBaq[j01]+dxD*aBaq[j11];
      Ba	=YF*f0+yF*f1+YD*fq0+yD*fq1;

      F		=XF*aF[i0]+xF*aF[i1]+XD*daF[i0]+xD*daF[i1];
      gYa	=XF*dgY[i0]+xF*dgY[i1]+XD*d2gY[i0]+xD*d2gY[i1];
      T		=XF*aT[i0]+xF*aT[i1]+XD*daT[i0]+xD*daT[i1];
      P		=XF*aP[i0]+xF*aP[i1]+XD*daP[i0]+xD*daP[i1];
      gYa	*=a;
    }
    else{	/* a == 0. */
      r		=sr[0];
      ra	=YF*sra[j0]+yF*sra[j1]+YD*sraq[j0]+yD*sraq[j1];
      rq	=dY*(sra[j1]-sra[j0])+dYD*sraq[j0]+dyD*sraq[j1];
      za	=YF*sza[j0]+yF*sza[j1]+YD*szaq[j0]+yD*szaq[j1];
      zq	=dY*(sza[j1]-sza[j0])+dYD*szaq[j0]+dyD*szaq[j1];
      B		=aB[0];
      Ba	=YF*aBa[j0]+yF*aBa[j1]+YD*aBaq[j0]+yD*aBaq[j1];
      Bq	=dY*(aBa[j1]-aBa[j0])+dYD*aBaq[j0]+dyD*aBaq[j1];
      F		=aF[0];
      gYa	=dgY[0];
      T		=aT[0];
      P		=aP[0];
      
      T		*=Xgr[k];
      D		=za*rq-zq*ra;
      a21	=1./(D*((F+T)/r+Xgr[k]*r*P));
#ifdef H
      r  =cgq0*Xgr[k]*B;
      B	*=r;
      r	=r*Xgr[k]+cgq0*Xgm[k];
#endif
      r  =Xgr[k]*B;
      B	*=r;
      r	=r*Xgr[k]+Xgm[k];

      Ba	*=r;
      Bq	*=-r;
      dXgr[k]	=0.;
      Bq	*=a21;
      Ba	*=a21;
      dXgf[k]	=B/F;
      T	=cos(Xgq[k]);
      P	=sin(Xgq[k]);
      dXa[k]	=Bq*T-Ba*P;
      dXgq[k]	=Bq*P+Ba*T;
    }
    a21	=F;
    T	*=Xgr[k];
    D	=za*rq-zq*ra;
    a12	=D*((F+T)/r+Xgr[k]*r*P);
    D	=-gYa/(r*D);
    a11	=(rq*rq+zq*zq)*D;
    D	*=ra*rq+za*zq;
    a22	=gYa*(1.+T/F);
    r	=1./(a11*a22-a12*a21);
    a11	*=r;
    a12	*=r;
    a21	*=r;
    a22	*=r;
    D	*=r;
#ifdef H
    r  =Xgr[k]*B*cgq0;
    B	*=r;
    r	=r*Xgr[k]+Xgm[k]*cgq0;
#endif
    r  =Xgr[k]*B;
    B	*=r;
    r	=r*Xgr[k]+Xgm[k];
    Ba	*=r;
    Bq	*=-r;
    dXgr[k]	=a22*Bq;
    dXgf[k]	=-a12*B+a11*Ba+D*Bq;
    dXa[k]	=-a21*Bq;
    dXgq[k]	=a22*B-a21*Ba;
    if(sw != NULL && sw[k]){
      Bq	=a*dXgq[k];
      Ba	=dXa[k];
      T	=cos(Xgq[k]);
      P	=sin(Xgq[k]);
      dXa[k]	=Ba*T-Bq*P;
      dXgq[k]	=Ba*P+Bq*T;
    }
  }
  return(0);
}

int MFieldLines(REAL *dXgq,REAL *dXgf,REAL *Xa,REAL *Xgq,int n) 
{
  int j0,j1,k;
  REAL r,ra,rq,za,zq,B,Ba,Bq,F,gYa;
  REAL a,q;
  REAL x,X,xx,XX;
  REAL f0,f1,fq0,fq1;
  REAL y,Y,yy,YY;
  REAL D;

  for(k=0; k < n; k++){
    a	=Xa[k];
    if(a < 1e-100){
      a	=0.;
    }
    q	=Xgq[k];

    j1	=q*crgq0;
    if(q < 0.){
      j1--;
    }
    q	-=cgq0*j1;
    j0	=q*rhq;
    if(j0 >= Np){
      j0	-=Np;
      q	-=cgq0;
    }
    if(j0 < 0){
      j0	+=Np;
      q	+=cgq0;
    }
    j1	=j0+1;
    y	=(q-gq[j0])*rhq;

    Y	=1.-y;
    YY	=Y*Y;
    yy	=y*y;
    YF	=YY*(3.-2.*Y);
    YD	=YY*y*hq;
    yF	=yy*(3.-2.*y);
    yD	=-yy*Y*hq;
    dY	=3.*y*Y;
    dYD	=Y-dY;
    dyD	=y-dY;
    dY	*=2.*rhq;

    if(a < 0.){
      a	=-a;
    }
    i0	=(a-sa[0])*rha;
    if(i0 >= Na){
      i0	=Na-1;
    }
    i1	=i0+1;
    x		=(a-sa[i0])*rha;
    X		=1.-x;
    XX	=X*X;
    xx	=x*x;
    XF	=XX*(3.-2.*X);
    XD	=XX*x*ha;
    xF	=xx*(3.-2.*x);
    xD	=-xx*X*ha;
    
    dX	=3.*x*X;
    dXD	=X-dX;
    dxD	=x-dX;
    dX	*=2.*rha;
    
    j00	=Np1*i0+j0;
    j01	=j00+1;
    j10	=j00+Np1;
    j11	=j10+1;
    
    /* r,r'_a,r'_q */
    f0	=XF*sr [j00]+xF*sr [j10]+XD*sra [j00]+xD*sra [j10];
    f1	=XF*sr [j01]+xF*sr [j11]+XD*sra [j01]+xD*sra [j11];
    fq0	=XF*srq[j00]+xF*srq[j10]+XD*sraq[j00]+xD*sraq[j10];
    fq1	=XF*srq[j01]+xF*srq[j11]+XD*sraq[j01]+xD*sraq[j11];
    r		=YF*f0+yF*f1+YD*fq0+yD*fq1;
    rq	=dY*(f1-f0)+dYD*fq0+dyD*fq1;
    f0	=dX*(sr [j10]-sr [j00])+dXD*sra [j00]+dxD*sra [j10];
    f1	=dX*(sr [j11]-sr [j01])+dXD*sra [j01]+dxD*sra [j11];
    fq0	=dX*(srq[j10]-srq[j00])+dXD*sraq[j00]+dxD*sraq[j10];
    fq1	=dX*(srq[j11]-srq[j01])+dXD*sraq[j01]+dxD*sraq[j11];
    ra	=YF*f0+yF*f1+YD*fq0+yD*fq1;
    
    /* z,z'_a,z'_q */
    f0	=XF*sz [j00]+xF*sz [j10]+XD*sza [j00]+xD*sza [j10];
    f1	=XF*sz [j01]+xF*sz [j11]+XD*sza [j01]+xD*sza [j11];
    fq0	=XF*szq[j00]+xF*szq[j10]+XD*szaq[j00]+xD*szaq[j10];
    fq1	=XF*szq[j01]+xF*szq[j11]+XD*szaq[j01]+xD*szaq[j11];
    zq	=dY*(f1-f0)+dYD*fq0+dyD*fq1;
    f0	=dX*(sz [j10]-sz [j00])+dXD*sza [j00]+dxD*sza [j10];
    f1	=dX*(sz [j11]-sz [j01])+dXD*sza [j01]+dxD*sza [j11];
    fq0	=dX*(szq[j10]-szq[j00])+dXD*szaq[j00]+dxD*szaq[j10];
    fq1	=dX*(szq[j11]-szq[j01])+dXD*szaq[j01]+dxD*szaq[j11];
    za	=YF*f0+yF*f1+YD*fq0+yD*fq1;
    
    /* B,B'_a,B'_q */
    f0	=XF*aB [j00]+xF*aB [j10]+XD*aBa [j00]+xD*aBa [j10];
    f1	=XF*aB [j01]+xF*aB [j11]+XD*aBa [j01]+xD*aBa [j11];
    fq0	=XF*aBq[j00]+xF*aBq[j10]+XD*aBaq[j00]+xD*aBaq[j10];
    fq1	=XF*aBq[j01]+xF*aBq[j11]+XD*aBaq[j01]+xD*aBaq[j11];
    B		=YF*f0+yF*f1+YD*fq0+yD*fq1;
    Bq	=dY*(f1-f0)+dYD*fq0+dyD*fq1;
    f0	=dX*(aB [j10]-aB [j00])+dXD*aBa [j00]+dxD*aBa [j10];
    f1	=dX*(aB [j11]-aB [j01])+dXD*aBa [j01]+dxD*aBa [j11];
    fq0	=dX*(aBq[j10]-aBq[j00])+dXD*aBaq[j00]+dxD*aBaq[j10];
    fq1	=dX*(aBq[j11]-aBq[j01])+dXD*aBaq[j01]+dxD*aBaq[j11];
    Ba	=YF*f0+yF*f1+YD*fq0+yD*fq1;
    
    F	=XF*aF[i0]+xF*aF[i1]+XD*daF[i0]+xD*daF[i1];
    gYa	=XF*dgY[i0]+xF*dgY[i1]+XD*d2gY[i0]+xD*d2gY[i1];
    gYa	*=a;

    D	=za*rq-zq*ra;
    dXgq[k]	=-gYa/(r*D*B);
    dXgf[k]	=F/(r*r*B);
  }
  return(0);
}

void esigetpressure_(double *p,double *dp,double *a, int *n)
{
  int k,i;
  double p0[Na1],dp0,dp1,d2p0,d2p1,x,X,p3x,p3X,XF,XD,xF,xD;

  p0[0]	=0.;
  dp1	=0.;
  d2p1	=aP[0]*dgY[0];
  for(i=1,k=0; i < Na1; i++,k++){
    dp0	=dp1;
    d2p0=d2p1;
    dp1	=aP[i]*dgY[i]*sa[i];
    d2p1=aP[i]*dgY[i]+(daP[i]*dgY[i]+aP[i]*d2gY[i])*sa[i];
    p0[i]=p0[k]+0.5*ha*((dp0+dp1)+ha*(d2p0-d2p1)/6.);
  }
  for(i=0; i < Na1; i++){
    p0[i]	-=p0[Na];
  }
  for(k=0; k < *n; k++){
    i	=(a[k]-sa[0])*rha;
    if(i >= Na){
      i	=Na-1;
    }
    x	=(a[k]-sa[i])/ha;
    X	=1.-x;
    XF	=X*X*(3.-2.*X);
    XD	=X*X*x*ha;
    xF	=x*x*(3.-2.*x);
    xD	=-x*x*X*ha;
    p3x	=x*x*x;
    p3X	=X*X*X;
    p[k]=p0[i]+0.5*ha*((1.+p3X*(X-2.))*aP[i]+p3x*(2.-x)*aP[i+1]
		       +ha*((1.-p3X*(X+4.*x))*daP[i]
			    -p3x*(4.*X+x)*daP[i+1])/6.);
    dp[k]	=(XF*aP[i]+xF*aP[i+1]+XD*daP[i]+xD*daP[i+1])*
      (XF*dgY[i]+xF*dgY[i+1]+XD*d2gY[i]+xD*d2gY[i+1])*a[k];
  }
  return;
}

static double cme=9.1094e-28; /* [g] - electron mass */
static double cmi=1.6726e-24; /* [g] - proton mass */
static double cc=2.9979e+10; /* [cm/sec] - speed of light */
static double ce=4.8032e-10; /* [CGS] - proton electric charge */
static double cEeV=1.6022e-12; /* erg energy of 1 eV */
static double cEkeV=1.6022e-9; /* erg energy of 1 keV */

static char tga,tgb;
static double sega,smga,gmga,Zga,pZei,segb,smgb,gmgb,Zgb;

void esigetdensity_(double *ne,double *dne,double *Te,double *Ti,double *a)
{
  int k,i;
  double s,x,X,XF,XD,xF,xD;

  i	=(*a-sa[0])*rha;
  if(i >= Na){
    i	=Na-1;
  }
  x	=(*a-sa[i])/ha;
  X	=1.-x;
  XF	=X*X*(3.-2.*X);
  XD	=X*X*x*ha;
  xF	=x*x*(3.-2.*x);
  xD	=-x*x*X*ha;
  s	=(*Te)+(*Ti) != 0. ? 30./((*Te)+(*Ti)) : 1.;
  *ne	=(XF*sne[i]+xF*sne[i+1]+XD*dsne[i]+xD*dsne[i+1])*s;
  if(dne != NULL){
    double dX,dXD,dxD;
    dX	=3.*x*X;
    dXD	=X-dX;
    dxD	=x-dX;
    dX	*=2.*rha;
    *dne	=(dX*(sne[i+1]-sne[i])+dXD*dsne[i]+dxD*dsne[i+1])*s;
  }
  return;
}

void GetEECollisions(double *gn,double *E,double *Te,double *ne)
{
  int i;
  double gl,t;
  double pva,pvb,v,x,gy;

  pva	=2.*cEkeV*(*E)/smga;
  pvb	=2.*cEkeV*(*Te)/smgb;

  x	=pva/pvb;
  v	=sqrt(x);
  if(v >= tF1){
    gy	=1.-cr2/x;
  }
  else{
    t	=v*rhF;
    i	=t;
    t	-=(double)i;
    gy	=(x-cr2)*(F0[i]+t*(F1[i]+t*(F2[i]+t*F3[i])));
  }
  gy	+=c2rqgp*v*exp(-x);

  gl	=*Te > 0.01 ? 24.-log(sqrt((*ne)*1e+8)/(*Te)) :
    23.-0.5*log((*ne)*1e+5/((*Te)*(*Te)*(*Te)));
  t	=sega*segb/smga;
  *gn	=8.*cgp*t*t*gl*(*ne)*1e+14*gy/(pva*sqrt(pva));
  return;
}

void GetIECollisions(double *gn,double *E,double *Te,double *ne)
{
  int i;
  double gl,t;
  double pva,pvb,v,x,gy;

  pva	=2.*cEkeV*(*E)/smga;
  pvb	=2.*cEkeV*(*Te)/smgb;

  x	=pva/pvb;
  v	=sqrt(x);
  if(v >= tF1){
    gy	=1.-cr2/x;
  }
  else{
    t	=v*rhF;
    i	=t;
    t	-=(double)i;
    gy	=(x-cr2)*(F0[i]+t*(F1[i]+t*(F2[i]+t*F3[i])));
  }
  gy	+=c2rqgp*v*exp(-x);

  if(*Te < 0.01*pZei){
    gl	=23.-0.5*log(pZei*(*ne)*1e+5/((*Te)*(*Te)*(*Te)));
  }
  else{
    gl	=24.-log(sqrt((*ne)*1e+8)/(*Te));
  }
  t	=sega*segb/smga;
  *gn	=8.*cgp*t*t*gl*(*ne)*1e+14*gy/(pva*sqrt(pva));
  return;
}


void GetEPCollisions(double *gn,double *E,double *Te,double *Ti,double *ne)
{
  int i;
  double gl,t;
  double pva,pvb,v,x,gy;

  pva	=2.*cEkeV*(*E)/smga;
  pvb	=2.*cEkeV*(*Te)/smgb;

  x	=pva/pvb;
  v	=sqrt(x);
  if(v >= tF1){
    gy	=1.-cr2/x;
  }
  else{
    t	=v*rhF;
    i	=t;
    t	-=(double)i;
    gy	=(x-cr2)*(F0[i]+t*(F1[i]+t*(F2[i]+t*F3[i])));
  }
  gy	+=c2rqgp*v*exp(-x);

  gl	=*Te > 0.01 ? 24.-log(sqrt((*ne)*1e+8)/(*Te)) :
    23.-0.5*log((*ne)*1e+5/((*Te)*(*Te)*(*Te)));
  t	=sega*segb/smga;
  *gn	=8.*cgp*t*t*gl*(*ne)*1e+14*gy/(pva*sqrt(pva));
  return;
}

#ifdef H

void GetIECollisions(double *gn,double *E,double *Te,double *ne)
{
  int i;
  double gl,t;
  double pva,pvb,v,x,gy;

  pva	=2.*cEkeV*(*E)/smga;
  pvb	=2.*cEkeV*(*Te)/smgb;

  x	=pva/pvb;
  v	=sqrt(x);
  if(v >= tF1){
    gy	=1.-cr2/x;
  }
  else{
    t	=v*rhF;
    i	=t;
    t	-=(double)i;
    gy	=(x-cr2)*(F0[i]+t*(F1[i]+t*(F2[i]+t*F3[i])));
  }
  gy	+=c2rqgp*v*exp(-x);

  if(*Te < 0.01*pZei){
    gl	=23.-0.5*log(pZei*(*ne)*1e+5/((*Te)*(*Te)*(*Te)));
  }
  else{
    gl	=24.-log(sqrt((*ne)*1e+8)/(*Te));
  }
  t	=sega*segb/smga;
  *gn	=8.*cgp*t*t*gl*(*ne)*1e+14*gy/(pva*sqrt(pva));
  return;
}
#endif

void GetIICollisions(double *gn,double *E,double *Ti,double *ne)
{
  int i;
  double gl,t,Tj;
  double pva,pvb,v,x,gy;
  
  Tj	=2.*(*E)/3.;
  pva	=2.*cEkeV*(*E)/smga;
  pvb	=2.*cEkeV*(*Ti)/smgb;

  x	=pva/pvb;
  v	=sqrt(x);
  if(v >= tF1){
    gy	=1.-cr2/x;
  }
  else{
    t	=v*rhF;
    i	=t;
    t	-=(double)i;
    gy	=(x-cr2)*(F0[i]+t*(F1[i]+t*(F2[i]+t*F3[i])));
  }
  gy	+=c2rqgp*v*exp(-x);
  gl	=23.-log(Zga*Zgb*Zgb*(gmga+gmgb)/(gmga*Tj+gmgb*(*Ti))
		 *sqrt((*ne)*1e+5/(*Ti)));
  t	=sega*segb/smga;
  *gn	=8.*cgp*t*t*gl*(*ne)*1e+14*gy/(pva*sqrt(pva));
  return;
}

void getcollisions_(
		       double *gn	/* Collision grequency */
		       ,double *E	/* Test particle Energy [keV] */
		       ,double *Tp	/* background temperature [keV] */
		       ,double *n	/* background density [10^{20}/m^3] */
		       )
{
  switch(tga){
  case 'e':	/* electrons */
    switch(tgb){
    case 'e':			/* electrons */
      GetEECollisions(gn,E,Tp,n);
      break;
    case 'P':			/* Proton */
    case 'D':			/* Deuteron */
    case 'T':			/* Triton */
    case 'A':			/* Alpha */
      pZei		=Zgb*Zgb;
      GetIECollisions(gn,E,Tp,n);
      break;
    }
    break;
  case 'P':	/* Proton */
  case 'D':	/* Deuteron */
  case 'T':	/* Triton */
  case 'A':			/* Alpha */
    switch(tgb){
    case 'e':			/* electrons */
      pZei		=Zga*Zga;
      GetIECollisions(gn,E,Tp,n);
      break;
    case 'P':			/* Proton */
    case 'D':			/* Deuteron */
    case 'T':			/* Triton */
    case 'A':			/* Alpha */
      GetIICollisions(gn,E,Tp,n);
      break;
    }
    break;
  }
  return;
}

int setparticletype_(
		     char *ega /* Type of test particles:
				  e   -electrons;
				  H,P -protons;
				  D   -deutrons;
				  T   -tritons;
				  A   -alphas.
				  */
		     ,char *egb /* Type of test particles:
				   e   -electrons;
				   H,P -protons;
				   D   -deutrons;
				   T   -tritons;
				   A   -alphas.
				   */
		     )
{
  tga	=*ega;
  tgb	=*egb;
  switch(tga){
  case 'E':	/* electrons */
    tga	='e';
    break;
  case 'H':	/* Proton */
  case 'p':	/* Proton */
  case 'h':	/* Proton */
    tga	='P';
    break;
  case 'd':	/* Deuteron */
    tga	='D';
    break;
  case 't':	/* Triton */
    tga	='T';
    break;
  case 'a':	/* Alpha */
    tga	='A';
    break;
  }
  switch(tgb){
  case 'E':	/* electrons */
    tgb	='e';
    break;
  case 'H':	/* Proton */
  case 'p':	/* Proton */
  case 'h':	/* Proton */
    tgb	='P';
    break;
  case 'd':	/* Deuteron */
    tgb	='D';
    break;
  case 't':	/* Triton */
    tgb	='T';
    break;
  case 'a':	/* Alpha */
    tgb	='A';
    break;
  }
  switch(tga){
  case 'e':	/* electrons */
    gmga	=cme/cmi;
    Zga		=-1.;
    break;
  case 'P':	/* Proton */
    gmga	=1.;
    Zga		=1.;
    break;
  case 'D':	/* Deuteron */
    gmga	=2.;
    Zga		=1.;
    break;
  case 'T':	/* Triton */
    gmga	=3.;
    Zga		=1.;
    break;
  case 'A':	/* Alpha */
    gmga	=4.;
    Zga		=2.;
    break;
  default:
    printf("'%c' -wrong type. Use from the set:\n'e','H','P','D','T','A' \n",
	   tga);
    return(1);
  }
  sega	=Zga*ce;
  smga	=gmga*cmi;

  switch(tgb){
  case 'e':	/* electrons */
    gmgb	=cme/cmi;
    Zgb		=-1.;
    break;
  case 'P':	/* Proton */
    gmgb	=1.;
    Zgb		=1.;
    break;
  case 'D':	/* Deuteron */
    gmgb	=2.;
    Zgb		=1.;
    break;
  case 'T':	/* Triton */
    gmgb	=3.;
    Zgb		=1.;
    break;
  case 'A':	/* Alpha */
    gmgb	=4.;
    Zgb		=2.;
    break;
  default:
    printf("'%c' -wrong type. Use from the set:\n'e','H','P','D','T','A' \n",
	   tgb);
    return(1);
  }
  segb	=Zgb*ce;
  smgb	=gmgb*cmi;

#ifdef H
  switch(tga){
  case 'e':	/* electrons */
    switch(tgb){
    case 'e':			/* electrons */
      getcollisions_	=GetEECollisions;
      break;
    case 'P':			/* Proton */
    case 'D':			/* Deuteron */
    case 'T':			/* Triton */
    case 'A':			/* Alpha */
      pZei		=Zgb*Zgb;
      getcollisions_	=GetIECollisions;
      break;
    }
    break;
  case 'P':	/* Proton */
  case 'D':	/* Deuteron */
  case 'T':	/* Triton */
  case 'A':			/* Alpha */
    switch(tgb){
    case 'e':			/* electrons */
      pZei		=Zga*Zga;
      getcollisions_	=GetIECollisions;
      break;
    case 'P':			/* Proton */
    case 'D':			/* Deuteron */
    case 'T':			/* Triton */
    case 'A':			/* Alpha */
      getcollisions_	=GetIICollisions;
      break;
    }
    break;
  }
#endif
  return(0);
}
#ifdef H
        XI
    0.000000000E+00    0.500000007E-01    0.100000001E+00    0.150000006E+00    0.200000003E+00
    0.250000000E+00    0.300000012E+00    0.349999994E+00    0.400000006E+00    0.449999988E+00
    0.500000000E+00    0.550000012E+00    0.600000024E+00    0.649999976E+00    0.699999988E+00
    0.750000000E+00    0.800000012E+00    0.850000024E+00    0.899999976E+00    0.949999988E+00
    0.100000000E+01
        Z_i
    0.100000000E+01    0.100000000E+01    0.100000000E+01    0.100000000E+01    0.100000000E+01
    0.100000000E+01    0.100000000E+01    0.100000000E+01    0.100000000E+01    0.100000000E+01
    0.100000000E+01    0.100000000E+01    0.100000000E+01    0.100000000E+01    0.100000000E+01
    0.100000000E+01    0.100000000E+01    0.100000000E+01    0.100000000E+01    0.100000000E+01
    0.100000000E+01
        A_i
    0.199892085E+01    0.199892085E+01    0.199892085E+01    0.199892085E+01    0.199892085E+01
    0.199892085E+01    0.199892085E+01    0.199892085E+01    0.199892085E+01    0.199892085E+01
    0.199892085E+01    0.199892085E+01    0.199892085E+01    0.199892085E+01    0.199892085E+01
    0.199892085E+01    0.199892085E+01    0.199892085E+01    0.199892085E+01    0.199892085E+01
    0.199892085E+01
        Z_I
    0.822000027E+01    0.822000027E+01    0.822000027E+01    0.822000027E+01    0.822000027E+01
    0.822000027E+01    0.822000027E+01    0.822000027E+01    0.822000027E+01    0.822000027E+01
    0.822000027E+01    0.822000027E+01    0.822000027E+01    0.822000027E+01    0.822000027E+01
    0.822000027E+01    0.822000027E+01    0.822000027E+01    0.822000027E+01    0.822000027E+01
    0.822000027E+01
        A_I
    0.164400005E+02    0.164400005E+02    0.164400005E+02    0.164400005E+02    0.164400005E+02
    0.164400005E+02    0.164400005E+02    0.164400005E+02    0.164400005E+02    0.164400005E+02
    0.164400005E+02    0.164400005E+02    0.164400005E+02    0.164400005E+02    0.164400005E+02
    0.164400005E+02    0.164400005E+02    0.164400005E+02    0.164400005E+02    0.164400005E+02
    0.164400005E+02
        n_e
    0.474103824E+14    0.470898998E+14    0.461185367E+14    0.447322815E+14    0.430942338E+14
    0.411247984E+14    0.390770469E+14    0.370484424E+14    0.352197720E+14    0.334929896E+14
    0.317389296E+14    0.295579733E+14    0.272194161E+14    0.250968193E+14    0.232451033E+14
    0.217424373E+14    0.205131350E+14    0.191996866E+14    0.160626943E+14    0.115832489E+14
    0.665311615E+13
        n_i
    0.353365792E+14    0.350976643E+14    0.343735686E+14    0.333401859E+14    0.321190767E+14
    0.306509729E+14    0.291245573E+14    0.276124764E+14    0.262495386E+14    0.249626332E+14
    0.236553057E+14    0.220297669E+14    0.202868029E+14    0.187048637E+14    0.173248630E+14
    0.162049605E+14    0.152886499E+14    0.143092998E+14    0.119703106E+14    0.863095153E+13
    0.495482128E+13
        n_I
    0.146883246E+13    0.145890938E+13    0.142882819E+13    0.138589967E+13    0.133517722E+13
    0.127418797E+13    0.121076512E+13    0.114792771E+13    0.109126924E+13    0.103775622E+13
    0.983409207E+12    0.915840168E+12    0.843383565E+12    0.777610142E+12    0.720223861E+12
    0.673658950E+12    0.635582114E+12    0.594937551E+12    0.497856880E+12    0.359160255E+12
    0.206605210E+12
        t_e
    0.597483691E+04    0.597050830E+04    0.580135400E+04    0.544564893E+04    0.507395898E+04
    0.460544434E+04    0.411372021E+04    0.365850415E+04    0.323481763E+04    0.284544556E+04
    0.249181152E+04    0.217255103E+04    0.188342664E+04    0.161730725E+04    0.137300513E+04
    0.114599939E+04    0.931960815E+03    0.728864441E+03    0.575167908E+03    0.454757294E+03
    0.339626434E+03
        t_i
    0.149032422E+05    0.149145977E+05    0.142802266E+05    0.130923320E+05    0.119312598E+05
    0.107058604E+05    0.944299707E+04    0.834354297E+04    0.739170508E+04    0.662283203E+04
    0.593573975E+04    0.527360352E+04    0.463863037E+04    0.403461084E+04    0.346309082E+04
    0.293445679E+04    0.247001392E+04    0.207305908E+04    0.168798596E+04    0.130814783E+04
    0.495946259E+03
        pot
    0.000000000E+00    0.000000000E+00    0.000000000E+00    0.000000000E+00    0.000000000E+00
    0.000000000E+00    0.000000000E+00    0.000000000E+00    0.000000000E+00    0.000000000E+00
    0.000000000E+00    0.000000000E+00    0.000000000E+00    0.000000000E+00    0.000000000E+00
    0.000000000E+00    0.000000000E+00    0.000000000E+00    0.000000000E+00    0.000000000E+00
    0.000000000E+00
#endif
