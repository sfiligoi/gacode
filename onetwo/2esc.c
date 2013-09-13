#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>
#include <string.h>

#if defined(_SEM_SEMUN_UNDEFINED)
union semun {
  int val;                    /* value for SETVAL */
  struct semid_ds *buf;       /* buffer for IPC_STAT, IPC_SET */
  unsigned short int *array;  /* array for GETALL, SETALL */
  struct seminfo *__buf;      /* buffer for IPC_INFO */
};
#endif

#define Sem0 0
#define Sem1 1

int ESNa,ESNa1,ESNp,ESNp1,ESFp,ESFp1;
double *ESaR0,*ESaZ0;
double ESgt[129];

static int Pid;
static int SemID;
static int ShmID;
static union semun Mysemun;
static struct sembuf bufP	={Sem0,-1,SEM_UNDO};
static struct sembuf bufV	={Sem1, 1,IPC_NOWAIT};
static void *lShm;
static int *lRetFl;
static int kD,kI;

static int *ShMemI;
static double *ShMemT;
static double *ShMemB;
static double *ShMemP;
static double *ShMemE;

int LaunchEsc(char *Opt)
{
  int i,k;
  char elem;
  char *lc,*ls,ln[64];
  double t,dt;

  kD	=sizeof(double);
  kI	=sizeof(int);

  /* Initialize */
  Pid	=getpid();

  /* Init semaphores */
  /* Two semaphores (Empty and Full) are created in one set */
  SemID	=semget((key_t)Pid, 2, 0666|IPC_CREAT);

  /* Init Empty to number of elements in shared memory */
  Mysemun.val	=1;
  semctl(SemID, Sem0, SETVAL, Mysemun);

  /* Init Full to zero, no elements are produced yet */
  Mysemun.val	=0;
  semctl(SemID, Sem1, SETVAL, Mysemun);  

  if(Opt != NULL && *Opt != '\0'){
    sprintf(ln,"cd ESC;Obj/esc2 %d %s&",Pid,Opt);
  }
  else{
    sprintf(ln,"cd ESC;Obj/esc2 %d&",Pid);
  }
  semop(SemID, &bufP, 1);   /*Lock ASTRA*/
  system(ln);
  semop(SemID, &bufV, 1);   /*Open ESC*/   		   
  semop(SemID, &bufP, 1);   /*Lock ASTRA*/   
  
  /* Init Shared memory */
  ShmID	=shmget((key_t)(Pid+1),0, 0666 | IPC_CREAT);
  /* attach shared memory to process */
  lShm=shmat(ShmID,NULL,0);

  ShMemI	=(int*)lShm;
  ShMemT	=(double*)(ShMemI+32);
  ShMemB	=ShMemT+16;
  ShMemP	=ShMemB+384;
  ShMemE	=(double*)(ShMemI+0x400*sizeof(double));

  lRetFl=(int*)lShm;
  lShm	+=kI;

  ESNa	=lRetFl[1];
  ESNa1	=ESNa+1;
  ESNp	=lRetFl[2];
  ESNp1	=ESNp+1;
  ESFp	=lRetFl[3];
  ESFp1	=ESFp+1;
  ESaR0	=ShMemE;
  ESaZ0	=ESaR0+1;

#ifdef H
  ESsr	=ESaZ0+1+15*ESNa1;
  i	=ESNp1*ESNa1;
  ESsz	=ESsr+i;
#endif
  ESISetMem(ShMemE+2+15*ESNa1+2*ESNp1*ESNa1+4*ESNa1*ESFp1,ESNa1,ESNp1);
  dt	=8.*atan(1.)/ESNp;
  for(i=0; i < ESNp1; i++){
    t	=dt*i;
    ESgt[i]	=t;
  }
  semop(SemID, &bufV, 1);   /*Open ESC*/   		   
  return(0);
}

void esc0_()
{
  int k;
  k	=0;
  if(k){
    LaunchEsc("-NoX");
  }
  else{
    LaunchEsc(NULL);
  }
  return;
}

void esc1_()
{
  int k;
  struct shmid_ds Myshmid_ds;

  esifree_();
  semop(SemID, &bufP, 1);   
  k	=4;
  memcpy(lShm,(void*)&k,kI);
  semop(SemID, &bufV, 1);   		   

  semop(SemID, &bufP, 1);   
  /* Remove semaphores */
  semctl(SemID,0, IPC_RMID, Mysemun);
  /* Remove shared memory */
  shmctl(ShmID, IPC_RMID, &Myshmid_ds);
  return;
}

void escin_(double	*pProf	/* p[] - pressure related profile */
	    ,double	*jProf	/* j[] - current related profile */
	    ,double	*aProf	/* a[] - normalized sqrt(Phi) square root 
				   of toroidal flux.
				   If a=NULL (dropped in FORTRAN) - uniform 
				   grid from 0 till 1,
				   otherwise it is considered as an array of 
				   grid points */
	    ,int	*nProf	/* number <= 129 of profile points including \
				   magnetic axis and plasma edge */
	    ,int	*nrhH	/* First digit
				   n - 0 - nomalized radial coordinate (0-1)
				       1 - absolute radial coordinate;
				   r - 0 - b (vertical semiaxis)
				           as the radial coordinate,
				       1 - V (volume),
				       2 - gF (toroidal flux),
				       3 - gY (poloidal flux, DO NOT use),
				   h specifies the meaning
				   of pProf:
				   0 - j_p;
				   1 - P =dp/dpsi;
				   2 -p [MPa].
				   Second digit H specifies the meaning
				   of jProf:
				   0 - j_s [MA/m^2];
				   1- j_|| [MA/m^2]=(j dot B)/(R_0 B grad phi),
				                   R_0 = magnetic axis;
				   2- j_||R_0 [MA/m] = (j dot B)/(B grad phi);
				   3 - T=FF';
				   6 - q;
				   7 - 1/q;
				   8 - \gY;
				   For Example:
				   26 - p[] and q[] profiles are supplied;
				   21 - p[] and j||[] profiles are supplied;
				   0 -  jp[] and js[] profiles are supplied.

				   Possible combinations are limited to:
				   0,1,2,6,7
				   10,11,12,16,17
				   21,22,26,27

				 */
	    ,double	*Rpv	/* R[m]- plasma boundary */
	    ,double	*Zpv	/* Z[m]- plasma boundary */
	    ,int	*Npv	/* number <= 257 of the plasma-vacuum 
				   points.
				   If *Npv >12, the first and last points
				   coincide */
	    ,double	*RBtor	/* RBtor [m Tesla] outside the 
				   plasma */
	    ,double	*Rext	/* Reference major radius [m] */
	    )
{
  int i;
  int k;
  void *lc;

  semop(SemID, &bufP, 1);   /*Lock ASTRA*/
  memcpy((void*)(ShMemI+4),(void*)nrhH,kI);
  memcpy((void*)(ShMemI+5),(void*)Npv,kI);
  memcpy((void*)(ShMemI+6),(void*)nProf,kI);
  memcpy((void*)(ShMemT+1),(void*)RBtor,kD);
  memcpy((void*)(ShMemT+2),(void*)Rext,kD);
  k	=(*Npv)*kD;
  memcpy((void*)ShMemB,(void*)Rpv,k);
  memcpy((void*)ShMemB+k,(void*)Zpv,k);

  k	=(*nProf)*kD;
  if(aProf != NULL){
    memcpy((void*)ShMemP,(void*)aProf,k);
  }
  memcpy((void*)ShMemP+k,(void*)pProf,k);
  memcpy((void*)ShMemP+2*k,(void*)jProf,k+2*sizeof(double));

  semop(SemID, &bufV, 1);   /*Open ESC for return confirmation*/
  semop(SemID, &bufP, 1);   /*Lock ASTRA*/   
  semop(SemID, &bufV, 1);   /*Open ESC*/   		   
  return;
}

void esc_(
	  int *Ffail 		/* Flag of failure:
				   0 - normal operation;
				   1 - some problems;
				   2 - problems near the boundary;
				   4 - problems near the axis;
				   8 - problems in the middle;
				   ....;
				   */
	  ,int *Fjob 		/* Flag for job assignment:
				   0 - nothing special was requested;
				   1 - save the geometry;
				   2 - take the geometry;
				   4 - write ESI data files;
				   ...;
				   */ 
	  ,int *Mpol		/* Working number of Fourier harmonics
				   in $\Psi$ */
	  ,double *sTime	/* time */
	  )
{
  int k;
  double *ld;

  semop(SemID, &bufP, 1);   /*Lock ASTRA*/   
#ifdef H
  memcpy((void*)ShMemT,(void*)sTime,kI);
#endif
  *ShMemT=*sTime;
  memcpy((void*)(ShMemI+8),(void*)Fjob,kI);

  semop(SemID, &bufV, 1);   		   
  semop(SemID, &bufP, 1);   /*Lock ASTRA*/   
  *Mpol	=ShMemI[7];
  *Ffail=ShMemI[9];

  ESIInit();
  semop(SemID, &bufV, 1);	/*Open ESC for output calculations*/
  return;
}
