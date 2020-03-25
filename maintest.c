#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* FMG solver for the 3d poisson model problem         */

#define pi  3.1415926535897931

typedef struct
{
double hx;             /* mesh size x                          */
double hy;             /* mesh size y                          */
double hz;			   /* mesh size z                          */
int    ii;             /* number of nodes x                    */
int    jj;             /* number of nodes y                    */
int    kk;             /* number of nodes y                    */
double ***u, ***f;       /* unknown and right hand side values   */
                       /* For convenience:                     */
double ***uconv;        /* extra u for convergence check in FMG */
double ***uold;         /* extra u for use in FAS coarsening    */
} Level;


typedef struct
{
int    nx0;             /* number of nodes coarsest grid     */
int    ny0;             /* number of nodes coarsest grid     */
int    maxlevel;        /* number of grid-levels             */
double xa,xb;           /* begin,end of computational domain */
double ya,yb;           /* begin,end of computational domain */
double za,zb;           /* begin,end of computational domain */
double wu;              /* work unit counter                 */
Level  *Lk ;            /* array of grid levels              */
} Stack;


/**********ROUTINES FOR DATASTRUCTURE**********/


double ***matrix3(int nx,int ny, int nz)
{
 int i,j;

 double ***m=(double ***)calloc(nx+1,sizeof(double**));

 for (i=0;i<=nx;i++)
	 {
	 m[i]=(double **)calloc(ny+1,sizeof(double *));
	 for (j=0;j<=ny;j++)
		{
		m[i][j]=(double *)calloc(nz+1,sizeof(double));
	 }
 }

 return m;
}


void initialize(Stack *U, int nx0, int ny0, int nz0, int maxlevel,
                double xa, double xb, double ya, double yb, double za, double zb)
{
/* initialize values in datastructure */

double  hx,hy,hz;
Level  *L;
int     i,ii,jj,kk;

U->xa=xa;
U->xb=xb;
U->ya=ya;
U->yb=yb;
U->za=za;
U->zb=zb;
U->maxlevel=maxlevel;
U->wu=0.0;

U->Lk=(Level *)calloc(maxlevel+1,sizeof(Level));

hx=(xb-xa)/nx0;
hy=(yb-ya)/ny0;
hz=(zb-za)/nz0;
ii=nx0;
jj=ny0;
kk=nz0;
for (i=1;i<=maxlevel;i++)
  {
  L=U->Lk+i;
  L->hx=hx;
  L->hy=hy;
  L->hz=hz;
  L->ii=ii;
  L->jj=jj;
  L->kk=kk;
  L->f    =matrix3(ii,jj,kk);
  L->u    =matrix3(ii,jj,kk);
  L->uold =matrix3(ii,jj,kk);
  L->uconv=matrix3(ii,jj,kk);
  if (i==maxlevel) printf("\n level%2d ii=%4d jj=%4d\n",i,ii,jj,kk);
  if ((L->f==NULL)||(L->u==NULL)||(L->uold==NULL)||(L->uconv==NULL))
    {
    printf("\nproblem allocating memory\n");
    exit(-1);
    }
  hx*=0.5;
  hy*=0.5;
  hz*=0.5;
  ii*=2;
  jj*=2;
  kk*=2;
  }
}


void finalize(Stack *U, int maxlevel)
{
/* free memory at end of program */

Level  *L;
int     i;
for (i=1;i<=maxlevel;i++)
  {
  L=U->Lk+i;
  free(L->f);
  free(L->u);
  free(L->uold);
  free(L->uconv);
  }
free(U->Lk);
}


/**********FUNCTIONS**********/


double U_a(double x, double y, double z)
{
/* analytical solution */

return(sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z));
}


double U_b(double x, double y, double z)
{
/* take analytical solution as boundary condition */

return(sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z));
}


double U_i(double x, double y, double z)
{
/* initial approximation */

return(0.0);
}


double F_i(double x, double y, double z)
{
/* right hand side function */

return(-8.0*pi*pi*sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z));
}


double Lu(double ***u, double rhx2, double rhy2, double rhz2, int i, int j, int n)
{
/* Computes operator in point i,j */

return(rhx2*(u[i-1][j  ][n  ]-2*u[i][j][n]+u[i+1][j  ][n  ])+
       rhy2*(u[i  ][j-1][n  ]-2*u[i][j][n]+u[i  ][j+1][n  ])+
	   rhz2*(u[i  ][j  ][n-1]-2*u[i][j][n]+u[i  ][j  ][n+1]));
}


/**********SINGLE GRID ROUTINES**********/


void init_uf(Stack *U, double (*u0)(double x, double y, double z),
                       double (*ub)(double x, double y, double z),
                       double (*f0)(double x, double y, double z), int k)

/* initializes u, sets boundary condition for u, and  */
/* initializes right hand side values on grid level k */
{
int    i,j,n,ii,jj,kk;
Level *L;
double x,y,z;

L =U->Lk+k;
ii=L->ii;
jj=L->jj;
kk=L->kk;

for (i=0;i<=ii;i++)
  {
  x=U->xa+i*L->hx;
  for (j=0;j<=jj;j++)
    {
    y=U->ya+j*L->hy;
	for (n=0;n<=kk;n++)
	 {
	 z=U->za+n*L->hz;
    if ((i==0)||(j==0)||(n==0)||(i==ii)||(j==jj)||(n==kk))
      L->u[i][j][n]=ub(x,y,z);
    else
      L->u[i][j][n]=u0(x,y,z);
      L->f[i][j][n]=f0(x,y,z);
    }
  }
}
}


void relax(Stack *U, int k)
{
/* perform Gauss-Seidel relaxation on gridlevel k */

Level *L;
double hx, hy, hz, rhx2, rhy2, rhz2, err, ***u, ***f;
int    i,j,n,ii,jj,kk,od;

L =U->Lk+k;
hx=L->hx;
hy=L->hy;
hz=L->hz;
ii=L->ii;
jj=L->jj;
kk=L->kk;
u =L->u;
f =L->f;

rhx2=1.0/(hx*hx);
rhy2=1.0/(hy*hy);
rhz2=1.0/(hz*hz);

for (i=1;i<=ii-1;i++)
  for (j=1;j<=jj-1;j++)
	  for (n=1;n<=kk-1;n++)
		u[i][j][n]+=(f[i][j][n]-Lu(u,rhx2,rhy2,rhz2,i,j,n))/(-2*rhx2-2*rhy2-2*rhz2);

err=0.0;
for (i=1;i<=ii-1;i++)
  for (j=1;j<=jj-1;j++)
	  for (n=1;n<=kk-1;n++)
		err+=fabs(f[i][j][n]-Lu(u,rhx2,rhy2,rhz2,i,j,n) );

U->wu += exp((U->maxlevel-k)*log(0.25));
printf("\nLevel %d Residual %8.5e Wu %8.5e",k,err/((ii-1)*(jj-1)*(kk-1)),U->wu);
}


double conver_a(Stack *U, int k)
{
/* computes L1 norm of difference between U and analytic solution */
/* on level k                                                     */

Level  *L;
double x,y,z,err,***u;
int    i,j,n;

L =U->Lk+k;
u =L->u;

err=0.0;

for (i=1;i<=L->ii-1;i++)
  {
  x=U->xa+i*L->hx;
  for (j=1;j<=L->jj-1;j++)
    {
    y=U->ya+j*L->hy;
	for (n=1;n<=L->kk-1;n++)
    {
	z=U->za+n*L->hz;
    err+=fabs(u[i][j][n]-U_a(x,y,z));
    }
  }
  }
return(err/((L->ii-1)*(L->jj-1)*(L->kk-1)));
}


/**********INTER GRID ROUTINES**********/


void coarsen_u(Stack *U, int k)
{
/* compute initial approximation on level k-1 */
/* in coarsening step from level k            */

int     ic,jc,nc,iic,jjc,kkc;
Level  *Lc, *L;
double ***uc, ***uco, ***u ;

L =U->Lk+k;
u =L->u;

Lc =U->Lk+k-1;
iic=Lc->ii; jjc=Lc->jj; kkc=Lc->kk;
uc =Lc->u;
uco=Lc->uold;

for (ic=1;ic<=iic-1;ic++)
  for (jc=1;jc<=jjc-1;jc++)
	  for (nc=1;nc<=kkc-1;nc++)
		uc[ic][jc][nc]=0.015625*(4.0*(u[2*ic-1][2*jc+1][2*nc-1]+u[2*ic-1][2*jc+1][2*nc+1]+u[2*ic-1][2*jc-1][2*nc-1]+u[2*ic-1][2*jc-1][2*nc+1]+
                                 u[2*ic+1][2*jc-1][2*nc-1]+u[2*ic+1][2*jc-1][2*nc+1]+u[2*ic+1][2*jc+1][2*nc-1]+u[2*ic+1][2*jc+1][2*nc+1])+
                            1.0*(u[2*ic-1][2*jc][2*nc-1]+u[2*ic+1][2*jc][2*nc-1]+u[2*ic][2*jc-1][2*nc-1]+u[2*ic][2*jc+1][2*nc-1]+
							     u[2*ic-1][2*jc][2*nc  ]+u[2*ic+1][2*jc][2*nc  ]+u[2*ic][2*jc-1][2*nc  ]+u[2*ic][2*jc+1][2*nc  ]+
								 u[2*ic-1][2*jc][2*nc+1]+u[2*ic+1][2*jc][2*nc+1]+u[2*ic][2*jc-1][2*nc+1]+u[2*ic][2*jc+1][2*nc+1])+
                            8.0*(u[2*ic  ][2*jc][2*nc-1]+
								 u[2*ic  ][2*jc][2*nc+1]+
								 u[2*ic-1][2*jc][2*nc]+u[2*ic+1][2*jc][2*nc]+u[2*ic][2*jc-1][2*nc]+u[2*ic][2*jc+1][2*nc])+
						    16.0*(u[2*ic][2*jc][2*nc]));

/* store coarse grid solution in uc0 array */
for (ic=0;ic<=iic;ic++)
  for (jc=0;jc<=jjc;jc++)
	  for (nc=0;nc<=kkc;nc++)
		uco[ic][jc][nc]=uc[ic][jc][nc];
}


void coarsen_f(Stack *U, int k)
{
/* compute coarse grid right hand side on level k-1 */
/* in coarsening step from level k                  */

int     ic,jc,nc,iic,jjc,kkc;
Level  *Lc,*L;
double ***u,***uc,***f,***fc ;
double hx,hy,hz,rh2x,rh2y,rh2z,hxc,hyc,hzc,rh2xc,rh2yc,rh2zc  ;
double rooo,ropo,romo,rpoo,rmoo,rmmo,rmpo,rpmo,rppo;
double roop,ropp,romp,rpop,rmop,rmmp,rmpp,rpmp,rppp;
double room,ropm,romm,rpom,rmom,rmmm,rmpm,rpmm,rppm;

L  =U->Lk+k;
hx =L->hx; hy=L->hy; hz=L->hz;
u  =L->u;
f  =L->f;

Lc  =U->Lk+k-1;
iic =Lc->ii; jjc=Lc->jj; kkc=Lc->kk;
hxc =Lc->hx; hyc=Lc->hy; hzc=Lc->hz;
fc  =Lc->f;
uc  =Lc->u;

rh2xc=1.0/(hxc*hxc);
rh2yc=1.0/(hyc*hyc);
rh2zc=1.0/(hzc*hzc);
rh2x =1.0/(hx*hx);
rh2y =1.0/(hy*hy);
rh2z =1.0/(hz*hz);

for (ic=1;ic<=iic-1;ic++)
  for (jc=1;jc<=jjc-1;jc++)
	  for (nc=1;nc<=kkc-1;nc++)
		{
		rooo=(f[2*ic][2*jc][2*nc] -Lu(u,rh2x,rh2y,rh2z,2*ic,2*jc,2*nc));
		rpoo=(f[2*ic+1][2*jc][2*nc] -Lu(u,rh2x,rh2y,rh2z,2*ic+1,2*jc,2*nc));
		rmoo=(f[2*ic-1][2*jc][2*nc] -Lu(u,rh2x,rh2y,rh2z,2*ic-1,2*jc,2*nc));
		ropo=(f[2*ic][2*jc+1][2*nc] -Lu(u,rh2x,rh2y,rh2z,2*ic,2*jc+1,2*nc));
		romo=(f[2*ic][2*jc-1][2*nc] -Lu(u,rh2x,rh2y,rh2z,2*ic,2*jc-1,2*nc));
		rmmo=(f[2*ic-1][2*jc-1][2*nc] -Lu(u,rh2x,rh2y,rh2z,2*ic-1,2*jc-1,2*nc));
		rpmo=(f[2*ic+1][2*jc-1][2*nc] -Lu(u,rh2x,rh2y,rh2z,2*ic+1,2*jc-1,2*nc));
		rmpo=(f[2*ic-1][2*jc+1][2*nc] -Lu(u,rh2x,rh2y,rh2z,2*ic-1,2*jc+1,2*nc));
		rppo=(f[2*ic+1][2*jc+1][2*nc] -Lu(u,rh2x,rh2y,rh2z,2*ic+1,2*jc+1,2*nc));
		roop=(f[2*ic][2*jc][2*nc+1] -Lu(u,rh2x,rh2y,rh2z,2*ic,2*jc,2*nc+1));
		rpop=(f[2*ic+1][2*jc][2*nc+1] -Lu(u,rh2x,rh2y,rh2z,2*ic+1,2*jc,2*nc+1));
		rmop=(f[2*ic-1][2*jc][2*nc+1] -Lu(u,rh2x,rh2y,rh2z,2*ic-1,2*jc,2*nc+1));
		ropp=(f[2*ic][2*jc+1][2*nc+1] -Lu(u,rh2x,rh2y,rh2z,2*ic,2*jc+1,2*nc+1));
        romp=(f[2*ic][2*jc-1][2*nc+1] -Lu(u,rh2x,rh2y,rh2z,2*ic,2*jc-1,2*nc+1));
        rmmp=(f[2*ic-1][2*jc-1][2*nc+1] -Lu(u,rh2x,rh2y,rh2z,2*ic-1,2*jc-1,2*nc+1));
        rpmp=(f[2*ic+1][2*jc-1][2*nc+1] -Lu(u,rh2x,rh2y,rh2z,2*ic+1,2*jc-1,2*nc+1));
        rmpp=(f[2*ic-1][2*jc+1][2*nc+1] -Lu(u,rh2x,rh2y,rh2z,2*ic-1,2*jc+1,2*nc+1));
        rppp=(f[2*ic+1][2*jc+1][2*nc+1] -Lu(u,rh2x,rh2y,rh2z,2*ic+1,2*jc+1,2*nc+1));
        room=(f[2*ic][2*jc][2*nc-1] -Lu(u,rh2x,rh2y,rh2z,2*ic,2*jc,2*nc-1));
        rpom=(f[2*ic+1][2*jc][2*nc-1] -Lu(u,rh2x,rh2y,rh2z,2*ic+1,2*jc,2*nc-1));
        rmom=(f[2*ic-1][2*jc][2*nc-1] -Lu(u,rh2x,rh2y,rh2z,2*ic-1,2*jc,2*nc-1));
        ropm=(f[2*ic][2*jc+1][2*nc-1] -Lu(u,rh2x,rh2y,rh2z,2*ic,2*jc+1,2*nc-1));
        romm=(f[2*ic][2*jc-1][2*nc-1] -Lu(u,rh2x,rh2y,rh2z,2*ic,2*jc-1,2*nc-1));
        rmmm=(f[2*ic-1][2*jc-1][2*nc-1] -Lu(u,rh2x,rh2y,rh2z,2*ic-1,2*jc-1,2*nc-1));
        rpmm=(f[2*ic+1][2*jc-1][2*nc-1] -Lu(u,rh2x,rh2y,rh2z,2*ic+1,2*jc-1,2*nc-1));
        rmpm=(f[2*ic-1][2*jc+1][2*nc-1] -Lu(u,rh2x,rh2y,rh2z,2*ic-1,2*jc+1,2*nc-1));
        rppm=(f[2*ic+1][2*jc+1][2*nc-1] -Lu(u,rh2x,rh2y,rh2z,2*ic+1,2*jc+1,2*nc-1));

    /* FAS coarse grid right hand side   */
    /* with full weighting of residuals  */

    fc[ic][jc][nc]=Lu(uc,rh2xc,rh2yc,rh2zc,ic,jc,nc)+
                  0.015625*(1.0*(rmom+rpom+romm+ropm+rmoo+rpoo+romo+ropo+rmop+rpop+romp+ropp)+
                  4.0*(rmpm+rmpp+rmmm+rmmp+rpmm+rpmp+rppm+rppp)+
                  8.0*(rpoo+rmoo+romo+ropo+roop+room)+
                  16.0*(rooo));
    }
}


void refine(Stack *U, int k)
{
/* Interpolation and addition of coarse grid correction from grid k-1 */
/* to grid k                                                          */

int     ic,jc,nc,iic,jjc,kkc;
Level  *Lc,*L;
double ***uc,***uco,***u;

L  =U->Lk+k;
u  =L->u  ;

Lc =U->Lk+k-1;
iic=Lc->ii; jjc=Lc->jj; kkc=Lc->kk;
uc =Lc->u ;
uco=Lc->uold ;

for (ic=1;ic<=iic;ic++)
  for (jc=1;jc<=jjc;jc++)
    for (nc=1;nc<=kkc;nc++)
        {
        if (ic<iic) u[2*ic  ][2*jc-1][2*nc  ]+=(uc[ic  ][jc  ][nc  ]-uco[ic  ][jc  ][nc  ]+
                                                uc[ic  ][jc-1][nc  ]-uco[ic  ][jc-1][nc  ])*0.5;

        if (ic<iic) u[2*ic  ][2*jc  ][2*nc-1]+=(uc[ic  ][jc  ][nc  ]-uco[ic  ][jc  ][nc  ]+
                                                uc[ic  ][jc  ][nc-1]-uco[ic  ][jc  ][nc-1])*0.5;

        if (ic<iic) u[2*ic  ][2*jc-1][2*nc-1]+=(uc[ic  ][jc  ][nc  ]-uco[ic  ][jc  ][nc  ]+
                                                uc[ic  ][jc-1][nc  ]-uco[ic  ][jc-1][nc  ]+
                                                uc[ic  ][jc  ][nc-1]-uco[ic  ][jc  ][nc-1]+
                                                uc[ic  ][jc-1][nc-1]-uco[ic  ][jc-1][nc-1])*0.25;

        if (jc<jjc) u[2*ic-1][2*jc  ][2*nc  ]+=(uc[ic  ][jc  ][nc  ]-uco[ic  ][jc  ][nc  ]+
                                                uc[ic-1][jc  ][nc  ]-uco[ic-1][jc  ][nc  ])*0.5;

        if (jc<jjc) u[2*ic-1][2*jc  ][2*nc-1]+=(uc[ic  ][jc  ][nc  ]-uco[ic  ][jc  ][nc  ]+
                                                uc[ic-1][jc  ][nc  ]-uco[ic-1][jc  ][nc  ]+
                                                uc[ic  ][jc  ][nc-1]-uco[ic  ][jc  ][nc-1]+
                                                uc[ic-1][jc  ][nc-1]-uco[ic-1][jc  ][nc-1])*0.25;

        if (nc<kkc) u[2*ic][2*jc][2*nc  ]+=(uc[ic  ][jc  ][nc  ]-uco[ic  ][jc  ][nc  ]);


        if (nc<kkc) u[2*ic-1][2*jc-1][2*nc  ]+=(uc[ic  ][jc  ][nc  ]-uco[ic  ][jc  ][nc  ]+
                                                uc[ic-1][jc  ][nc  ]-uco[ic-1][jc  ][nc  ]+
                                                uc[ic  ][jc-1][nc  ]-uco[ic  ][jc-1][nc  ]+
                                                uc[ic-1][jc-1][nc  ]-uco[ic-1][jc-1][nc  ])*0.25;

        u[2*ic-1][2*jc-1][2*nc-1]+=(uc[ic  ][jc  ][nc  ]-uco[ic  ][jc  ][nc  ]+
                                    uc[ic-1][jc  ][nc  ]-uco[ic-1][jc  ][nc  ]+
                                    uc[ic-1][jc-1][nc  ]-uco[ic-1][jc-1][nc  ]+
                                    uc[ic-1][jc-1][nc-1]-uco[ic-1][jc-1][nc-1]+
                                    uc[ic  ][jc-1][nc  ]-uco[ic  ][jc-1][nc  ]+
                                    uc[ic  ][jc-1][nc-1]-uco[ic  ][jc-1][nc-1]+
                                    uc[ic  ][jc  ][nc-1]-uco[ic  ][jc  ][nc-1]+
                                    uc[ic-1][jc  ][nc-1]-uco[ic  ][jc  ][nc-1])*0.125;
        }
}

void fmg_interpolate(Stack *U, int k)
{
/* interpolation of coarse grid k-1 solution to fine grid k */
/* to serve as first approximation. bi-cubic interpolation  */

int    ic,jc,nc,iic,jjc,kkc,i,j,n,ii,jj,kk;
Level  *Lc,*L;
double ***uc,***u, ***uconv;

L  =U->Lk+k;
ii =L ->ii; jj =L->jj; kk=L->kk;
u  =L->u  ;

Lc    =U->Lk+k-1;
iic   =Lc->ii; jjc=Lc->jj; kkc=Lc->kk;
uc    =Lc->u ;
uconv =Lc->uconv;


/* store grid k-1 solution for later use in convergence check */

for (ic=1;ic<=iic-1;ic++)
  for (jc=1;jc<=jjc-1;jc++)
    for (nc=1;nc<=kkc-1;nc++)
        uconv[ic][jc][nc]=uc[ic][jc][nc];

/* first inject to points coinciding with coarse grid points */

for (ic=1;ic<=iic-1;ic++)
  for (jc=1;jc<=jjc-1;jc++)
    for (nc=1;nc<=kkc-1;nc++)
        u[2*ic][2*jc][2*nc]=uc[ic][jc][nc];

/* interpolate intermediate x direction */
for (i=2;i<=ii-2;i+=2)
  {
  u[i][1][0]=(u[i][0][0]+u[i][2][0])*0.5;
  u[i][jj-1][0]=(u[i][jj-2][0]+u[i][jj][0])*0.5;

  u[i][1][kk]=(u[i][0][kk]+u[i][2][kk])*0.5;
  u[i][jj-1][kk]=(u[i][jj-2][kk]+u[i][jj][kk])*0.5;

  for (j=3;j<=jj-3;j+=2)
    for (n=1;n<=kk-1;n+=1)
    u[i][j][n]=(u[i][j-1][n]+u[i][j+1][n])*0.5;

  }

/* interpolate in x direction */

for (j=1;j<=jj-1;j++)
  {
  u[1][j][0]=(u[0][j][0]+u[2][j][0])*0.5;
  u[ii-1][j][0]=(u[ii-2][j][0]+u[ii][j][0])*0.5;

  u[1][j][kk]=(u[0][j][kk]+u[2][j][kk])*0.5;
  u[ii-1][j][kk]=(u[ii-2][j][kk]+u[ii][j][kk])*0.5;

  for (i=3;i<=ii-3;i+=2)
    for (n=1;n<=kk-1;n+=1)
    u[i][j][n]=(u[i-1][j][n]+u[i+1][j][n])*0.5;

   }

  /* interpolate in z direction */

for (j=1;j<=jj-1;j++)
  {
  u[0][1][n]=(u[0][0][n]+u[0][2][n])*0.5;
  u[0][jj-1][n]=(u[0][jj-2][0]+u[0][jj][n])*0.5;

  u[ii][1][n]=(u[ii][0][n]+u[ii][2][n])*0.5;
  u[ii][jj-1][n]=(u[ii][jj-2][n]+u[ii][jj][n])*0.5;

   }
}


double conver(Stack *U, int k)
{
/* convergence check using converged solutions on level k */
/* and on next coarser grid k-1                           */

int    ic,jc,nc;
Level  *Lc,*L;
double ***uc,***u;
double err;

L =U->Lk+k;

if (k==U->maxlevel) u=L->u;
  else u=L->uconv;

Lc =U->Lk+k-1;
uc =Lc->uconv;

err=0.0;
for (ic=1;ic<=Lc->ii-1;ic++)
  for (jc=1;jc<=Lc->jj-1;jc++)
    for (nc=1;nc<=Lc->kk-1;nc++)
        err+=fabs(uc[ic][jc][nc]-u[2*ic][2*jc][2*nc]);

return(err/((Lc->ii-1)*(Lc->jj-1)*(Lc->kk-1)));
}


/**********MULTIGRID DRIVING ROUTINES**********/


void cycle(Stack *U, int k, int nu0, int nu1, int nu2, int gamma)

/* performs coarse grid correction cycle starting on level k */
/* nu1 pre-relaxations, nu2 postrelaxations, nu0 relaxations */
/* on the coarsest grid, cycleindex gamma=1 for Vcycle,      */
/* gamma=2 for Wcycle                                        */

{
int i,j,n;

if (k==1)
  for (i=1;i<=nu0;i++) relax(U,k);
else
  {
  for (i=1;i<=nu1;i++) relax(U,k);
  coarsen_u(U,k);
  coarsen_f(U,k);
  for (j=1;j<=gamma;j++) cycle(U,k-1,nu0,nu1,nu2,gamma);
  refine(U,k);
  for (i=1;i<=nu2;i++) relax(U,k);
  }
}


void fmg(Stack *U, int k, int nu0, int nu1, int nu2, int gamma, int ncy)
{
/* performs FMG with k levels and ncy cycles per level */

int i,j;

if (U->maxlevel==1)
  for (j=1;j< ncy;j++)
    for (i=1;i<=nu0;i++) relax(U,k);
if (k==1) for (i=1;i<=nu0;i++) relax(U,k);
else
  {
  fmg(U,k-1,nu0,nu1,nu2,gamma,ncy);
  fmg_interpolate(U,k);
  for (j=1;j<=ncy;j++)
    {
    cycle(U,k,nu0,nu1,nu2,gamma);
    printf("\n");
    }
  }
}


/**********MAIN PROGRAM**********/


int main()
{
Stack  U;
int j,maxlev,ncy;

printf("\ngive maxlevel"); scanf("%d",&maxlev);
printf("\ngive ncy     "); scanf("%d",&ncy);

initialize(&U,4,4,4,maxlev,-1.0,1.0,-1.0,1.0,-1.0,1.0);

for (j=maxlev;j>=1;j--) init_uf(&U, U_i, U_b, F_i, j);

fmg(&U,maxlev,10,2,1,1,ncy);

printf("\n\nLevel %d: er=%10.8e\n\n",maxlev, conver_a(&U,maxlev));

if (maxlev>1)
  {
  printf("\n");
  for (j=2;j<=maxlev;j++) printf("\naen(%2d,%2d)=%8.5e",j,j-1,conver(&U,j));
  }

printf("\n");
finalize(&U,maxlev);
}
