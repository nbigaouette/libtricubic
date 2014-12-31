#include <stdio.h>
#include <stdlib.h>
#include <tricubic.h>
#include <math.h>

#define MAX 10.0
#define Nt 101

int test1(void) {
  int i;
  double v1,v2;
  double x,y,z;
  double f[8],dfdx[8],dfdy[8],dfdz[8],d2fdxdy[8],d2fdxdz[8],d2fdydz[8],d3fdxdydz[8];
  double a[64];

  for (i=0;i<8;i++) {
    f[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    dfdx[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    dfdy[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    dfdz[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    d2fdxdy[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    d2fdxdz[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    d2fdydz[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    d3fdxdydz[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
  }

  tricubic_get_coeff(a,f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz);

  printf("TESTING F VALUES...\n");
  for (i=0;i<8;i++) {
    tricubic_pointID2xyz(i,&x,&y,&z);
    v1=f[i];
    v2=tricubic_eval(a,x,y,z);
    printf("%i\t%lf\t%lf\tError=%lg\n",i,v1,v2,fabs(v1-v2));
    if (fabs(v1-v2)>1.0e-10)
      return(1);
  }
  printf("TESTING DFDX VALUES...\n");
  for (i=0;i<8;i++) {
    tricubic_pointID2xyz(i,&x,&y,&z);
    v1=dfdx[i];
    v2=tricubic_eval(a,x,y,z,1,0,0);
    printf("%i\t%lf\t%lf\tError=%lg\n",i,v1,v2,fabs(v1-v2));
    if (fabs(v1-v2)>1.0e-10)
      return(1);
  }
  printf("TESTING DFDY VALUES...\n");
  for (i=0;i<8;i++) {
    tricubic_pointID2xyz(i,&x,&y,&z);
    v1=dfdy[i];
    v2=tricubic_eval(a,x,y,z,0,1,0);
    printf("%i\t%lf\t%lf\tError=%lg\n",i,v1,v2,fabs(v1-v2));
    if (fabs(v1-v2)>1.0e-10)
      return(1);
  }
  printf("TESTING DFDZ VALUES...\n");
  for (i=0;i<8;i++) {
    tricubic_pointID2xyz(i,&x,&y,&z);
    v1=dfdz[i];
    v2=tricubic_eval(a,x,y,z,0,0,1);
    printf("%i\t%lf\t%lf\tError=%lg\n",i,v1,v2,fabs(v1-v2));
    if (fabs(v1-v2)>1.0e-10)
      return(1);
  }
  printf("TESTING D2FDXDY VALUES...\n");
  for (i=0;i<8;i++) {
    tricubic_pointID2xyz(i,&x,&y,&z);
    v1=d2fdxdy[i];
    v2=tricubic_eval(a,x,y,z,1,1,0);
    printf("%i\t%lf\t%lf\tError=%lg\n",i,v1,v2,fabs(v1-v2));
    if (fabs(v1-v2)>1.0e-10)
      return(1);
  }
  printf("TESTING D2FDXDZ VALUES...\n");
  for (i=0;i<8;i++) {
    tricubic_pointID2xyz(i,&x,&y,&z);
    v1=d2fdxdz[i];
    v2=tricubic_eval(a,x,y,z,1,0,1);
    printf("%i\t%lf\t%lf\tError=%lg\n",i,v1,v2,fabs(v1-v2));
    if (fabs(v1-v2)>1.0e-10)
      return(1);
  }
  printf("TESTING D2FDYDZ VALUES...\n");
  for (i=0;i<8;i++) {
    tricubic_pointID2xyz(i,&x,&y,&z);
    v1=d2fdydz[i];
    v2=tricubic_eval(a,x,y,z,0,1,1);
    printf("%i\t%lf\t%lf\tError=%lg\n",i,v1,v2,fabs(v1-v2));
    if (fabs(v1-v2)>1.0e-10)
      return(1);
  }
  printf("TESTING D3FDXDYDZ VALUES...\n");
  for (i=0;i<8;i++) {
    tricubic_pointID2xyz(i,&x,&y,&z);
    v1=d3fdxdydz[i];
    v2=tricubic_eval(a,x,y,z,1,1,1);
    printf("%i\t%lf\t%lf\tError=%lg\n",i,v1,v2,fabs(v1-v2));
    if (fabs(v1-v2)>1.0e-10)
      return(1);
  }

  return(0);

}

int test2(void) {
  int i,j;
  double v1,v2;
  double x,y,z;
  double f1[8],df1dx[8],df1dy[8],df1dz[8],d2f1dxdy[8],d2f1dxdz[8],d2f1dydz[8],d3f1dxdydz[8];
  double f2[8],df2dx[8],df2dy[8],df2dz[8],d2f2dxdy[8],d2f2dxdz[8],d2f2dydz[8],d3f2dxdydz[8];
  double a1[64],a2[64];
  //FILE *f;
  double rho;

  for (i=0;i<8;i++) {
    f1[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    df1dx[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    df1dy[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    df1dz[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    d2f1dxdy[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    d2f1dxdz[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    d2f1dydz[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    d3f1dxdydz[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    f2[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    df2dx[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    df2dy[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    df2dz[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    d2f2dxdy[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    d2f2dxdz[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    d2f2dydz[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
    d3f2dxdydz[i]=-MAX+2*MAX*((double)(rand()))/((double)(RAND_MAX));
  }

  /* f1 on the upper face mut be equal to f2 on the lower face */

  int iarr1[4]={4,5,6,7};
  int iarr2[4]={0,1,2,3};

  for (i=0;i<4;i++) {
    f1[iarr1[i]]=f2[iarr2[i]];
    df1dx[iarr1[i]]=df2dx[iarr2[i]];
    df1dy[iarr1[i]]=df2dy[iarr2[i]];
    df1dz[iarr1[i]]=df2dz[iarr2[i]];
    d2f1dxdy[iarr1[i]]=d2f2dxdy[iarr2[i]];
    d2f1dxdz[iarr1[i]]=d2f2dxdz[iarr2[i]];
    d2f1dydz[iarr1[i]]=d2f2dydz[iarr2[i]];
    d3f1dxdydz[iarr1[i]]=d3f2dxdydz[iarr2[i]];
  }

  tricubic_get_coeff(a1,f1,df1dx,df1dy,df1dz,d2f1dxdy,d2f1dxdz,d2f1dydz,d3f1dxdydz);
  tricubic_get_coeff(a2,f2,df2dx,df2dy,df2dz,d2f2dxdy,d2f2dxdz,d2f2dydz,d3f2dxdydz);

  printf("CONTINUITY CHECK...\n");
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on F:         %lg\n",sqrt(rho));
  if (sqrt(rho)>1.0e-8)
    return(1);
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,1,0,0);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,1,0,0);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on DFDX:      %lg\n",sqrt(rho));
  if (sqrt(rho)>1.0e-8)
    return(1);
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,0,1,0);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,0,1,0);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on DFDY:      %lg\n",sqrt(rho));
  if (sqrt(rho)>1.0e-8)
    return(1);
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,0,0,1);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,0,0,1);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on DFDZ:      %lg\n",sqrt(rho));
  if (sqrt(rho)>1.0e-8)
    return(1);
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,1,1,0);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,1,1,0);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D2FDXDY:   %lg\n",sqrt(rho));
  if (sqrt(rho)>1.0e-8)
    return(1);
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,1,0,1);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,1,0,1);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D2FDXDZ:   %lg\n",sqrt(rho));
  if (sqrt(rho)>1.0e-8)
    return(1);
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,0,1,1);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,0,1,1);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D2FDXDY:   %lg\n",sqrt(rho));
  if (sqrt(rho)>1.0e-8)
    return(1);
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,1,1,1);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,1,1,1);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D3FDXDYDZ: %lg\n",sqrt(rho));
  if (sqrt(rho)>1.0e-8)
    return(1);

  printf("THESE ARE NOT NECESSARILY CONTINUOUS...\n");

  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,2,0,0);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,2,0,0);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D2FDX2:    %lg\n",sqrt(rho));
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,0,2,0);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,0,2,0);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D2FDY2:    %lg\n",sqrt(rho));
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,0,0,2);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,0,0,2);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D2FDZ2:    %lg\n",sqrt(rho));
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,3,0,0);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,3,0,0);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D3FDX3:    %lg\n",sqrt(rho));
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,0,3,0);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,0,3,0);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D3FDY3:    %lg\n",sqrt(rho));
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,0,0,3);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,0,0,3);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D3FDZ3:    %lg\n",sqrt(rho));
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,2,1,0);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,2,1,0);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D3FDX2DY:  %lg\n",sqrt(rho));
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,2,0,1);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,2,0,1);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D3FDX2DZ:  %lg\n",sqrt(rho));
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,1,2,0);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,1,2,0);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D3FDXDY2:  %lg\n",sqrt(rho));
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,1,0,2);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,1,0,2);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D3FDXDZ2:  %lg\n",sqrt(rho));
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,0,2,1);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,0,2,1);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D3FDY2DZ:  %lg\n",sqrt(rho));
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=tricubic_eval(a1,x,y,z,0,1,2);
      z=(double)(0.0);
      v2=tricubic_eval(a2,x,y,z,0,1,2);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D3FDYDZ2:  %lg\n",sqrt(rho));

  return(0);

}


int main(int narg, char *arg[]) {

  printf("LIBTRICUBIC TEST PROGRAM\n");

  printf("************************\n");

  printf("libtricubic v%s\n",tricubic_version());

  printf("************************\n");

  test1();

  printf("************************\n");

  test2();

  printf("************************\n");

  printf("LIBTRICUBIC END OF TESTS\n");

  return(0);
}
