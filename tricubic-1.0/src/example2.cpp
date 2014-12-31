#include <stdio.h>
#include <stdlib.h>
#include <tricubic.h>
#include <math.h>

#define MAX 10.0
#define Nt 101

double monoA[4][4]={{1,0,0,0},{0,0,1,0},{-3,3,-2,-1},{2,-2,1,1}};

double monocubic_run(double a[4], double x) {
  double ret=(double)(0.0);
  int i;
  for (i=0;i<4;i++) {
    if (i!=0) {
      ret+=a[i]*pow(x,(double)(i));
    } else {
      ret+=a[i];
    }
  }
  return(ret);
}

void monocubic_mat(double a[4], double x[4]) {
  int i,j;
  for (i=0;i<4;i++) {
    a[i]=(double)(0.0);
    for (j=0;j<4;j++) {
      a[i]+=monoA[i][j]*x[j];
    }
  }
}

double monocubic(double f[8],double dfdx[8],double dfdy[8],double dfdz[8],double d2fdxdy[8], double d2fdxdz[8], double d2fdydz[8], double d3fdxdydz[8], double x, double y, double z, int dx, int dy, int dz) {
  double fx[4];
  double dfdyx[4];
  double dfdzx[4];
  double fy[2];
  double dfdzy[2];
  double d2fdydzx[2];
  int i;
  double coeff[4];
  double val[4];

  double dq = 1.0e-5;
  if (dx>0) {
    return((monocubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x+dq,y,z,dx-1,dy,dz)-monocubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x-dq,y,z,dx-1,dy,dz))/(2*dq));
  }
  if (dy>0) {
    return((monocubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y+dq,z,dx,dy-1,dz)-monocubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y-dq,z,dx,dy-1,dz))/(2*dq));
  }
  if (dz>0) {
    return((monocubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z+dq,dx,dy,dz-1)-monocubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z-dq,dx,dy,dz-1))/(2*dq));
  }

  //first we interpolate 4 values wrt x
  for (i=0;i<4;i++) {
    val[0]=f[2*i];
    val[1]=f[2*i+1];
    val[2]=dfdx[2*i];
    val[3]=dfdx[2*i+1];
    monocubic_mat(coeff,val);
    fx[i]=monocubic_run(coeff,x);
  }
  //next we interpolate 4 values dfdy wrt x
  for (i=0;i<4;i++) {
    val[0]=dfdy[2*i];
    val[1]=dfdy[2*i+1];
    val[2]=d2fdxdy[2*i];
    val[3]=d2fdxdy[2*i+1];
    monocubic_mat(coeff,val);
    dfdyx[i]=monocubic_run(coeff,x);
  }
  //next we interpolate 4 values dfdz wrt x
  for (i=0;i<4;i++) {
    val[0]=dfdz[2*i];
    val[1]=dfdz[2*i+1];
    val[2]=d2fdxdz[2*i];
    val[3]=d2fdxdz[2*i+1];
    monocubic_mat(coeff,val);
    dfdzx[i]=monocubic_run(coeff,x);
  }
  //next we interpolate 4 values d2fdydz wrt x
  for (i=0;i<4;i++) {
    val[0]=d2fdydz[2*i];
    val[1]=d2fdydz[2*i+1];
    val[2]=d3fdxdydz[2*i];
    val[3]=d3fdxdydz[2*i+1];
    monocubic_mat(coeff,val);
    d2fdydzx[i]=monocubic_run(coeff,x);
  }
  //next we interpolate 2 values of f wrt y
  for (i=0;i<2;i++) {
    val[0]=fx[2*i];
    val[1]=fx[2*i+1];
    val[2]=dfdyx[2*i];
    val[3]=dfdyx[2*i+1];
    monocubic_mat(coeff,val);
    fy[i]=monocubic_run(coeff,y);
  }
  //next we interpolate 2 values of dfdz wrt y
  for (i=0;i<2;i++) {
    val[0]=dfdzx[2*i];
    val[1]=dfdzx[2*i+1];
    val[2]=d2fdydzx[2*i];
    val[3]=d2fdydzx[2*i+1];
    monocubic_mat(coeff,val);
    dfdzy[i]=monocubic_run(coeff,y);
  }
  //finally interpolation of f wrt z
  val[0]=fy[0];
  val[1]=fy[1];
  val[2]=dfdzy[0];
  val[3]=dfdzy[1];
  monocubic_mat(coeff,val);
  return(monocubic_run(coeff,z));

  return(0.0);
}

int test1(void) {
  int i;
  double v1,v2;
  double x,y,z;
  double f[8],dfdx[8],dfdy[8],dfdz[8];
  double d2fdxdy[8], d2fdxdz[8];
  double d2fdydz[8], d3fdxdydz[8];
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

  printf("TESTING F VALUES...\n");
  for (i=0;i<8;i++) {
    tricubic_pointID2xyz(i,&x,&y,&z);
    v1=f[i];
    v2=monocubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z,0,0,0);
    printf("%i\t%lf\t%lf\tError=%lg\n",i,v1,v2,fabs(v1-v2));
  }

  printf("TESTING DFDX VALUES...\n");
  for (i=0;i<8;i++) {
    tricubic_pointID2xyz(i,&x,&y,&z);
    v1=dfdx[i];
    v2=monocubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z,1,0,0);
    printf("%i\t%lf\t%lf\tError=%lg\n",i,v1,v2,fabs(v1-v2));
  }

  printf("TESTING DFDY VALUES...\n");
  for (i=0;i<8;i++) {
    tricubic_pointID2xyz(i,&x,&y,&z);
    v1=dfdy[i];
    v2=monocubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z,0,1,0);
    printf("%i\t%lf\t%lf\tError=%lg\n",i,v1,v2,fabs(v1-v2));
  }

  printf("TESTING DFDZ VALUES...\n");
  for (i=0;i<8;i++) {
    tricubic_pointID2xyz(i,&x,&y,&z);
    v1=dfdz[i];
    v2=monocubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z,0,0,1);
    printf("%i\t%lf\t%lf\tError=%lg\n",i,v1,v2,fabs(v1-v2));
  }

  printf("TESTING D2FDXDY VALUES...\n");
  for (i=0;i<8;i++) {
    tricubic_pointID2xyz(i,&x,&y,&z);
    v1=d2fdxdy[i];
    v2=monocubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z,1,1,0);
    printf("%i\t%lf\t%lf\tError=%lg\n",i,v1,v2,fabs(v1-v2));
  }

  printf("TESTING D2FDXDZ VALUES...\n");
  for (i=0;i<8;i++) {
    tricubic_pointID2xyz(i,&x,&y,&z);
    v1=d2fdxdz[i];
    v2=monocubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z,1,0,1);
    printf("%i\t%lf\t%lf\tError=%lg\n",i,v1,v2,fabs(v1-v2));
  }

  printf("TESTING D2FDYDZ VALUES...\n");
  for (i=0;i<8;i++) {
    tricubic_pointID2xyz(i,&x,&y,&z);
    v1=d2fdydz[i];
    v2=monocubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z,0,1,1);
    printf("%i\t%lf\t%lf\tError=%lg\n",i,v1,v2,fabs(v1-v2));
  }

  printf("TESTING D3FDXDYDZ VALUES...\n");
  for (i=0;i<8;i++) {
    tricubic_pointID2xyz(i,&x,&y,&z);
    v1=d3fdxdydz[i];
    v2=monocubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z,1,1,1);
    printf("%i\t%lf\t%lf\tError=%lg\n",i,v1,v2,fabs(v1-v2));
  }

  return(0);
}

int test2(void) {
  int i,j;
  double v1,v2;
  double x,y,z;
  double f1[8],df1dx[8],df1dy[8],df1dz[8],d2f1dxdy[8],d2f1dxdz[8],d2f1dydz[8],d3f1dxdydz[8];
  double f2[8],df2dx[8],df2dy[8],df2dz[8],d2f2dxdy[8],d2f2dxdz[8],d2f2dydz[8],d3f2dxdydz[8];
  FILE *f;
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

  printf("CONTINUITY CHECK...\n");
  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=monocubic(f1,df1dx,df1dy,df1dz,d2f1dxdy,d2f1dxdz,d2f1dydz,d3f1dxdydz,x,y,z,0,0,0);
      z=(double)(0.0);
      v2=monocubic(f2,df2dx,df2dy,df2dz,d2f2dxdy,d2f2dxdz,d2f2dydz,d3f2dxdydz,x,y,z,0,0,0);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on F:         %lg\n",sqrt(rho));

  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=monocubic(f1,df1dx,df1dy,df1dz,d2f1dxdy,d2f1dxdz,d2f1dydz,d3f1dxdydz,x,y,z,1,0,0);
      z=(double)(0.0);
      v2=monocubic(f2,df2dx,df2dy,df2dz,d2f2dxdy,d2f2dxdz,d2f2dydz,d3f2dxdydz,x,y,z,1,0,0);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on DFDX:      %lg\n",sqrt(rho));

  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=monocubic(f1,df1dx,df1dy,df1dz,d2f1dxdy,d2f1dxdz,d2f1dydz,d3f1dxdydz,x,y,z,0,1,0);
      z=(double)(0.0);
      v2=monocubic(f2,df2dx,df2dy,df2dz,d2f2dxdy,d2f2dxdz,d2f2dydz,d3f2dxdydz,x,y,z,0,1,0);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on DFDY:      %lg\n",sqrt(rho));

  rho=(double)(0.0);
  f=fopen("testmonoc.dat","w");
  fprintf(f,"VARIABLE=\"x\"\"y\"\"v1\"\"v2\"\"d\"\n");
  fprintf(f,"ZONE I=%i J=%i\n",Nt,Nt);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=monocubic(f1,df1dx,df1dy,df1dz,d2f1dxdy,d2f1dxdz,d2f1dydz,d3f1dxdydz,x,y,z,0,0,1);
      z=(double)(0.0);
      v2=monocubic(f2,df2dx,df2dy,df2dz,d2f2dxdy,d2f2dxdz,d2f2dydz,d3f2dxdydz,x,y,z,0,0,1);
      fprintf(f,"%g\t%g\t%g\t%g\t%g\n",x,y,v1,v2,v1-v2);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  fclose(f);
  printf("  C1-error on DFDZ:      %lg\n",sqrt(rho));

  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=monocubic(f1,df1dx,df1dy,df1dz,d2f1dxdy,d2f1dxdz,d2f1dydz,d3f1dxdydz,x,y,z,1,1,0);
      z=(double)(0.0);
      v2=monocubic(f2,df2dx,df2dy,df2dz,d2f2dxdy,d2f2dxdz,d2f2dydz,d3f2dxdydz,x,y,z,1,1,0);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D2FDXDY:   %lg\n",sqrt(rho));

  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=monocubic(f1,df1dx,df1dy,df1dz,d2f1dxdy,d2f1dxdz,d2f1dydz,d3f1dxdydz,x,y,z,1,0,1);
      z=(double)(0.0);
      v2=monocubic(f2,df2dx,df2dy,df2dz,d2f2dxdy,d2f2dxdz,d2f2dydz,d3f2dxdydz,x,y,z,1,0,1);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D2FDXDZ:   %lg\n",sqrt(rho));

  rho=(double)(0.0);
  for (j=0;j<Nt;j++) {
    for (i=0;i<Nt;i++) {
      x=0.0+1.0*i/((double)(Nt-1));
      y=0.0+1.0*j/((double)(Nt-1));
      z=(double)(1.0);
      v1=monocubic(f1,df1dx,df1dy,df1dz,d2f1dxdy,d2f1dxdz,d2f1dydz,d3f1dxdydz,x,y,z,0,1,1);
      z=(double)(0.0);
      v2=monocubic(f2,df2dx,df2dy,df2dz,d2f2dxdy,d2f2dxdz,d2f2dydz,d3f2dxdydz,x,y,z,0,1,1);
      rho=(rho<(v1-v2)*(v1-v2))?(v1-v2)*(v1-v2):rho;
    }
  }
  printf("  C1-error on D2FDYDZ:   %lg\n",sqrt(rho));


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
