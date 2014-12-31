/************************************************************/
/* example1.cpp : illustrates the use of libtricubic        */
/* Francois Lekien <lekien@mit.edu> 2004-01-20              */
/************************************************************/

///Required Include File: tricubic.h
///If tricubic.h has not been installed in a directory
///  accessible by the compiler, use -I/path/to/tricubic
#include <tricubic.h>

///Define the box
///Tricubic is written for cubes of side 1. Multiplications and
///  divisions are needed along the way for rectangular box with
///  arbitrary sides. See below for details.
#define dx 0.2
#define dy 2.0
#define dz 0.3

///These are the 8 functions that need to be known at the 8 corners.
///The functions are f, the three first derivatives dfdx, dfdy, dfdz,
///  the 3 mixed 2nd order derivatives and the mixed 3rd order derivative.
///The derivatives can be obtained by numercal differentiation of f at the
///  other corners of the grid. See Numerical Recipies for examples in 2D
///The order of the points is as follow:
///  0: x=0; y=0; z=0;
///  1: x=1; y=0; z=0;
///  2: x=0; y=1; z=0;
///  3: x=1; y=1; z=0;
///  4: x=0; y=0; z=1;
///  5: x=1; y=0; z=1;
///  6: x=0; y=1; z=1;
///  7: x=1; y=1; z=1;
///For convenience, the ordering of the points is available at run time
///   using tricubic_pointID2xyz(). See man page for details
double fval[8]={1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9};
double dfdxval[8]={1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9};
double dfdyval[8]={1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9};
double dfdzval[8]={1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9};
double d2fdxdyval[8]={1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9};
double d2fdxdzval[8]={1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9};
double d2fdydzval[8]={1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9};
double d3fdxdydzval[8]={1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9};

int main(int narg, char *arg[]) {
  double a[64];
  double f1, f2, dfdx, d2fdxdy;
  int i;
  double x,y,z;
  ///First we compute the 64 coefficients for the cube.
  ///The first step is to scale the derivatives that have been computed
  ///  in a rectangular box instead of the cube of side 1.
  ///Notice that this step can be avoided by computing the numerical
  ///  derivatives without dividing by the length of the boxes.
  for (i=0;i<8;i++) {
    fval[i]*=1.0;
    dfdxval[i]*=dx;
    dfdyval[i]*=dy;
    dfdzval[i]*=dz;
    d2fdxdyval[i]*=dx*dy;
    d2fdxdzval[i]*=dx*dz;
    d2fdydzval[i]*=dy*dz;
    d3fdxdydzval[i]*=dx*dy*dz;
  } 
  ///Next we get the set of coefficient for this cube
  tricubic_get_coeff(a,fval,dfdxval,dfdyval,dfdzval,d2fdxdyval,d2fdxdzval,d2fdydzval,d3fdxdydzval);
  ///To get the value in the middle of the cube, we always use (.5,.5,.5)
  ///  i.e. relative coordinates
  f1=tricubic_eval(a,.5,.5,.5);
  ///To get the value at a point x,y,z (with referrence to corner ID 0
  ///  we devide by each length
  f2=tricubic_eval(a,x/dx,y/dy,z/dz);
  ///Derivatives can be computed similarly but need to be scaled by dx,dy,dz
  dfdx=tricubic_eval(a,x/dx,y/dy,z/dz,1,0,0)/dx;
  d2fdxdy=tricubic_eval(a,x/dx,y/dy,z/dz,1,1,0)/(dx*dy);
}
