#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include"libs/matfun.h"
int main() {
    double theta1 = M_PI;
    double theta2 = M_PI / 3;
    
    double x1 = cos(theta1);
    double y1 = sin(theta1);
    double x2 = cos(theta2);
    double y2 = sin(theta2);
   
    int m=2,n=1;
    
    double **a = createMat(m,n); 
    	a[0][0]=x1;
	a[1][0]=y1;

    double **b = createMat(m,n); 
    	b[0][0]=x2;
	b[1][0]=y2;

    double **c = Matadd(a,b,m,n);
    double norm = Matnorm(c,m);
    //printing the matrices
    printMat(a,m,n); //printing Mat a
    printMat(b,m,n); //printing Mat b
    printMat(c,m,n); //printing Mat c i.e (a+b)
    printf("||a+b||= %lf\n",norm);
    return 0;
}

