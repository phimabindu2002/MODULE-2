#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"matfun.h"
//#include"coeffs-mat.h"
//#include"matfun.h"
int main() {
    int m = 2;  // You can set the dimensions of your matrix
    int n = 3;

// Create a matrix
	double **G_v = createMat(m, n);
	double **C_m= createMat(3,3);
	double **C_mid= createMat(3,3);
	double **R_o=createMat(2,2);
	double **C_mid_dir= createMat(3,3);
	double **C_alt= createMat(3,3);
	double **C_in= createMat(3,3);
	double **cont_mat= createMat(3,3);
	double **i_con= createMat(3,3);

// Reading matrix data from .dat file
G_v=loadMat("vert.dat",m,n);
C_m=loadMat("C.dat",3,3);
C_mid=loadMat("C_mid.dat",3,3);
C_mid_dir=loadMat("C_mid_dir.dat",3,3);
R_o=loadMat("R.dat",2,2);
C_alt=loadMat("C_alt.dat",3,3);
C_in=loadMat("C_in.dat",3,3);



	// **********************VECTORS*****************************
	double **G_dir=Matmul(G_v,C_m,2,3,3);
	double **G_n=Matmul(R_o,G_dir,2,2,3);
	double **G_con=diag(Matmul(transposeMat(G_n,2,3),G_v,3,2,3),3,3);
	double **G_dis=sqrt_diag(Matmul(transposeMat(G_dir,2,3),G_dir,3,2,3),3,3);
	double **G_line=Mathstack(transposeMat(G_n,2,3),transposeMat(G_con,1,3),3,2,1);
	//***********************MEDIANS*******************************
	double **G_mid=Matscale(Matmul(G_v,C_mid,2,3,3),2,3,0.5);
	double **G_med_dir=Matmul(G_v,C_mid_dir,2,3,3);
	double **G_n_med = Matmul(R_o,G_med_dir,2,2,3);
	double **cmat_med=diag(Matmul(transposeMat(G_n_med,2,3),G_v,3,2,3),3,3);
	double **linemat_med=Mathstack(transposeMat(G_n_med,2,3),transposeMat(cmat_med,1,3),3,2,1);	
	double **G_G=line_intersect(linemat_med,3,3);
	//***********************ALTITUDE*******************************
	double **G_n_alt=Matmul(G_v,C_alt,2,3,3);
	double **cmat_alt=diag(Matmul(transposeMat(G_n_alt,2,3),G_v,3,2,3),3,3);
	double **linemat_alt=Mathstack(transposeMat(G_n_alt,2,3),transposeMat(cmat_alt,1,3),3,2,1);	
	double **G_H=line_intersect(linemat_alt,3,3);
	//******************PERPENDICULAR BISECTOR**********************
	double **cmat_perp_bis=diag(Matmul(transposeMat(G_n_alt,2,3),G_mid,3,2,3),3,3);
	double **linemat_perp_bis=Mathstack(transposeMat(G_n_alt,2,3),transposeMat(cmat_perp_bis,1,3),3,2,1);	 
	double **G_O=line_intersect(linemat_perp_bis,3,3);
	//********************ANGULAR  BISECTOR************************
        double **secvec=Matscale(Matmul(G_dis,C_in, 1,3,3),1,3,0.5);
        i_con[0][0]=  -1 ;      i_con[0][1]=secvec[0][2]/G_dis[0][2];     i_con[0][2]= secvec[0][1]/G_dis[0][0] ;
        i_con[1][0]=secvec[0][2]/G_dis[0][1] ;       i_con[1][1]= -1 ;       i_con[1][2]=secvec[0][0]/G_dis[0][0];
        i_con[2][0]= secvec[0][1]/G_dis[0][1];      i_con[2][1]= secvec[0][0]/G_dis[0][2] ;    i_con[2][2]= -1 ;
        double **a_dir=Matmul(G_v,i_con,2,3,3);
        double **a_nor=Matmul(R_o,a_dir,2,2,3);
        double **a_cof=diag(Matmul(transposeMat(a_nor,2,3),G_v,3,2,3),3,3);
        double **a_line=Mathstack(transposeMat(a_nor,2,3),transposeMat(a_cof,1,3),3,2,1);
        double **a_i = line_intersect(a_line,3,3);


        cont_mat[0][0]=  0 ;    cont_mat[0][1]=secvec[0][2]/G_dis[0][2];     cont_mat[0][2]= secvec[0][1]/G_dis[0][0] ;
        cont_mat[1][0]=secvec[0][2]/G_dis[0][1] ;       cont_mat[1][1]= 0 ;      cont_mat[1][2]=secvec[0][0]/G_dis[0][0];
        cont_mat[2][0]= secvec[0][1]/G_dis[0][1];      cont_mat[2][1]= secvec[0][0]/G_dis[0][2] ;    cont_mat[2][2]= 0 ;
        double **G_i=Matmul(G_v,cont_mat,2,3,3);
        double **G_imid=Matscale(Matmul(G_i,C_mid,2,3,3),2,3,0.5);
        double **G_idir_alt=Matmul(G_i,C_alt,2,3,3);
        double **cmat_iperp_bis=diag(Matmul(transposeMat(G_idir_alt,2,3),G_imid,3,2,3),3,3);
	double **linemat_imed=Mathstack(transposeMat(G_idir_alt,2,3),transposeMat(cmat_iperp_bis,1,3),3,2,1);
        double **G_I=line_intersect(linemat_imed,3,3);
	//*******************Eigen vector approach to find the contact points******************************
         double **h=createMat(2,1);  
        h[0][0]=G_v[0][0]-G_I[0][0];    h[1][0]=G_v[1][0]-G_I[0][1];
        double **V=createMat(2,2);
        V[0][0]=1;      V[0][1]=0;      V[1][0]=0;       V[1][1]=1;
        double **u=createMat(2,1);
        u[0][0]=G_I[0][0]-G_I[0][0];    u[1][0]=G_I[0][1]-G_I[0][1];
        double **f=createMat(1,1);
        f[0][0]=sqrt(pow(G_I[0][0]-G_i[0][0],2)+pow(G_I[0][1]-G_i[1][0],2));
        f[0][0]=-f[0][0]*f[0][0];
double **gh=Matadd(Matadd(Matmul(Matmul(transposeMat(h,2,1),V,1,2,2),h,1,2,1), Matscale( Matmul(transposeMat(u,2,1),h,1,2,1),1,1,2) ,1,1),f,1,1);
double **sigmat=Matsub(Matmul(Matadd(Matmul(V,h,2,2,1),u,2,1),transposeMat(Matadd(Matmul(V,h,2,2,1),u,2,1),2,1),2,1,2) ,Matscale(V,2,2,gh[0][0]) ,2,2);





double **E_val=Mateigval(sigmat);
double **P=Mateigvec(sigmat);
double **u1=createMat(2,1);
u1[0][0]=sqrt(fabs(E_val[1][0]));  u1[1][0]=sqrt(fabs(E_val[0][0]));
double **u2=createMat(2,1);
u2[0][0]=sqrt(fabs(E_val[1][0]));  u2[1][0]=-sqrt(fabs(E_val[0][0]));


double **m1=Matmul(P,u1,2,2,1);
double **m2=Matmul(P,u2,2,2,1);

double **mu1n=Matmul(transposeMat(m1,2,1),Matadd(Matmul(V,h,2,2,1),u,2,1),1,2,1);  
double **mu1d=Matmul(transposeMat(m1,2,1),Matmul(V,m1,2,2,1),1,2,1);
double mu1=-mu1n[0][0]/mu1d[0][0];

double **mu2n=Matmul(transposeMat(m2,2,1),Matadd(Matmul(V,h,2,2,1),u,2,1),1,2,1);  
double **mu2d=Matmul(transposeMat(m2,2,1),Matmul(V,m2,2,2,1),1,2,1);
double mu2=-mu2n[0][0]/mu2d[0][0];

double **t1=createMat(2,1); double **t2=createMat(2,1);
t1[0][0]=mu1*m1[0][0]; t1[1][0]=mu1*m1[1][0];
t2[0][0]=mu2*m2[0][0]; t2[1][0]=mu2*m2[1][0];

double **E=Matadd(h,t1,2,1);
double **F=Matadd(h,t2,2,1);
	E[0][0]=E[0][0]+G_I[0][0];	E[1][0]=E[1][0]+G_I[0][1];
	F[0][0]=F[0][0]+G_I[0][0];	F[1][0]=F[1][0]+G_I[0][1];


// Printing matrices . 	
printf("\n Vectors \n");
printf("vertices matrix= \n");
printMat(G_v,2,3);
//printf("direction matrix= \n");
//printMat(G_dir,2,3);
//printf("normal matrix= \n");
//printMat(G_n,2,3);
//printf("constant matrix= \n");
//printMat(G_con,1,3);
//printf("distance matrix= \n");
//printMat(G_dis,1,3);
//printf("line  matrix= \n");
//printMat(G_line,3,3);
//
//printf("\n Medians \n");
//printf("midpoint matrix= \n");
//printMat(G_mid,2,3);
//printf("median direction  matrix= \n");
//printMat(G_med_dir,2,3);
//printf("median normal matrix= \n");
//printMat(G_n_med,2,3);
//printf("median constant matrix= \n");
//printMat(cmat_med,1,3);
//printf("median line matrix= \n");
//printMat(linemat_med,3,3);
//printf("Centroid = \n");
//printMat(G_G,1,2);
//
//printf("\n Altitude \n");
//printf("altitude normal matrix= \n");
//printMat(G_n_alt,2,3);
//printf("altitude constant matrix= \n");
//printMat(cmat_alt,1,3);
//printf("Altitude line matrix= \n");
//printMat(linemat_alt,3,3);
//printf("Orthocentre = \n");
//printMat(G_H,1,2);
//
//printf("\n Perpendicular Bisector \n");
//printf("perp_bisect constant  matrix= \n");
//printMat(cmat_perp_bis,1,3);
//printf("perp_bis line matrix= \n");
//printMat(linemat_perp_bis,3,3);
//printf("Circumcentre = \n");
//printMat(G_O,1,2);
//
//printf("\n Angular Bisector \n");
//printf("m,n,p values = \n");
//printMat(secvec,1,3);
//printf("ang_bis direction matrix= \n");
//printMat(a_dir,2,3);
//printf("ang_bis normal matrix= \n");
//printMat(a_nor,2,3);
//printf("ang_bis constant matrix= \n");
//printMat(a_cof,1,3);
//printf("ang_bis line  matrix= \n");
//printMat(a_line,3,3);
//printf("ang_bis intersection points= \n");
//printMat(a_i,1,2);
printf("Eigen vector Approach\n");
printf("Incentre = \n");
printMat(G_I,1,2);
printf("contact points = \n");
printMat(G_i,2,3);
printf(" h = \n");
printMat(h,2,1);
printf("V = \n");
printMat(V,2,2);
printf("u = \n");
printMat(u,2,1);
printf("iradius = \n");
printMat(f,1,1);
printf("gh = \n");
printMat(gh,1,1);
printf("sigmat = \n");
printMat(sigmat,2,2);
printf("Eigen values = \n");
printMat(E_val,2,1);
printf("Eigen vectors = \n");
printMat(P,2,2);
//printf("gauss mat = \n");
//printMat(ga,2,2);
printf("u1 = \n");
printMat(u1,2,1);
printf("u2 = \n");
printMat(u2,2,1);
printf("m1 = \n");
printMat(m1,2,1);
printf("m2 = \n");
printMat(m2,2,1);
printf("E = \n");
printMat(E,2,1);
printf("F = \n");
printMat(F,2,1);
printf("%lf/n",sqrt(abs(E_val[1][0])));
}
