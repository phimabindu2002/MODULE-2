//Code by G V V Sharma
//October 27, 2023
//matrix operations using arrays
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "him.h"
//printing the matrix

int  main()
{

FILE *fp; //file pointer
double **A,**B,**C,**G_v,**D,**E,**F;//declare vertices
double **m1,**m2,**m3,**m4,**m5;//direction vectors
double **n1,**n2,**n3,**n4,**n11,**n12,**n13;//normal vectors
double **temp;//temporary array
double c1,c2,c3;//direction vectors
double a,b,c; //triangle sides
int m =2, n=3, i, j;
double **mat =createMat(m,n);//vertices matrix 
mat = loadMat("vertices.dat",m, n);//loading matrix from file
//temp= loadMat("circ.dat",4, 1);//loading matrix from file
temp= loadMat("matex.dat",2, 2);//loading matrix from file

//Extracting vertices
A = Matcol(mat,m,0);
B = Matcol(mat,m,1);
C = Matcol(mat,m,2);

//printf("point matrix:\n");
//printMat(mat,2,3);

//*****************************************VECTORS*****************************************
printf("\nVECTORS:       \n");
printf("\ndirection matrix:\n");

//Direction vectors
m1 = Matsub(A,B,m,1); //m1=A-B
m2 = Matsub(B,C,m,1);  //m2=B-C
m3 = Matsub(C,A,m,1); //m3=C-A
m4 = p_m(m1,m2,m3,2,3);
printMat(m4,2,3);
printf("\n");

//Normal vectors
n1 = normVec(m1);
n2 = normVec(m2);
n3 = normVec(m3);
printf("normal matrix:\n");
n4 = p_m(n1,n2,n3,2,3);
printMat(n4,2,3);


c1 = Matdot(n1,A,2);
c2 = Matdot(n2,B,2);
c3 = Matdot(n3,C,2);

//Sides
a = Matnorm(m1,m);
b = Matnorm(m2,m);
c = Matnorm(m3,m);

printf("\nconstant matrix :\n %lf , %lf , %lf \n",c1,c2,c3);//AB line coefficient

//Triangle sides
//printf("\nDistance:\n %lf , %lf ,%lf \n",a,b,c);
double **G_dis=matrix(a,b,c);
printf("\nDistane matrix\n");
printMat(G_dis,1,3);
//***************************************END OF VECTORS*************************************

//******************************************MEDIAN******************************************
printf("\n\nMEDIANS \n");

//mid points
double **M_P=Mid_point(A,B,C);
printf("\nMID POINTS\n");
printMat(M_P,2,3);

D=Matcol(M_P,2,0);
E=Matcol(M_P,2,1);
F=Matcol(M_P,2,2);

//mid direction matrix
double **M_P_D = p_m(Matsub(A,D,2,1),Matsub(B,E,2,1),Matsub(C,F,2,1),2,3);
printf("\nMID POINT DIRECTION MATRIX\n");
printMat(M_P_D,2,3);
double **M_P_D_1 , **M_P_D_2 , **M_P_D_3 ;
M_P_D_1=Matcol(M_P_D,m,0);
M_P_D_2=Matcol(M_P_D,m,1);
M_P_D_3=Matcol(M_P_D,m,2);

//normal matrix
n1=normVec(M_P_D_1);
n2=normVec(M_P_D_2);
n3=normVec(M_P_D_3);
n4 = p_m(n1,n2,n3,2,3);
printf("\nmedian normal matrix\n");
printMat(n4,2,3);

//median constant matrix
c1 = Matdot(n1,D,2);
c2 = Matdot(n2,E,2);
c3 = Matdot(n3,F,2);

double **const_mat=matrix(c1,c2,c3);
printf("\nmedian Constant matrix\n");
printMat(const_mat,1,3);

//centroid
n11=Matcol(transposeMat(n4,2,3),m,0);
n12=Matcol(transposeMat(n4,2,3),m,1);

double **linemat_med=h_concat(transposeMat(n4,2,3),transposeMat(const_mat,1,3),3,2,3,1);	
double **G_G=line_intersect(linemat_med,3,3);
//printMat(linemat_med,3,3);
printf("\ncentroid:\n");
printMat(G_G,1,2);
//*************************************END OF MEDIAN***************************************


//***************************************ALTITUDE******************************************

printf("\n\nALTITUDE\n");
printf("\nAltitude normal matrix\n");
m5 = p_m(m2,m3,m1,2,3);
printMat(m5,2,3);
printf("\nAltitude constant matrix\n");

c1 = Matdot(m2,A,2);
c2 = Matdot(m3,B,2);
c3 = Matdot(m1,C,2);
double **A_c=matrix(c1,c2,c3);
printMat(A_c,1,3);
printf("\nORTHOCENTRE\n");
double **linemat_alt=h_concat(transposeMat(m5,2,3),transposeMat(A_c,1,3),3,2,3,1);
double **G_H=line_intersect(linemat_alt,3,3);
printMat(G_H,1,2);
//***************************************END OF ALTITUDE************************************

//***************************************PERPENDICULAR BISECTOR*****************************
printf("\n\nPERPENDICULAR BISECTOR\n");
c1 = Matdot(m2,D,2);
c2 = Matdot(m3,E,2);
c3 = Matdot(m1,F,2);

double **P_D_C=matrix(c1,c2,c3);
printf("\nperpendicular Constant matrix\n");
printMat(P_D_C,1,3);
double **linemat_perp_bis=h_concat(transposeMat(m5,2,3),transposeMat(P_D_,1,3),3,2,3,1);	printf("\nPerpendicular line matrix\n");
printMat(linemat_perp_bis,3,3);
double **G_O=line_intersect(linemat_perp_bis,3,3);
printf("\nCIRCUMCENTRE:\n");
printMat(G_O,1,2);

//********************************END OF PERPENDICULAR BISECTOR*****************************
//************************************ANGULAR BISECTOR**************************************
double **C_in= createMat(3,3);
C_in=loadMat("C_in.dat",3,3);

//M,N,P values
printf("\n\nANGULAR BISECTOR\n");
printf("m,n,p\n");
double **secvec=Matscale(Matmul(G_dis,C_in, 1,3,3),1,3,0.5);
printMat(secvec,1,3);

//angular direction matrix
printf("\nANGULAR BISECTOR Direction matrix\n");
double **i_con= createMat(3,3);
i_con[0][0]=  -1 ;      
i_con[0][1]=secvec[0][2]/G_dis[0][2];    
i_con[0][2]= secvec[0][1]/G_dis[0][0] ;
i_con[1][0]=secvec[0][2]/G_dis[0][1] ;
i_con[1][1]= -1 ;
i_con[1][2]=secvec[0][0]/G_dis[0][0];
i_con[2][0]= secvec[0][1]/G_dis[0][1];
i_con[2][1]= secvec[0][0]/G_dis[0][2] ;
i_con[2][2]= -1 ;

double **a_dir=Matmul(mat,i_con,2,3,3);
printMat(a_dir,2,3);

//Normal matrix
printf("\nAngular normal matrix\n");
double **a_n1,**a_n2,**a_n3,**a_nor;
a_n1=normVec(Matcol(a_dir,2,0));
a_n2=normVec(Matcol(a_dir,2,1));
a_n3=normVec(Matcol(a_dir,2,2));
a_nor=p_m(a_n1,a_n2,a_n3,2,3);
printMat(a_nor,2,3);

//constant matrix
printf("\nAngular constant matrix= \n");
double a_c1,a_c2,a_c3,**a_c_m;
a_c1 = Matdot(a_n1,A,2);
a_c2 = Matdot(a_n2,B,2);
a_c3 = Matdot(a_n3,C,2);
a_c_m=matrix(a_c1,a_c2,a_c3);
printMat(a_c_m,1,3);

//Angular intersection points
double **a_line=h_concat(transposeMat(a_nor,2,3),transposeMat(a_c_m,1,3),3,2,3,1);
double **a_i = line_intersect(a_line,3,3);
printf("\nAngular Intersection points\n");
printMat(a_i,1,2);


//Incentre
double **C_mid= createMat(3,3);
C_mid=loadMat("C_mid.dat",3,3);
double **C_alt= createMat(3,3);
C_alt=loadMat("C_alt.dat",3,3);
double **cont_mat= createMat(3,3);

cont_mat[0][0]=  0 ;
cont_mat[0][1]=secvec[0][2]/G_dis[0][2];
cont_mat[0][2]= secvec[0][1]/G_dis[0][0] ;
cont_mat[1][0]=secvec[0][2]/G_dis[0][1] ;
cont_mat[1][1]= 0 ;
cont_mat[1][2]=secvec[0][0]/G_dis[0][0];
cont_mat[2][0]= secvec[0][1]/G_dis[0][1];
cont_mat[2][1]= secvec[0][0]/G_dis[0][2] ;
cont_mat[2][2]= 0 ;

double **G_i=Matmul(mat,cont_mat,2,3,3);
double **G_imid=Mid_point(Matcol(G_i,2,0),Matcol(G_i,2,1),Matcol(G_i,2,2));
double **G_idir_alt=Matmul(G_i,C_alt,2,3,3);
double a_cc1,a_cc2,a_cc3,**cmat_iperp_bis;
a_cc1 = Matdot(Matcol(G_idir_alt,2,0),Matcol(G_imid,2,0),2);
a_cc2 = Matdot(Matcol(G_idir_alt,2,1),Matcol(G_imid,2,1),2);
a_cc3 = Matdot(Matcol(G_idir_alt,2,2),Matcol(G_imid,2,2),2);

cmat_iperp_bis=matrix(a_cc1,a_cc2,a_cc3);
double **linemat_imed=h_concat(transposeMat(G_idir_alt,2,3),transposeMat(cmat_iperp_bis,1,3),3,2,3,1);
double **G_I=line_intersect(linemat_imed,3,3);

printf("\nIncentre = \n");
printMat(G_I,1,2);
//contact points
printf("\ncontact points = \n");
printMat(G_i,2,3);


//********************************END OF ANGULAR BISECTOR*****************************

return 0;
}C

