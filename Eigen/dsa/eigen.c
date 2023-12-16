#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib/listgen.h"
#include "lib/listfun.h"
int  main()
{
	//******************Eigen value approach to find E,F********************************
	//Vertex-A
	avyuh *h= loadList("A.dat", 2, 1);
	//Incentre
	avyuh *u=loadList("I.dat",2,1);	
	//Inradius
	double r= 1.376;
	//U=-I
	u->vector->data=-u->vector->data;
	u->next->vector->data=-u->next->vector->data;  
	//2X2 Identity matrix
	avyuh *V=Listeye(2);			  
	avyuh *f=createList(1,1);
	f->vector->data=pow(Listnorm(u),2)-r*r;	//     from circle equatin x^2+y^2=f

	// finding circle equation   (gh = h.T@V@h+2*u.T@h+f )
	avyuh *gh=Listadd( Listadd(  Listmul( Listmul(transposeList(h),V ),h ),Listscale( Listmul( transposeList(u),h ) ,2) ),f );
	//finding sigma matrix (sigmat = (V@h+u)@(V@h+u).T-gh*V)
	avyuh *sigmat=Listsub(Listmul( Listadd(Listmul(V,h),u) , transposeList(Listadd(Listmul(V,h),u)) ) ,Listscale(V,gh->vector->data));
	// finding Eigen values and Eigen vectors
	avyuh *E_val=Listeigval(sigmat);
	avyuh *P=Listeigvec(sigmat);
	//finding u1,u2
	avyuh *u1=createList(2,1);
	u1->vector->data=sqrt(fabs(E_val->next->vector->data));
	u1->next->vector->data=sqrt(fabs(E_val->vector->data));
	avyuh *u2=createList(2,1);
	u2->vector->data=sqrt(fabs(E_val->next->vector->data));
	u2->next->vector->data=-sqrt(fabs(E_val->vector->data));
	//finding direction vectors m1,m2     (m = P@u)
	avyuh *m1=Listmul(P,u1);
	avyuh *m2=Listmul(P,u2);
	//finding mu1,mu2		 (mu = -(m.T@(V@h+u))/(m.T@V@m)
	double mu1=-( (Listmul(transposeList(m1),Listadd(Listmul(V,h),u)) )->vector->data) / ( (Listmul(transposeList(m1),Listmul(V,m1)) )->vector->data);
	double mu2=-( (Listmul(transposeList(m2),Listadd(Listmul(V,h),u)) )->vector->data) / ( (Listmul(transposeList(m2),Listmul(V,m2)) )->vector->data);
	//finding E,F			(x = h + mu*m)
	avyuh *E=Listadd(h,Listscale(m1,mu1));
	avyuh *F=Listadd(h,Listscale(m2,mu2));

	// Printing lists . 	
	printf("V = \n");	
	printList(V);
	printf("h = \n");	
	printList(h);
	printf("u = \n");	
	printList(u);
	printf("f = \n");	
	printList(f);

	printf("E= \n");             
	printList(E);
	printf("F= \n");             
	printList(F);
return 0;
}
