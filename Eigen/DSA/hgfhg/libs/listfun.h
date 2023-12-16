
//Function declaration
double Listnorm(avyuh *a);//norm of a vector                        //
double Listdot(avyuh *a, avyuh * b);//inner product                 //
double ListVecdot(sadish *a, sadish * b);//inner product            //
avyuh *Listsub(avyuh *a, avyuh *b);//subtract two matrices          //
avyuh *Listadd(avyuh *a, avyuh *b);//add two matrices               //
avyuh *Listscale(avyuh *a,double k);//scale matrix
sadish *Listvecscale(sadish *a,double k);//scale vector
//avyuh *Listscale(avyuh *a, int m, int n, double k);//scale matrix
avyuh *Listinv(avyuh *mat, int m);//invert an m  x m matrix, m <=3
avyuh *Listmul(avyuh *a, avyuh *b);//multiply matrices a and b      // 
avyuh *rotList(double theta); //rotation matrix                     //
avyuh *normVec(avyuh *a); //normal vector                           //
avyuh *circulantList(avyuh *c);
avyuh *Listsec(avyuh *a, avyuh * b, int m, double k);//section formula

double Listtrace(avyuh *a);
double Listdet(avyuh *a);
avyuh *Listquad(double a,double b,double c);
avyuh *Listeigval(avyuh *a);
avyuh *Listeigvec(avyuh *a);
avyuh *VertToList(avyuh *a,avyuh *b,avyuh *c);// points to list conversion
avyuh *Listdiag(avyuh *a);//diagonal matrix
avyuh *Listsqrtdiag(avyuh *a);//square root of elements in a diagonal matrix
sadish *xtrtdiag(avyuh *a);//to extract diagonal vector from a matrix
sadish *xtrtsqrtdiag(avyuh *a);//to extract sqrt diagonal vector from a matrix
avyuh *Liststack(avyuh *a,avyuh *b);//concates list a,b horizontally and transposes the result
avyuh *line_intersect(avyuh *a);//to find intersection point of lines (solving pair of linear equations)
avyuh *Listunit(avyuh *a);//unit vector
avyuh *Listrow(avyuh *a, int k);// kth row
avyuh *Listeye(int k);//identity matrix
avyuh *Listbasis(int k);//standard basis vector of length k
sadish *ListVecShift(sadish *a);//circulalry right shift vector
double readidx(avyuh *a,int m,int n);//to read the value of an element at index (m,n)
avyuh *ListContactPts(avyuh *V,avyuh *h,avyuh *u,avyuh *f);
//End function declaration
avyuh *ListContactPts(avyuh *V,avyuh *h,avyuh *u,avyuh *f)
{
	// finding circle equation   (gh = h.T@V@h+2*u.T@h+f )
	avyuh *gh=Listadd( Listadd(  Listmul( Listmul(transposeList(h),V ),h ),Listscale( Listmul( transposeList(u),h ) ,2) ),f );
	//finding sigma matrix (sigmat = (V@h+u)@(V@h+u).T-gh*V)
	avyuh *sigmat=Listsub(Listmul( Listadd(Listmul(V,h),u) , transposeList(Listadd(Listmul(V,h),u)) ) ,Listscale(V,gh->vector->data));
	// finding Eigen values and Eigen vectors
	avyuh *E_val=Listeigval(sigmat);		avyuh *P=Listeigvec(sigmat);
	//finding u1,u2
	avyuh *u1=createList(2,1);
	u1->vector->data=sqrt(fabs(E_val->next->vector->data));		u1->next->vector->data=sqrt(fabs(E_val->vector->data));
	avyuh *u2=createList(2,1);
	u2->vector->data=sqrt(fabs(E_val->next->vector->data));      	u2->next->vector->data=-sqrt(fabs(E_val->vector->data));
	//finding direction vectors m1,m2     (m = P@u)
	avyuh *m1=Listmul(P,u1); 	avyuh *m2=Listmul(P,u2);
	//finding mu1,mu2		 (mu = -(m.T@(V@h+u))/(m.T@V@m)
	double mu1=-( (Listmul(transposeList(m1),Listadd(Listmul(V,h),u)) )->vector->data) / ( (Listmul(transposeList(m1),Listmul(V,m1)) )->vector->data);
	double mu2=-( (Listmul(transposeList(m2),Listadd(Listmul(V,h),u)) )->vector->data) / ( (Listmul(transposeList(m2),Listmul(V,m2)) )->vector->data);
	//finding E,F			(x = h + mu*m)
	avyuh *E=Listadd(h,Listscale(m1,mu1)); 	avyuh *F=Listadd(h,Listscale(m2,mu2));
	return Liststack(transposeList(E),transposeList(F));	

}
//End function declaration

// kth row
avyuh *Listrow(avyuh *a, int k){
	avyuh *c= (avyuh *)malloc(sizeof(avyuh)),*tempa=a,*head;
	c->next=NULL;
	head=c;
	for(int i=0;i<k;i++){
		tempa=tempa->next;  }
	c->vector=tempa->vector;

	return transposeList(head);

}
//unit vector
avyuh *Listunit(avyuh *a){
	double k= Listnorm(a);
	avyuh *c=Listscale(a, 1/k);//scale vector
	return c;
}

//standard basis vector of length k
avyuh *Listbasis(int k){
	avyuh *head=(avyuh *)malloc(sizeof(avyuh));
	sadish *c = createVec(k);
	head->vector = c;
	head->next = NULL;
	for(int i = 0; i < k; i++){
		c->data= 0;
		c= c->next;
	}
	head->vector->data = 1;
	return head;
}
//identity matrix
avyuh *Listeye(int k){
	avyuh *c=Listbasis(k);
return circulantList(c);
}

//circulant matrix
avyuh *circulantList(avyuh *a){
	avyuh *c=(avyuh *)malloc(sizeof(avyuh));
	avyuh *head = c;
	sadish *ctemp;
	ctemp=a->vector;
	head->vector=a->vector;
	for(sadish *temp=a->vector;temp->next!=NULL;temp=temp->next){
			c->next = (avyuh *)malloc(sizeof(avyuh));
			c->next->next=NULL;
			c= c->next;
		c->vector = ListVecShift(ctemp);
		ctemp=c->vector;
	}
	return head;
}

//circulalry right shift vector
sadish *ListVecShift(sadish *a){
	sadish *tempa, *temp;
	sadish *head = ListVecopy(a);
	for(temp=head;temp->next->next!=NULL;temp=temp->next);
temp->next->next = head;
tempa = temp->next;
temp->next = NULL;
return tempa;
}






//function to find the trace of matrix
double Listtrace(avyuh *a)
{
avyuh *di=Listdiag(a);
double sumofdiag=0;

for(sadish *temp=di->vector;temp!=NULL;temp=temp->next)
	sumofdiag+=temp->data;

return sumofdiag;
}
// end of function to find the trace of matrix




//function to find det of 2x2 matrix
double Listdet(avyuh *a)
{
return ((a->vector->data*a->next->vector->next->data)-(a->vector->next->data*a->next->vector->data));
}
// end of function to find the det of 2x2 matrix




//function to find the roots of quadratic equation
avyuh *Listquad(double a, double b,double c)
{
avyuh *lam=createList(2,1);
double D = sqrt(pow(b,2.0)-4*a*c);
        double den =2.0*a;
lam->vector->data= (-b+D)/den;
lam->next->vector->data= (-b-D)/den;
return lam;

}
//end of function to find the roots of quadratic equation




//function to find the eigen values of a matrix
avyuh *Listeigval(avyuh *a)
{
double b=-Listtrace(a);
double c=Listdet(a);
return Listquad(1,b,c);
}
//end of functoin to find the eigen values




// function to find eigen vectors 
avyuh *Listeigvec(avyuh *a)
{
avyuh *lam=Listeigval(a);
printList(lam);
avyuh *omat=rotList(M_PI/2);
//A-lambdaI   approach 2
avyuh *b1=Listadd(a,Listscale(Listeye(2),-lam->vector->data));
avyuh *b2=Listadd(a,Listscale(Listeye(2),-lam->next->vector->data));
// unit vector approach 3
avyuh *c1 = Listunit(Listrow(b1,0));
avyuh *c2 = Listunit(Listrow(b2,0));
// find eigen vector
avyuh *p1=transposeList(Listmul(omat,c1));
avyuh *p2=transposeList(Listmul(omat,c2));
return Liststack(p1,p2);
}
//end of function to find eigen vectors

	
//function to scale the vector with value k
sadish *Listvecscale(sadish *a,double k)
{
sadish *head = (sadish *)malloc(sizeof(sadish)), *c, *tempa=a;
	c = head; 
	head->next = NULL;
	for(sadish *tempa=a;tempa!=NULL;tempa=tempa->next){
                c->data =k*tempa->data;
	if(tempa->next!=NULL){
		c->next = (sadish *)malloc(sizeof(sadish));
	c->next->next=NULL;
		c= c->next;
	}
	}
return head;
}
//end of  function to scale a vector




// function to scale a list with value k
avyuh *Listscale(avyuh *a,double k)
{
avyuh *c= (avyuh *)malloc(sizeof(avyuh)), *head; 
	c->next = NULL;
	head = c;
	 for(avyuh *tempa=a;tempa!=NULL;tempa=tempa->next){
                c->vector = Listvecscale(tempa->vector,k);
	if(tempa->next!=NULL){
		c->next = (avyuh *)malloc(sizeof(avyuh));
		c->next->next=NULL;
		c= c->next;
	}
	}
	return head;


}
//end of function to scale  a list 





//inner product
double ListVecdot(sadish *a, sadish *b){
	double val = 0;
	sadish  *tempb=b;
	 for(sadish *tempa=a;tempa!=NULL;tempa=tempa->next){
		val += tempa->data*tempb->data;
		//tempa = tempa->next;
		tempb = tempb->next;
	}
	return val;
}



double Listdot(avyuh *a, avyuh *b){
	return Listmul(transposeList(a),b)->vector->data ;
}



sadish *ListVecsub(sadish *a, sadish *b){
	sadish *head = (sadish *)malloc(sizeof(sadish)), *c,  *tempb=b;
	c = head; 
	head->next = NULL;
	for(sadish *tempa=a;tempa!=NULL;tempa=tempa->next){
		c->data = tempa->data-tempb->data;
		//tempa = tempa->next;
		tempb = tempb->next;
	if(tempa->next!=NULL){
		c->next = (sadish *)malloc(sizeof(sadish));
		c->next->next=NULL;
		c= c->next;
	}
	}
	return head;
}




//subtract two matrices
avyuh *Listsub(avyuh *a, avyuh *b){
	avyuh *c= (avyuh *)malloc(sizeof(avyuh)),  *tempb = b, *head; 
	c->next = NULL;
	head = c;
	 for(avyuh *tempa=a;tempa!=NULL;tempa=tempa->next){
		c->vector = ListVecsub(tempa->vector,tempb->vector);
		//tempa = tempa->next;
		tempb = tempb->next;
	if(tempa->next!=NULL){
		c->next = (avyuh *)malloc(sizeof(avyuh));
		c->next->next=NULL;
		c= c->next;
	}
	}
	return head;
}




sadish *ListVecadd(sadish *a, sadish *b){
        sadish *head = (sadish *)malloc(sizeof(sadish)), *c, *tempb=b;
        c = head; 
        head->next = NULL;
         for(sadish *tempa=a;tempa!=NULL;tempa=tempa->next){
                c->data = tempa->data+tempb->data;
                //tempa = tempa->next;
                tempb = tempb->next;
        if(tempa->next!=NULL){
                c->next = (sadish *)malloc(sizeof(sadish));
                c->next->next=NULL;
                c= c->next;
        }
        }
        return head;
}





//add two matrices
avyuh *Listadd(avyuh *a, avyuh *b){
        avyuh *c= (avyuh *)malloc(sizeof(avyuh)), *tempb = b, *head; 
        c->next = NULL;
        head = c;
         for(avyuh *tempa=a;tempa!=NULL;tempa=tempa->next){
                c->vector = ListVecadd(tempa->vector,tempb->vector);
                //tempa = tempa->next;
                tempb = tempb->next;
        if(tempa->next!=NULL){
                c->next = (avyuh *)malloc(sizeof(avyuh));
                c->next->next=NULL;
                c= c->next;
        }
        }
        return head;
}





//norm of a vector
double Listnorm(avyuh *a){
	return sqrt(Listdot(a,a));
}



//rotation matrix
avyuh *rotList(double theta){ 
avyuh *head = createList(2,2), *temp;//create empty 2 x 2 matrix 
sadish *row1, *row2;
double c = cos(theta), s = sin(theta);

row1 = head->vector;
row2 = head->next->vector;
row1->data = c;
row1->next->data = -s;
row2->data = s;

row2->next->data = c;
return head;

}




avyuh *VertToList(avyuh *a,avyuh *b,avyuh *c)
{
	int i;
avyuh *head = createList(2,3);
sadish *r1,*r2;
r1      =head->vector;
r2      =head->next->vector;

r1->data             =Vecind(a->vector,0)->data;
r1->next->data       =Vecind(b->vector,0)->data;
r1->next->next->data =Vecind(c->vector,0)->data;

r2->data             =Vecind(a->vector,1)->data;
r2->next->data       =Vecind(b->vector,1)->data;
r2->next->next->data =Vecind(c->vector,1)->data;

return head;
}




//Matrix multiplication
avyuh *Listmul(avyuh *a, avyuh *b){
	avyuh *c=(avyuh *)malloc(sizeof(avyuh)),*atemp, *btemp;
	avyuh *head = c;
	sadish *cvec=(sadish *)malloc(sizeof(sadish));
	head->next = NULL;
	for(atemp=a; atemp !=NULL; atemp=atemp->next){
		cvec=(sadish *)malloc(sizeof(sadish));
		cvec->next = NULL;
		c->vector = cvec;  
		for(btemp=transposeList(b); btemp !=NULL; btemp=btemp->next){
			cvec->data=ListVecdot(atemp->vector, btemp->vector);//inner product
		if(btemp->next!=NULL){
			cvec->next = (sadish *)malloc(sizeof(sadish));
			cvec->next->next=NULL;
			cvec= cvec->next;
		}
		}

		if(atemp->next!=NULL){
			c->next = (avyuh *)malloc(sizeof(avyuh));
			c->next->next=NULL;
			c= c->next;
		}
	}
return head;
}





//function to find the diagonal matrix
avyuh *Listdiag(avyuh *a){
	avyuh *head= (avyuh *)malloc(sizeof(avyuh)), *alist;
	head->vector= xtrtdiag(a);
	head->next= NULL;
return head;
}



//function to extract diagonal vector from the matrix
sadish *xtrtdiag(avyuh *a){
	int i = 0,j = 0;//dummy integers
	int m = 2;
	//avyuh *head= (avyuh *)malloc(sizeof(avyuh)), *alist;
	sadish *head=(sadish *)malloc(sizeof(sadish)),*temp, *btemp;
	btemp = head;
	head->next= NULL;
	for (avyuh *alist = a; alist != NULL; alist= alist->next){
		temp = Vecind(alist->vector,i);//getting address of the nth column
		btemp->data = temp->data;
		i++;
	if(alist->next !=NULL){
		btemp->next = (sadish *)malloc(sizeof(sadish));
		btemp->next->next=NULL;
		btemp= btemp->next;
	}
	}
return head;
}



// function to find sqrt diagonal matrix
avyuh *Listsqrtdiag(avyuh *a){
        avyuh *head= (avyuh *)malloc(sizeof(avyuh)), *alist;
        head->vector= xtrtsqrtdiag(a);
        head->next= NULL;
return head;   
} 



//function to extract sqrt diagonal vector from a matrix
sadish *xtrtsqrtdiag(avyuh *a){
        int i = 0,j = 0;//dummy integers
        int m = 2;
        //avyuh *head= (avyuh *)malloc(sizeof(avyuh)), *alist;
        sadish *head=(sadish *)malloc(sizeof(sadish)),*temp, *btemp;
        btemp = head;
        head->next= NULL;
        for (avyuh *alist = a; alist != NULL; alist= alist->next){
                temp = Vecind(alist->vector,i);//getting address of the nth column
                btemp->data = sqrt(temp->data);
                i++;
        if(alist->next !=NULL){
                btemp->next = (sadish *)malloc(sizeof(sadish));
                btemp->next->next=NULL;
                btemp= btemp->next;
        }
        }
return head;
}



//function to concate lists a and b horizontally and transpose the result
avyuh *Liststack(avyuh *a, avyuh *b){
	avyuh *c= (avyuh *)malloc(sizeof(avyuh)), *tempb = b, *head; 
	c->next = NULL;
	head = c;
	 for(avyuh *tempa=a;tempa!=NULL;tempa=tempa->next){
		c->vector = tempa->vector;
		//tempa = tempa->next;
	
		c->next = (avyuh *)malloc(sizeof(avyuh));
		c->next->next=NULL;
		c= c->next;

	
	}
	c->vector=b->vector;
	return transposeList(head);
}



// function to find the point of intersection of lines
avyuh  *line_intersect(avyuh *a)
{
avyuh *temp = transposeList(a);
double a1,a2,b1,b2,c1,c2;
a1=temp->vector->data;
a2=temp->vector->next->data;
b1=temp->next->vector->data;
b2=temp->next->vector->next->data;
c1=-temp->next->next->vector->data;
c2=-temp->next->next->vector->next->data;
//printf("%lf,%lf,%lf,%lf,%lf,%lf,",a1,a2,b1,b2,c1,c2);
avyuh *head = createList(1,2);
head->vector->data=(b1*c2-b2*c1)/(a1*b2-a2*b1);
head->vector->next->data=(c1*a2-c2*a1)/(a1*b2-a2*b1);
return head;

}


/*
//function to read the value of element at index (m,n)
double readidx(avyuh *a,int m,int n)
{
avyuh *head=a;

for(int i=0;i<m;i++)
{
head=head->next;
}
//for (int j=0;j<n;j++)
//{
//head->vector=head->vector->next;
//}
double val = Vecind(head->vector,n)->data;
return val;
}
*/
