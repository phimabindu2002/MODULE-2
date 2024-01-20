#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libs/listgen.h"
#include "libs/listfun.h"


int  main()
{
avyuh *a,*b,*c; //lists a,b,c
double x1, y1, x2, y2, norm;
int m =2, n=1;
double theta2 = M_PI/3;
double theta1 = theta2+(2*M_PI/3);
    x1 = cos(theta1);
    y1 = sin(theta1);
    x2 = cos(theta2);
    y2 = sin(theta2);
   
//loading values of 'a' into list    
    a = createList(m,n); 
    a->vector->data=x1;
    a->next->vector->data=y1;

//loading values of 'b' into list
    b = createList(m,n); 
    b->vector->data=x2;
    b->next->vector->data=y2;

//adding a,b storing result into 'c'
c = Listadd(a, b);

//finding the ||a+b||
norm =Listnorm(transposeList(c));

//printing the lists a,b,c
printf("a =\n"),printList(a),printf("\n");
printf("b =\n"),printList(b),printf("\n");
printf("c =\n"),printList(c),printf("\n");
printf("||a+b|| = %lf \n",norm);
    return 0;
}
