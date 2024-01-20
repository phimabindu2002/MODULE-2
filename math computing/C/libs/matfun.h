//Functions created by
// G V V Sharma
// October 27, 2023

//Function declaration
double **createMat(int m,int n);//create m x n matrix array
void printMat(double **p,int m,int n);//print matrix
double **loadMat(char *str,int m,int n);//load matrix from file
double Matnorm(double **a, int m);//norm of a vector
double Matdot(double **a, double ** b, int m);//inner product
double **Matsub(double **a, double **b, int m, int n);//subtract two matrices
double **Matadd(double **a, double **b, int m, int n);//add two matrices
double **Matscale(double **a, int m, int n, double k);//scale matrix
double **Matinv(double **mat, int m);//invert an m  x m matrix, m <=3
double **Matmul(double **a, double **b, int m, int n, int p);//multiply matrices a and b
double **transposeMat(double **a,  int m, int n);//transpose of a
double **rotMat(double theta); //rotation matrix
double **normVec(double **a); //normal vector
void circulantMat(double **c, int m);
double **Matsec(double **a, double ** b, int m, double k);//section formula
double **matrix_merge(double **m1,double **m2,double **m3,int m, int n);
double **matrix(double c1,double c2,double c3);

//section formula
double **Matsec(double **a, double ** b, int m, double k){
	double **temp=createMat(m,1);
	temp = Matscale(Matadd(a,Matscale(b,m,1,k),m,1),m,1,1/(k+1));
	return temp;
}

double **matrix(double c1,double c2,double c3){
     double **c;
     c=createMat(1,3);
     c[0][0]=c1;
     c[0][1]=c2;
     c[0][2]=c3;
     return c;
}

//add matrices
double **Matadd(double **a,double **b, int m, int n){
int i, j;
double **c;
c = createMat(m,n);

 for(i=0;i<m;i++)
 {
  for(j=0;j<n;j++)
  {
c[i][j]= a[i][j]+b[i][j];
  }
 }
return c;
}

//scale matrix
double **Matscale(double **a, int m, int n, double k){
int i, j;
double **c;
c = createMat(m,n);

 for(i=0;i<m;i++)
 {
  for(j=0;j<n;j++)
  {
c[i][j]= k*a[i][j];
  }
 }
return c;
}

//Generating a circulant matrix from a vector
void circulantMat(double **c, int m){
    int i,j,k;
 
    // Forming the circulant matrix
    for (int i = 1; i <= m - 1; i++) {
        for (int j = 0; j <= m - 1; j++) {
            if (j - 1 >= 0)
                c[j][i] = c[j - 1][i - 1];
            else
                c[j][i] = c[m - 1][i - 1];
        }
    }
}
//Obtaining the normal vector
double **normVec(double **m){
	double **temp;
	temp = Matmul(rotMat(M_PI/2),m,2,2,1);
	return temp;
}

//Defining the function for matrix creation
double **createMat(int m,int n)
{
 int i;
 double **a;
 
 //Allocate memory to the pointer
a = (double **)malloc(m * sizeof( *a));
    for (i=0; i<m; i++)
         a[i] = (double *)malloc(n * sizeof( *a[i]));

 return a;
}
//End function for matrix creation

//Extract column
//
double **Matcol(double **a,int m, int n){
	int i = 0;
	double **b = createMat(m,1);//create column with m rows

//extract column vector
	for (i = 0; i < m; i++){
		b[i][0] = a[i][n];
	}
return b;
}



//Read  matrix from file
double **loadMat(char *str,int m,int n)
{
FILE *fp;
double **a;
int i,j;


a = createMat(m,n);
fp = fopen(str, "r");

 for(i=0;i<m;i++)
 {
  for(j=0;j<n;j++)
  {
   fscanf(fp,"%lf",&a[i][j]);
  }
 }
//End function for reading matrix from file

fclose(fp);
 return a;

}


//Defining the function for printing
void printMat(double **p, int m,int n)
{
 int i,j;

 for(i=0;i<m;i++)
 {
  for(j=0;j<n;j++)
  printf("%lf ",p[i][j]);
 printf("\n");
 }
}
//End function for printing

//Rotation matrix

double **rotMat(double theta){
double **temp=createMat(2,2);//creating the matrix
double c = cos(theta), s = sin(theta);
temp[0][0] = c;
temp[0][1] = -s;
temp[1][0] = s;
temp[1][1] = c;

return temp;
}
//inner product
double Matdot(double **a, double ** b, int m){
	double **temp= Matmul(transposeMat(a,m,1),b,1,m,1);
	return temp[0][0];
}
//Defining the function for norm

double Matnorm(double **a, int m){
	return sqrt(Matdot(a, a, m));
}
//Defining the function for difference of matrices

double **Matsub(double **a, double **b, int m, int n)
{
int i, j;
double **c;
c = createMat(m,n);

 for(i=0;i<m;i++)
 {
  for(j=0;j<n;j++)
  {
c[i][j]= a[i][j]-b[i][j];
  }
 }
return c;

}
//End function for difference of matrices

//Defining the function for inverse of a matrix
//code adapted from the internet


double **Matinv(double **a, int m)
{
double **c, det=0;
int i,j;
c = createMat(m,m);
printMat(c,m,m);
if (m==2){
det = a[0][0]*a[1][1]-a[0][1]*a[1][0];

c[0][0] = a[1][1]/det;
c[0][1] = -a[1][0]/det;
c[1][0] = -a[0][1]/det;
c[1][1] = a[0][0]/det;
}
else if(m==3){
for(i=0;i<m;i++)
      det += a[0][i]*(a[1][(i+1)%3]*a[2][(i+2)%3] - a[1][(i+2)%3]*a[2][(i+1)%3]);
 
   for(i=0;i<m;i++){
      for(j=0;j<m;j++)
	   c[i][j]=((a[(i+1)%3][(j+1)%3] * a[(i+2)%3][(j+2)%3]) - (a[(i+1)%3][(j+2)%3]*a[(i+2)%3][(j+1)%3]))/det;
   }
}
else {
	printf("Invalid input \n");
	exit(0);
}
 
return c;
}
// End  function for inverse of 2x2 matrix


//Defining the function for product of matrices

double **Matmul(double **a, double **b, int m, int n, int p)
{
int i, j, k;
double **c, temp =0;
c = createMat(m,p);

 for(i=0;i<m;i++)
 {
  for(k=0;k<p;k++)
  {
    for(j=0;j<n;j++)
    {
	temp= temp+a[i][j]*b[j][k];
    }
	c[i][k]=temp;
	temp = 0;
  }
 }
return c;

}
//End function for difference of matrices

//Defining the function for transpose of matrix

double **transposeMat(double **a,  int m, int n)
{
int i, j;
double **c;
//printf("I am here");
c = createMat(n,m);

 for(i=0;i<n;i++)
 {
  for(j=0;j<m;j++)
  {
c[i][j]= a[j][i];
//  printf("%lf ",c[i][j]);
  }
 }
return c;

}
//End function for transpose of matrix


//printing of three matrices
double **matrix_merge(double **m1,double **m2,double **m3,int m, int n){
        double **resultant;
       resultant = createMat(2,3);
        for (int i=0; i<m;i++){ 
		for (int j=0;j<n;j++){ if(j==0)
                	resultant[i][j]=*m1[i];
                        if(j==1)
                        resultant[i][j]=*m2[i];
                        if(j==2)
                        resultant[i][j]=*m3[i];
                        }
                }
        return resultant;
        }

//merging Two_matrices 
double **matrix_2Merge(double **mat1, int rows1, int cols1, double **mat2, int rows2, int cols2){
    int mergedCols = cols1 + cols2;
    double **resultant = createMat(rows1, mergedCols);
    for (int i = 0; i < rows1; i++) {
        for (int j = 0; j < cols1; j++) {
            resultant[i][j] = mat1[i][j];
        }
        for (int j = 0; j < cols2; j++) {
            resultant[i][j + cols1] = mat2[i][j];
        }
    }
    return resultant;
}

// defining a function for finding the intersection point between two lines
double **line_intersect(double **a,int m,int n)
{
double **c=createMat(1,2);
a[0][2]=-a[0][2];
a[1][2]=-a[1][2];
c[0][0]= (a[0][1]*a[1][2]-a[1][1]*a[0][2])/(a[0][0]*a[1][1]-a[1][0]*a[0][1]);
c[0][1]= (a[0][2]*a[1][0]-a[1][2]*a[0][0])/(a[0][0]*a[1][1]-a[1][0]*a[0][1]);

return c;
}

