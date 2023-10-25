import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

A=np.array([1,3])
B=np.array([-3,0])
C=np.array([0,4])
#finding the side lengths
a=np.linalg.norm(C-B)
b=np.linalg.norm(C-A)
c=np.linalg.norm(A-B)
print("a=",a,"b=",b,"c=",c)

#Orthogonal matrix
omat = np.array([[0,1],[-1,0]]) 

#Triangle vertices
def tri_vert(a,b,c):
  p = (a**2 + c**2-b**2 )/(2*a)
  q = np.sqrt(c**2-p**2)
  A = np.array([p,q]) 
  B = np.array([0,0]) 
  C = np.array([a,0]) 
  return  A,B,C

def dir_vec(A,B):
  return B-A

def norm_vec(A,B):
  return omat@dir_vec(A,B)
  #return np.matmul(omat, dir_vec(A,B))

#Generate line points
def line_gen(A,B):
  len =10
  dim = A.shape[0]
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(0,1,len)
  for i in range(len):
    temp1 = A + lam_1[i]*(B-A)
    x_AB[:,i]= temp1.T
  return x_AB

#finding Incentre 
t = norm_vec(B,C) 
n1 = t/np.linalg.norm(t) #unit normal vector
t = norm_vec(C,A)
n2 = t/np.linalg.norm(t)
t = norm_vec(A,B)
n3 = t/np.linalg.norm(t)

#Intersection of two lines
def line_intersect(n1,A1,n2,A2):
  N=np.block([[n1],[n2]])
  p = np.zeros(2)
  p[0] = n1@A1
  p[1] = n2@A2
  #Intersection
  P=np.linalg.inv(N)@p
  return P

I=line_intersect(n1-n3,B,n1-n2,C)  #Incentre
#Incentre
print("I = ",I)


def icircle(A,B,C):
  k1 = 1
  k2 = 1
  p = np.zeros(2)
  t = norm_vec(B,C)
  n1 = t/np.linalg.norm(t)
  t = norm_vec(C,A)
  n2 = t/np.linalg.norm(t)
  t = norm_vec(A,B)
  n3 = t/np.linalg.norm(t)
  p[0] = n1@B- k1*n2@C
  p[1] = n2@C- k2*n3@A
  r = n1@(I-B)
  return r
  
def circ_gen(O,r):
	len = 50
	theta = np.linspace(0,2*np.pi,len)
	x_circ = np.zeros((2,len))
	x_circ[0,:] = r*np.cos(theta)
	x_circ[1,:] = r*np.sin(theta)
	x_circ = (x_circ.T + O).T
	return x_circ      

#Generating all lines
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)

#Generating the incircle
r = icircle(A,B,C)
x_icirc= circ_gen(I,r)

#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')

#plt.plot(x_icirc[0,:],x_icirc[1,:],label='$incircle$')

#finding k for E_3 and F_3 and D_3
k1=((I-A)@(A-B))/((A-B)@(A-B))
k2=((I-A)@(A-C))/((A-C)@(A-C))
k3=((I-B)@(C-B))/((C-B)@(C-B))
#finding E_3 and F_3 and D_3
E3=A+(k1*(A-B))
F3=A+(k2*(A-C))
D3=B+(k3*(C-B))
#x_ID3 = line_gen(I,D3)
#plt.plot(x_ID3[0,:],x_ID3[1,:],label='$ID3$')
x_AI = line_gen(A,I)
#plt.plot(x_AI[0,:],x_AI[1,:],'--',label='$AI$')
x_BI = line_gen(B,I)
#plt.plot(x_BI[0,:],x_BI[1,:],'--',label='$BI$')
x_CI = line_gen(C,I)
#plt.plot(x_CI[0,:],x_CI[1,:],'--',label='$CI$')
print("k1 = ",k1)
print("k2 = ",k2)
print("k3 = ",k3)
print("D3 = ",D3)
print("E3 = ",E3)
print("F3 = ",F3)

r=a=np.linalg.norm(I-D3)
print("r=",r)
print(A-I)
#Labeling the coordinates
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
I = I.reshape(-1,1)
D3 = D3.reshape(-1,1)
E3 = E3.reshape(-1,1)
F3 = F3.reshape(-1,1)
tri_coords = np.block([[A,B,C,D3,E3,F3]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D3','E3','F3']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center


plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid(True) # minor
plt.axis('equal')
plt.savefig("Angle1.jpg",bbox_inches='tight')
