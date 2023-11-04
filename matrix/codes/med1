#Code by GVV Sharma
#September 7, 2023
#Revised October 1, 2023
#Revised October 2, 2023
#released under GNU GPL
#Matrix Algebra


import sys                                          #for path to external scripts
sys.path.insert(0, './CoordGeo')        #path to my scripts
import numpy as np
import numpy.linalg as LA
import scipy.linalg as SA
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pandas as pd

#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen


#if using termux
import subprocess
import shlex
#end if



#-----------------Vectors-------------------------------
print("vectors begin","\n")
#Input parameters from excel file
df= pd.read_excel('vertices.xlsx')
#print(df ,"\n")

#Triangle Vertices
G_v= df.to_numpy()[:,:]
print("G_v=",G_v,"\n")

#Direction vector circulant matrix
C_m= SA.circulant([1,0,-1]).T
print("C_m=", C_m ,"\n")


#Direction vector Matrix
G_dir = G_v@C_m
print("G_dir=",G_dir,"\n")

#Normal vector matrix
G_n = R_o@G_v@C_m
print("G_n=",G_n,"\n")

#Find the line constants
cmat = np.diag(G_n.T@G_v).reshape(-1,1)

#line matrix
linmat = np.block([G_n.T,cmat])
print("linmat=",linmat,"\n")

#sides vector
dis = np.linalg.norm(G_dir, axis=0).reshape(-1,1)
print("dis=",dis,"\n")

'''
#Finding the angles of the triangle
dmat = np.diag(1/d)
G_dnorm = G_dir@dmat
G_dgram = G_dnorm.T@G_dnorm
#print(np.degrees(np.arccos(G_dgram)))
'''
print("vector ends","\n")
#-----------------Vectors Ends-------------------------------

#-----------------Medians-------------------------------
print("Medians begin","\n")
#Median circulant matrix
C_mid = SA.circulant([0,1,1]).T
print("C_mid=",C_mid,"\n")

#Mid point matrix
G_mid = 0.5*G_v@C_mid
print("G_mid=",G_mid,"\n")

#Median direction circulant matrix
C_mid_dir = SA.circulant([1,-0.5,-0.5])

#Median direction matrix
G_med_dir = G_v @ C_mid_dir
print("G_med_dir=",G_med_dir,"\n")

#Normal vector matrix
G_n_med = R_o@G_med_dir
print("G_n_med=",G_n_med,"\n") 
#Find the line constants
cmat_med = np.diag(G_n_med.T@G_v).reshape(-1,1)

#median  matrix
linmat_med = np.block([G_n_med.T,cmat_med])
print("linemat_med=",linmat_med ,"\n")

#Find the centroid
G_G=LA.lstsq(G_n_med.T,cmat_med,rcond=None)
print("G_G=",G_G,"\n")
#print(LA.lstsq(G_n_med.T,cmat_med),"\n")
print("median ends","\n")
#-----------------Median Ends-------------------------------

#-----------------Altitude-------------------------------
print("Altitude begins","\n")
#Circulant matrix
C_alt= SA.circulant([0,-1,1]).T
print("C_alt=",C_alt ,"\n")

#Normal Matrix
G_dir_alt = G_v@C_alt
print("G_dir_alt=",G_dir_alt,"\n")

#Find the line constants
cmat_alt = np.diag(G_dir_alt.T@G_v).reshape(-1,1)
print( np.block([G_dir_alt.T,cmat_alt]),"\n")

#altitude matrix
linmat_alt= np.block([G_dir_alt.T,cmat_alt])
print("linemat_alt",linmat_alt,"\n")

#Find the orthocentre
G_H=LA.lstsq(G_dir_alt.T,cmat_alt,rcond=None)
print("G_H=",G_H,"\n")
print("altitude ends","\n")


#-----------------Altitude Ends-------------------------------

#-----------------Perpendicular Bisector-------------------------------
print("perpendicular bisector begins","\n")
#Find the line constants
cmat_perp_bis= np.diag(G_dir_alt.T@G_mid).reshape(-1,1)
print( np.block([G_dir_alt.T,cmat_perp_bis]))


#Find the Circumcentre
#G_O = A.lstsq(G_dir_alt.T,cmat_perp_bis)
G_O = LA.lstsq(G_dir_alt.T, cmat_perp_bis, rcond=None)
print("G_O=",G_O,"\n")
#print("\n",A.lstsq(G_dir_alt.T,cmat_perp_bis))
print("Perpendicular bisector ends","\n")
#-----------------Perpendicular Bisector  Ends-------------------------------

#-----------------Angular Bisector-------------------------------
print("Angular bisector begins","\n")
#Incircle circulant matrix
C_in = SA.circulant([1,1,0]).T
print("c_in=",C_in,"\n")
#m,n,p
secvec = LA.inv(C_in)@dis
print("secvec",secvec,"\n")
#orignal
#cont_mat =np.array([np.block([secvec[1]/dis[1],secvec[0]/dis[2],0]), np.block([0, secvec[2]/dis[2],secvec[1]/dis[0]]),np.block([secvec[2]/dis[1],0,secvec[0]/dis[0]])])
#my
cont_mat =np.array([np.block([0,secvec[2]/dis[2],secvec[1]/dis[0]]), np.block([ secvec[2]/dis[1],0,secvec[0]/dis[0]]),np.block([secvec[1]/dis[1],secvec[0]/dis[2],0])])
print("cont_mat",cont_mat,"\n")
G_incir = G_v @ cont_mat
print("incircle=" ,G_incir)
#np.block(np.block([secvec[1]/dis[1],secvec[0]/dis[2],0]), np.block([0, secvec[2]/dis[2],secvec[1]/dis[0]]),np.block([secvec[2]/dis[1],0,secvec[0]/dis[0]]))
#print(np.array([np.block([secvec[1]/dis[1],secvec[0]/dis[2],0]), np.block([0, secvec[2]/dis[2],secvec[1]/dis[0]]),np.block([secvec[2]/dis[1],0,secvec[0]/dis[0]])]))
#cont_mat = np.array([secvec[1]/dis[1],secvec[0]/dis[2],0],dtype=object)
#print("\n",cont_mat)
#print("\n",secvec[1]/dis[0])
#print("\n",C_in,"\n","\n",C_in.T)
tvec = np.array([1,2,3]).reshape(-1,1)
tC = SA.circulant([0,1,0])
#print("\n",tC,"\n",tvec)
#print(tC@tvec)

#for Incentre
G_iv= G_incir
C_imid = SA.circulant([0,1,1]).T
G_imid = 0.5*G_iv@C_imid
C_ialt= SA.circulant([0,-1,1]).T
G_idir_alt = G_iv@C_ialt
cmat_iperp_bis= np.diag(G_idir_alt.T@G_imid).reshape(-1,1)
print (cmat_iperp_bis,"\n")
G_iO = LA.lstsq(G_idir_alt.T, cmat_iperp_bis, rcond=None)
print("Angular bisector ends","\n")

#i_cont_mat =np.array([np.block([-1,secvec[2]/dis[2],secvec[1]/dis[0]]), np.block([ secvec[2]/dis[1],-1,secvec[0]/dis[0]]),np.block([secvec[1]/dis[1],secvec[0]/dis[2],-1])])
##i_con =np.array([np.block([secvec[2]/dis[2],secvec[1]/dis[0],-1]), np.block([ -1,secvec[0]/dis[0],secvec[2]/dis[1]]),np.block([secvec[0]/dis[2],-1,secvec[1]/dis[1]])])
i_con =np.array([np.block([-1,secvec[2]/dis[2],secvec[1]/dis[0]]), np.block([ secvec[2]/dis[1],-1,secvec[0]/dis[0]]),np.block([secvec[1]/dis[1],secvec[0]/dis[2],-1])])
print("i_con",i_con,"\n")
#i_dir=G_v @i_con
i_dir=G_v@i_con
print("i_dir",i_dir,"\n")
i_nor=R_o@i_dir
print("i_nor",i_nor,"\n")
i_cof=np.diag(i_nor.T@G_v).reshape(-1,1)
print("i_cof",i_cof,"\n")
G_I = LA.lstsq(i_nor.T,i_cof, rcond=None)
print(G_I,"\n")



#------------------------------Angular Bisector Ends-------------------------------


#vertices
A = G_v[:, 0]
B = G_v[:, 1]
C = G_v[:, 2]
#midpoints
D = G_mid[:, 0]
E = G_mid[:, 1]
F = G_mid[:, 2]
#centroid
G = G_G[0]
#Orthocentre
H = G_H[0]
#circumcentre
O = G_O[0]
#incentre
I = G_iO[0]
Ic= G_I[0]
#incircle contact points
D3=G_incir[:,0]
E3=G_incir[:,1]
F3=G_incir[:,2]
print("\n","A=",A,"B=",B,"C=",C,"D=",D,"E=",E,"F=",F,"G=",G,"H=",H,"I=",I,"O=",O,"D3=",D3,"E3=",E3,"F3=",F3)

def midpoint(P, Q):
    return (P + Q) / 2  
#normal vector 
def norm_vec(A,B):
  omat = np.array([[0,1],[-1,0]]) 
  return omat.T@(A-B)
#to find the coefficients and constant of the equation of perpendicular bisector of BC
def perpendicular_bisector(B, C):
    midBC=midpoint(B,C)
    dir=B-C
    constant = -dir.T @ midBC
    return dir,constant
equation_coeff1,const1 = perpendicular_bisector(A, B)
equation_coeff2,const2 = perpendicular_bisector(B, C)
equation_coeff3,const3 = perpendicular_bisector(C, A)
def line_dir_pt(m, A, k1=0, k2=1):
    len = 10
    dim = A.shape[0]
    x_AB = np.zeros((dim, len))
    lam_1 = np.linspace(k1, k2, len)
    for i in range(len):
        temp1 = A + lam_1[i] * m
        x_AB[:, i] = temp1.T
    return x_AB
# Calculate the perpendicular vector and plot arrows
def perpendicular(B, C, label):
    perpendicular=norm_vec(B,C)
    mid = midpoint(B, C)
    x_D = line_dir_pt(perpendicular, mid, 0, 1)
    plt.arrow(mid[0], mid[1], perpendicular[0], perpendicular[1], color='blue', head_width=0.4, head_length=0.4, label=label)
    plt.arrow(mid[0], mid[1], -perpendicular[0], -perpendicular[1], color='blue', head_width=0.4, head_length=0.4)
    return x_D

# for Perpendicular Bisector
#x_D = perpendicular(A, B, 'OD')
#x_E = perpendicular(B, C, 'OE')
#x_F = perpendicular(C, A, 'OF')

A=A.reshape(-1,1)
B=B.reshape(-1,1)
C=C.reshape(-1,1)
D=D.reshape(-1,1)
E=E.reshape(-1,1)
F=F.reshape(-1,1)
G=G.reshape(-1,1)
H=H.reshape(-1,1)
Ic=Ic.reshape(-1,1)
D3=D3.reshape(-1,1)
E3=E3.reshape(-1,1)
F3=F3.reshape(-1,1)
#print("A=",A,"B=",B,"C=",C,"D=",D,"E=",E,"F=",F,"G=",G,"H=",H,"O=",O,"D3=",D3,"E3=",E3,"F3=",F3)

radius = np.linalg.norm(A-O)
iradius=np.linalg.norm(D3-I)

#Generating the circumcirclecircle
[O,r] = ccircle(A,B,C)
x_ccirc= circ_gen(O,radius)
[I,ir] = ccircle(D3,E3,F3)
x_icirc= circ_gen(I,iradius)


#Generating all lines 
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A) 
x_AD = line_gen(A,D)
x_BE = line_gen(B,E)
x_CF = line_gen(C,F)
x_OA = line_gen(O,A)
x_IF3= line_gen(I,F3)
x_AD3= line_gen(A,D3)
x_BE3= line_gen(B,E3)
x_CF3= line_gen(C,F3)




#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
#plt.plot(x_AD[0,:],x_AD[1,:],label='$AD$')
#plt.plot(x_BE[0,:],x_BE[1,:],label='$BE$')
#plt.plot(x_CF[0,:],x_CF[1,:],label='$CF$')
#plt.plot(x_OA[0,:],x_OA[1,:],label='$OA$')
#plt.plot(x_IF3[0,:],x_IF3[1,:],label='$IF3$')
plt.plot(x_AD3[0,:],x_AD3[1,:],'--',label='$AD3$')
plt.plot(x_BE3[0,:],x_BE3[1,:],'--',label='$BE3$')
plt.plot(x_CF3[0,:],x_CF3[1,:],'--',label='$CF3$')


#Plotting the circumcircle
#plt.plot(x_ccirc[0,:],x_ccirc[1,:],label='$circumcircle$')
plt.plot(x_icirc[0,:],x_icirc[1,:],label='$incircle$')



#Labeling the coordinates
tri_coords = np.block([[A,B,C,D3,E3,F3,I,Ic]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D3','E3','F3','I','Ic']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')

#if using termux
plt.savefig('incentre.png')
#subprocess.run(shlex.split("termux-open ./figs/tri_sss.pdf"))
#else
plt.show()
