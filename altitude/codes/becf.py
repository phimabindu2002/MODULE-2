import sys                                          #for path to external scripts
sys.path.insert(0, '/sdcard/digital-design/geometry/Codes/CoordGeo')        #path to my scripts
import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen

#if using termux
import subprocess
import shlex
#end if

#Triangle vertices
A = np.array([1,3]).reshape(-1,1)
B = np.array([-3,0]).reshape(-1,1) 
C = np.array([0,4]).reshape(-1,1) 

#Side parameters
m1 = dir_vec(A,B)
m2 = dir_vec(B,C)
m3 = dir_vec(C,A)
n1 = norm_vec(A,B)
n2 = norm_vec(B,C)
n3 = norm_vec(C,A)

#Altitude foot
D = line_intersect(m2,A,n2,B)
E = line_intersect(m3,B,n3,C)
F = line_intersect(m1,C,n1,A)
#print(D,E,F)

#Intersecton of BE and CF
H = line_intersect(m3,B,m1,C)

#Orthocentre
print(H)
#Orthogonality check
print((A-H).T@(B-C))

#Generating all lines
x_AD = line_gen(A,D)
x_BE = line_gen(B,E)
x_CF = line_gen(C,F)
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_AH = line_gen(H,A)
x_AE = line_gen(E,A)
x_AF = line_gen(A,F)
x_EH = line_gen(E,H)
x_FH = line_gen(F,H)

#Plotting all lines
#plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
#plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
#plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
#plt.plot(x_AD[0,:],x_AD[1,:],label='$AD$')
plt.plot(x_BE[0,:],x_BE[1,:],label='$BE$')
plt.plot(x_CF[0,:],x_CF[1,:],label='$CF$')
#plt.plot(x_AH[0,:],x_AH[1,:],'--',label='$AH$')
#plt.plot(x_AE[0,:],x_AE[1,:],'--',label='$BE$')
#plt.plot(x_AF[0,:],x_AF[1,:],'--',label='$CF$')
#plt.plot(x_EH[0,:],x_EH[1,:],'--',label='$EH$')
#plt.plot(x_FH[0,:],x_FH[1,:],'--',label='$FH$')

#Labeling the coordinates
tri_coords = np.block([[B,E,C,F]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['B','E','C','F']
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
plt.savefig('/sdcard/digital-design/geometry/figs/becf.jpg')
#subprocess.run(shlex.split("termux-open ./figs/tri_sss.pdf"))
#else
plt.show()
