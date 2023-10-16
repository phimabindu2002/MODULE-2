import sys                                          #for path to external scripts
sys.path.insert(0, '/sdcard/hima/geometry/codes/CoordGeo')   #path to my scripts
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
A = np.array([-3,-5]).reshape(-1,1)
B = np.array([3,-5]).reshape(-1,1) 
C = np.array([-4,-3]).reshape(-1,1) 

#Triangle  mid points
D = (B+C)/2
E = (C+A)/2
F = (A+B)/2
#print(D,E,F)

#Median parameters
n1 = norm_vec(A,D)
c1 = n1.T@A
n2 = norm_vec(B,E)
c2 = n2.T@B
n3 = norm_vec(C,F)
c3 = n3.T@C
print(n1,c1,n2,c2,n3,c3)

#Intersecton of BE and CF
G = line_intersect(n2,B,n3,C)

#Centroid
#print(G, (A+B+C)/3)

#Collinearity check
mat=np.block([[1,1,1],[A,G,D]]).T
print(LA.matrix_rank(mat))

#Generating all lines
x_AD = line_gen(A,D)
x_BE = line_gen(B,E)
x_CF = line_gen(C,F)
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)

#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_AD[0,:],x_AD[1,:],label='$AD$')
plt.plot(x_BE[0,:],x_BE[1,:],label='$BE$')
plt.plot(x_CF[0,:],x_CF[1,:],label='$CF$')

#Labeling the coordinates
tri_coords = np.block([[A,B,C,G]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','G']
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
plt.savefig('/sdcard/hima/gemoetry/figs/median2.pdf')
#subprocess.run(shlex.split("termux-open ./figs/tri_sss.pdf"))
#else
plt.show()
