import sys
sys.path.insert(0,'/sdcard/geometry/geometry/codes/CoordGeo') #for path to external scripts

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
A = np.array([-3,0]).reshape(-1,1)
B = np.array([4,2]).reshape(-1,1)
C = np.array([1,4]).reshape(-1,1)

#Triangle sides
c = LA.norm(A-B)
a = LA.norm(B-C)
b = LA.norm(C-A)
#print(a,b,c)

#Direction Vectors
m1 = dir_vec(A,B)
m2 = dir_vec(B,C)
m3 = dir_vec(C,A)
#print(m1,m2,m3)

#Line parameters
n1 = norm_vec(A,B)
c1 = n1.T@A
n2 = norm_vec(B,C)
c2 = n2.T@B
n3 = norm_vec(C,A)
c3 = n3.T@C
#print(n1,c1,n2,c2,n3,c3)

#Angles
angA = np.degrees(np.arccos((-m1.T@m3)/(c*b)))
angB = np.degrees(np.arccos((-m1.T@m2)/(c*a)))
angC = np.degrees(np.arccos((-m2.T@m3)/(a*b)))
#print(angA,angB,angC)


#Generating all lines
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)

#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
#plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
#plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')

#Labeling the coordinates
tri_coords = np.block([[A,B]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B']
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
plt.savefig('/sdcard/geometry/geometry/figs/tri.jpg')
#subprocess.run(shlex.split("termux-open ./figs/tri_sss.pdf"))
#else
plt.show()
