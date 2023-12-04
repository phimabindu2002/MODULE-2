#Code by GVV Sharma
#June 12, 2022
#Revised November 14, 2023
#released under GNU GPL
#https://www.gnu.org/licenses/gpl-3.0.en.html

#Drawing a pair of tangents to a conic

import sys                                          #for path to external scripts
sys.path.insert(0, '/sdcard/Download/geometry/codes/CoordGeo')        #path to my scripts
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import *

#if using termux
import subprocess
import shlex
#end if



#Triangle vertices
A = np.array([1,3]).reshape(-1,1)
B = np.array([-3,0]).reshape(-1,1) 
C = np.array([0,4]).reshape(-1,1) 

#Incircle parameters
h  = A
V = np.eye(2)
[I,r] = icircle(A,B,C)
u = -I
f = LA.norm(I)**2-r**2

#eigen values
[x1,x2] = contact(V,u,f,h)

#printing Eigen VAlues
print("E=",x1)
print("F=",x2)

#Generating all lines
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)

#generating incircle
x_circ = circ_gen(I,r)

#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_circ[0,:],x_circ[1,:],label='$incircle$')

#Labeling the coordinates
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
E = x1.reshape(-1,1)
F = x2.reshape(-1,1)

tri_coords = np.block([[A, B, C,E,F]])
plt.scatter(tri_coords[0, :], tri_coords[1, :])
vert_labels = ['A', 'B', 'C','E','F']
for i, txt in enumerate(vert_labels):
    offset = 10 if txt == 'F' else -10
    plt.annotate(txt,
                 (tri_coords[0, i], tri_coords[1, i]),
                 textcoords="offset points",
                 xytext=(0, offset),
                 ha='center')
plt.xlabel('$X-axis$')
plt.ylabel('$Y-axis$')
plt.legend(loc='upper left')
plt.grid() # minor
plt.axis('equal')
plt.savefig('/sdcard/Download/geometry/codes/CoordGeo/eigen.png')
    
