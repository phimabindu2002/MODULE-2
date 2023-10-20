import sys                                          #for path to external scripts
sys.path.insert(0, '/sdcard/digital-design/geometry/codes/CoordGeo')        #path to my scripts
import numpy as np
import mpmath as mp
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

D,E,F = tri_mid_pt(A,B,C)
G = (C+B)/2

#Generating the circumcircle
[O,R] = ccircle(A,B,C)
x_circ= circ_gen(O,R)

#angle BOC = 2 angle BAC
m1 = dir_vec(A,B)
m2 = dir_vec(A,C)
m3 = dir_vec(B,C)
v1 = dir_vec(O,B)
v2 = dir_vec(O,C)

print(ang_vec(m1,m2),2*mp.pi-ang_vec(v1,v2))

#OD perpendicular to BC
m3 = dir_vec(B,C)
print(m3.T@(O-D))


#Plotting the circumcircle
#plt.plot(x_circ[0,:],x_circ[1,:],label='$circumcircle$')

#Generating all lines
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_OD = line_gen(O,D)
x_OE = line_gen(E,O)
x_OF = line_gen(O,F)
x_OG = line_gen(O,G)

#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='AB')
plt.plot(x_BC[0,:],x_BC[1,:],label='BC')
plt.plot(x_CA[0,:],x_CA[1,:],label='CA')
plt.plot(x_OD[0,:],x_OD[1,:],'--',label='$OD$')
plt.plot(x_OE[0,:],x_OE[1,:],'--',label='$OE$')
plt.plot(x_OF[0,:],x_OF[1,:],'--',label='$OF$')
#plt.plot(x_OG[0,:],x_OG[1,:],label='OG')

#Labeling the coordinates
tri_coords = np.block([[A,B,C,D,O,E,F]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D','O','E','F']
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
plt.savefig('/sdcard/hima/geometry/figs/abcdefo.jpg')
#subprocess.run(shlex.split("termux-open ./figs/tri_sss.pdf"))
#else
plt.show()
