import sys                                          #for path to external scripts
sys.path.insert(0, './CoordGeo')
import numpy as np
import matplotlib.pyplot as plt
import math

# local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen

theta2 = math.pi/3
theta1 = theta2+(2*math.pi/3)

x1, y1 = round(np.cos(theta1),3), round(np.sin(theta1),3)
x2, y2 = round(np.cos(theta2),3), round(np.sin(theta2),3)

a = np.array([x1, y1])
b = np.array([x2, y2])
O = np.array([0, 0])
c = a + b
norm = round(LA.norm(c),3)
print("a=",a)
print("b=",b)
print("c=",c)
#print("o=",O)
print("||a+b|| = ",norm)

# Plotting the unit circle
theta = np.linspace(0, 2 * np.pi, 100)
x_circle = np.cos(theta)
y_circle = np.sin(theta)
plt.plot(x_circle, y_circle, label='Unit Circle')

#Generating all lines
x_Oa = line_gen(O,a)
x_Ob = line_gen(O,b)
x_Oc = line_gen(O,c)


#Plotting all lines
plt.plot(x_Oa[0,:],x_Oa[1,:],label='$Oa$')
plt.plot(x_Ob[0,:],x_Ob[1,:],label='$Ob$')
plt.plot(x_Oc[0,:],x_Oc[1,:],label='$Oc$')

#Labeling the coordinates
a = a.reshape(-1,1)
b = b.reshape(-1,1)
c = c.reshape(-1,1)
O = O.reshape(-1,1)

tri_coords = np.block([[a, b, c, O]])
plt.scatter(tri_coords[0, :], tri_coords[1, :])
vert_labels = ['a', 'b', 'c','O']
for i, txt in enumerate(vert_labels):
    offset = 10 if txt == 'O' else -10
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
# Adding grid and legend
plt.legend()
plt.savefig("python_plot.png")
plt.show()
