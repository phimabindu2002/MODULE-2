import sys                                          #for path to external scripts
sys.path.insert(0, './CoordGeo')   
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

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

# Function to calculate the contact points
#def contact(V, u, f, h):
    #A = V @ u
    #a = u.T @ A
    #b = u.T @ h + h.T @ A
    #c = h.T @ A - f
    #delta = b**2 - a * c
    #lambda1 = (-b + np.sqrt(b**2 - a * c)) / a
    #lambda2 = (-b - np.sqrt(b**2 - a * c)) / a
    #P = h + lambda1 * A
    #Q = h + lambda2 * A
    #return P, Q

# Given parameters

a1 = 0
a2 = np.sqrt(1 - a1**2)
discriminant = (a1**2 - 4 * (a1**2 - 3/4))
b1_plus = (-a1 + np.sqrt(discriminant)) / 2
b1_minus = (-a1 - np.sqrt(discriminant)) / 2
b2 = np.sqrt(1/4 - a1/2 * (2*a1 + np.sqrt(3 - 3*a1**2)))
O = np.array([0,0]).reshape(-1,1)
A = np.array([a1,a2]).reshape(-1,1)
B1 = np.array([b1_plus,b2]).reshape(-1,1) 
B2 = np.array([b1_minus,b2]).reshape(-1,1) 

# Printing the calculated values
print("a1 =",a1)
print("a2 =",a2)
print("b1_plus =", b1_plus)
print("b1_minus =", b1_minus)
print("b2 =", b2)

#Generating all lines
x_OA = line_gen(O,A)
x_OB1 = line_gen(O,B1)
x_OB2 = line_gen(O,B2)

#Plotting all lines
plt.plot(x_OA[0,:],x_OA[1,:],label='$OA$')
plt.plot(x_OB1[0,:],x_OB1[1,:],label='$OB1$')
plt.plot(x_OB2[0,:],x_OB2[1,:],label='$OB2$')

#plt.plot(x_circ[0,:],x_circ[1,:],label='$incircle$')

#Labeling the coordinates
O = O.reshape(-1,1)
A = A.reshape(-1,1)
B1 = B1.reshape(-1,1)
B2 = B2.reshape(-1,1)
tri_coords = np.block([[A,B1,B2,O]])
plt.scatter(tri_coords[0, :], tri_coords[1, :])
vert_labels = ['A','B1','B2','O']
for i, txt in enumerate(vert_labels):
    offset = 10 if txt == 'F' else -10
    plt.annotate(txt,
                 (tri_coords[0, i], tri_coords[1, i]),
                 textcoords="offset points",
                 xytext=(0, offset),
                 ha='center')

# Save the figure
plt.xlabel('$X-axis$')
plt.ylabel('$Y-axis$')
plt.legend(loc='upper left')
plt.grid() # minor
plt.axis('equal')
plt.savefig('/home/himabindu/module2/geometry/codes/CoordGeo/pythoncopy.png')
plt.show()

