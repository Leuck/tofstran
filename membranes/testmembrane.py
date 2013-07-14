#! /usr/bin/env python3
import numpy as np
from membranesolver import solvemembrane
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
# data to define a test membrane
# structure geometry
nodes = np.array([
    #[ x, y]
    [ 0,-1 ],
    [ 1,-1 ],
    [ 0,0 ],
    [ 1,0 ],
    [ 0,1 ],
    [ 1,1 ]
    ])
elements = np.array([
    #[ node a, node b, node c, property set number]
    [ 0, 1, 3, 0 ],
    [ 0, 3, 2, 1 ],
    [ 2, 3, 5, 0 ],
    [ 2, 5, 4, 0 ]
    ])

# material and section properties
properties = np.array([
    #[ k ] 
    [ 2e1 ],
    [ 4e1 ]
    ])

# contour conditions
loads = np.array([
    #[ node number, fz]
    [ 2, 1 ],
    [ 3, 3 ]
    ])
nan = float("nan")
fixities = np.array([
    #[ node number, z]
    # nan means free
    [ 0, -0.02 ],
    [ 1, 0 ],
    [ 4, 0.01 ],
    [ 5, 0 ]
    ])

R = solvemembrane(nodes, elements, properties, loads, fixities )

print("\n== DISPLACEMENTS ==\n")
Dis = np.split(R, len(nodes) )
for i in range(len(Dis)):
    print("Node", i, "\t", Dis[i][0][0] )

print("\n== MAX/MIN DISPLACEMENTS ==\n")
print(R.max(),"\t",R.min(),"\n")

# show results in graphical window
x=nodes[:,0]
y=nodes[:,1]
z=R[:,0]
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.4)
fig.colorbar(surf)
plt.show()

#print("\n== ELEMENT FORCES ==\n")
#for i in list(range(len(R[1]))):
    #print("Element", i, "\n", R[1][i] )

#print("\n== NORMAL STRESSES ==\nA\t\t    B\n")
#for i in list(range(len(R[2]))):
    #print("Element", i, "\n", R[2][i] )
    
# vim: ts=4 et sw=4 sts=4 ai
