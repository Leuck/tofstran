#! /usr/bin/env python3
import numpy as np
from membranesolver import solvemembrane
from membranesolver import readmsh
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
# data to define a test membrane
# structure geometry
#nodes = np.array([
    #[ x, y]
    #[ 0,0 ],
    #[ 1,0 ],
    #[ 0,1 ],
    #[ 1,1 ],
    #[ 0,2 ],
    #[ 1,2 ]
    #])
#elements = np.array([
    #[ node a, node b, node c, property set number]
    #[ 0, 1, 3, 0 ],
    #[ 0, 3, 2, 1 ],
    #[ 2, 3, 5, 0 ],
    #[ 2, 5, 4, 0 ]
    #])

nodes, elements, groups = readmsh('rectangle.msh')
# material and section properties
properties = np.array([
    #[ k ] 
    [ 2e1 ],
    [ 4e1 ]
    ])

loads = 10*np.ones((len(nodes),2))
for i in np.arange(len(loads)):
    loads[i,0]=i
    if nodes[i,0]**2+nodes[i,1]**2>=1:
        loads[i,1]=0
# contour conditions
#loads = np.array([
    ##[ node number, fz]
    #[ 200, 1 ],
    #[ 201, 3 ]
    #])
nan = float("nan")
fixities = np.zeros((len(groups),2))
fixities[:,0]=groups
    #[ node number, x, y, theta]
    # nan means free
for g in groups:
    loads[g,1]=0


R = solvemembrane(nodes, elements, properties, loads, fixities )

#print("\n== DISPLACEMENTS ==\n")
#Dis = np.split(R, len(nodes) )
#for i in list(range(len(Dis))):
    #print("Node", i, "\t", Dis[i][0][0] )
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
