#! /usr/bin/env python3
import numpy as np
from membranesolver import solvemembrane
# data to define a test frame
# structure geometry
nodes = np.array([
    #[ x, y]
    [ 0,0 ],
    [ 1,0 ],
    [ 0,1 ],
    [ 1,1 ],
    [ 0,2 ],
    [ 1,2 ]
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
    #[ node number, x, y, theta]
    # nan means free
    [ 0, 0 ],
    [ 1, 0 ],
    [ 4, 0 ],
    [ 5, 0 ]
    ])

R = solvemembrane(nodes, elements, properties, loads, fixities )

print("\n== DISPLACEMENTS ==\n")
Dis = np.split(R[0], len(nodes) )
for i in list(range(len(Dis))):
    print("Node", i, "\n", Dis[i] )

#print("\n== ELEMENT FORCES ==\n")
#for i in list(range(len(R[1]))):
    #print("Element", i, "\n", R[1][i] )

#print("\n== NORMAL STRESSES ==\nA\t\t    B\n")
#for i in list(range(len(R[2]))):
    #print("Element", i, "\n", R[2][i] )
    
# vim: ts=4 et sw=4 sts=4 ai
