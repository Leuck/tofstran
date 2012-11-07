#! /usr/bin/env python3
import numpy as np
from tofstran import solvetruss2d
# data to define a test frame
# structure geometry
nodes = np.array([
    #[ x, y]
    [ 0,0 ],
    [ 1,0 ],
    [ 2,0 ],
    [ .5,1 ],
    [ 1.5,1 ]
    ])
elements = np.array([
    #[ node a, node b, property set number]
    [ 0, 1, 0 ],
    [ 1, 2, 0 ],
    [ 0, 3, 0 ],
    [ 1, 3, 0 ],
    [ 1, 4, 0 ],
    [ 2, 4, 0 ],
    [ 3, 4, 0 ]
    ])

# material and section properties
properties = np.array([
    #[ E, A, I, h ]
    [ 2e11, 2.5e-3,  2.0833e-6, .1 ],
    [ 2e11, 6.25e-4, 3.2552e-8, .025 ]
    # rectangular sections, 25x100 and 25x25 mm
    ])

# contour conditions
loads = np.array([
    #[ node number, Fx, Fy, Mz]
    [ 1, 0, -1700, 0 ],
    [ 4, -100, 0, 0 ]
    ])
nan = float("nan")
fixities = np.array([
    #[ node number, x, y, rotation]
    # nan means free
    [ 0, 0, 0, nan],
    [ 2, nan, 0, nan]
    ])

R = solvetruss2d(nodes, elements, properties, loads, fixities )

print("\n== DISPLACEMENTS ==\n")
Dis = np.split(R[0], len(nodes) )
for i in list(range(len(Dis))):
    print("Node", i, "\n", Dis[i] )

print("\n== ELEMENT FORCES ==\n")
for i in list(range(len(R[1]))):
    print("Element", i, "\n", R[1][i][2][0] )

print("\n== NORMAL STRESSES ==\nA\t\t    B\n")
for i in list(range(len(R[2]))):
    print("Element", i, "\n", R[2][i] )
    
# vim: ts=4 et sw=4 sts=4 ai
