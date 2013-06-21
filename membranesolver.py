#! /usr/bin/env python3

# Copyright (C) 2012 Ricardo Frederico Leuck Filho
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import sys

############################
# solver for membranes
############################

def solvemembrane(nodes, elements, properties, loads, fixities ):
    """Finds node displacements in a
     2-dimensional linear elastic membrane."""

    # element's stiffness
    def Kle(k, A):
        """Returns a membrane element's stiffness matrix."""
        invA = np.linalg(A)
        K=np.zeros((3,3))
        for i in list(range(3)):
            for j in list(range(3)):
                K[i,j] = invA[1,i]*invA[1,j]+invA[2,i]*invA[2,j]

        return K

    # degrees of freedom per node
    dof = 1
    # number of elements
    noe = len(elements)

    KGsize = np.size(nodes,0)*dof
    KG = np.zeros((KGsize,KGsize))

    # list containing all elements' stiffness in global coordinates
    kg = []
    for i in list(range(noe)):
        k = properties[elements[i,2], 0]
        NC = np.ones((3,3))
        coord[0] = nodes[elements[i,0] ]
        coord[1] = nodes[elements[i,1] ]
        coord[2] = nodes[elements[i,2] ]
        L = ((coord2[0] -coord1[0])**2 +(coord2[1] -coord1[1])**2)**(1/2)
        kg.append(Kge(Kle(E,A,L), coord1, coord2 ))

    # assemble global stiffness matrix
    for i in list(range(noe)):
        na = elements[i,0]
        nb = elements[i,1]
        nc = elements[i,2]
        for j in list(range(dof)):
            for k in list(range(dof)):
                KG[na*dof+j, na*dof+k] += kg[i][j,k]
                KG[nb*dof+j, nb*dof+k] += kg[i][j+dof, k+dof]
                KG[na*dof+j, nb*dof+k] += kg[i][j, k+dof]
                KG[nb*dof+j, na*dof+k] += kg[i][j+dof, k]

    # assemble loads vector
    F = np.zeros((KGsize,1))
    for i in list(range(len(loads))):
        for j in list(range(dof)):
            F[loads[i,0]*dof + j ] = loads[i, 1+j]

    # applies fixity conditions to global stiffness matrix
    bignumber = 8**(sys.float_info.max_10_exp/2)
    for i in list(range(len(fixities))):
        for j in list(range(dof)):
            if fixities[i,j+1] == 0:
                KG[fixities[i,0]*dof+j, fixities[i,0]*dof+j ] = float("inf")
            elif not np.isnan(fixities[i,j+1]):
                KG[fixities[i,0]*dof+j, fixities[i,0]*dof+j ] = \
                        KG[fixities[i,0]*dof+j, fixities[i,0]*dof+j ] * bignumber
                F[fixities[i,0]*dof+j ] = fixities[i,1+j]* \
                        KG[fixities[i,0]*dof+j, fixities[i,0]*dof+j ]

    # solves system
    X = np.linalg.solve(KG,F)

    # forces on elements (in local coordinates)
    Xge = np.split(X, len(nodes) )
    Fe = []
    for i in list(range(noe)):
        na = elements[i,0]
        nb = elements[i,1]
        coord1 = nodes[na]
        coord2 = nodes[nb]
        Fe.append(np.dot(kg[i], np.vstack((Xge[na], Xge[nb] ))))
        Fe[i] = np.dot( rotationMatrix(coord1, coord2 ), Fe[i] )

    # normal stresses
    S = []
    for i in list(range(noe)):
        N = Fe[i][0,0]
        A = properties[elements[i,2], 1 ]
        Sx = N/A
        S.append(np.array([Sx ]))

    return [X, Fe, S ]
# end of function solvetruss2d()

# vim: ts=4 et sw=4 sts=4 ai