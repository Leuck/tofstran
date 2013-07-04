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

    # element's stiffness matrix
    def Ke(k, A):
        """Returns a membrane element's stiffness matrix."""
        invA = np.linalg.inv(A)
        K=np.zeros((3,3))
        for i in list(range(3)):
            for j in list(range(3)):
                K[i,j] = k*(invA[i,0]*invA[j,0]+invA[i,1]*invA[j,1])
        return K

    # degrees of freedom per node
    dof = 1
    # number of elements
    noe = len(elements)
    # nodes per element
    nne = len(elements[0])-1

    KGsize = np.size(nodes,0)*dof
    KG = np.zeros((KGsize,KGsize))
    np.disp(KG)

    # list containing all elements' stiffness in global coordinates
    kg = []
    for i in list(range(noe)):
        k = properties[elements[i,3],0]
        coord = np.array([
            [ nodes[elements[i,0],0], nodes[elements[i,0],1], 1 ],
            [ nodes[elements[i,1],0], nodes[elements[i,1],1], 1 ],
            [ nodes[elements[i,2],0], nodes[elements[i,2],1], 1 ]
            ])
        #np.disp(coord)
        kg.append( Ke(k,coord) )

    for i in list(range(noe)):
        np.disp(kg[i])

    # assemble global stiffness matrix
    # for every element
    kint=0
    for i in list(range(noe)):
        #dof = elements[i][0:-1]
        jint=0
        # for every line
        for j in elements[i][0:-1]:
            # and every column
            for k in elements[i][0:-1]:
                KG[j,k]+=kg[i][jint,kint]
                kint += 1
            jint += 1
            kint=0


    #np.disp(KG)

    # assemble loads vector
    F = np.zeros((KGsize,1))
    for i in list(range(len(loads))):
        for j in list(range(dof)):
            F[loads[i,0]*dof + j ] = loads[i, 1+j]

    np.disp(F)

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

    return X
# end of function solvetruss2d()

# vim: ts=4 et sw=4 sts=4 ai
