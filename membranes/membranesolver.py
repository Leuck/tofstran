#! /usr/bin/env python3

# Copyright (C) 2013 Ricardo Frederico Leuck Filho
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

def readmsh(filename):
    """Get node, element and node group lists."""
    f = open(filename,'r')
    lines = f.readlines()
    for line in np.arange(len(lines)):
        # get mesh groups
        if lines[line]=='$PhysicalNames\n':
            numberOfGroups = int(lines[line+1])
            groups = np.empty([numberOfGroups, 3])
            for i in np.arange(numberOfGroups):
                groups[i] = np.fromstring(lines[line+2+i],sep=' ')
            #np.disp(groups)

        #get node coordinate list
        elif lines[line]=='$Nodes\n':
            numberOfNodes = int(lines[line+1])
            nodes = np.empty([numberOfNodes, 4])
            for i in np.arange(numberOfNodes):
                nodes[i] = np.fromstring(lines[line+2+i],sep=' ')
        # get elements
        elif lines[line]=='$Elements\n':
            numberOfElements = int(lines[line+1])
            elements = []
            lineElements = []
            for i in np.arange(numberOfElements):
                elementType = np.fromstring(lines[line+2+i],sep=' ')[1]
                elementGroup = np.fromstring(lines[line+2+i],sep=' ')[3]
                if elementType==1: # if line element
                    lineElements.append(np.fromstring(lines[line+2+i],sep=' '))
                elif elementType==2: # if triangular element
                    elements.append(np.fromstring(lines[line+2+i],sep=' '))
    lineelements = np.zeros((len(lineElements),2))
    for i in np.arange(len(lineElements)):
        lineelements[i]=lineElements[i][-2:]-1
    lineelements = lineelements.flatten()
    lineelements = np.sort(lineelements)
                
    triaelements = np.zeros((len(elements),4))
    for i in np.arange(len(elements)):
        triaelements[i,0:3]=elements[i][-3:]-1

    nodes = nodes[:,1:4]

    return nodes,triaelements,lineelements

############################
# solver for membranes
############################

def solvemembrane(nodes, elements, properties, loads, fixities ):
    """Finds node displacements in a
     2-dimensional linear elastic membrane."""

    # element's stiffness matrix
    def Ke(k, A):
        """Returns a membrane element's stiffness matrix."""
        #np.disp(A)
        invA = np.linalg.inv(A)
        area = np.linalg.det(A)*0.5
        K=np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                K[i,j] = area*k*(invA[0,i]*invA[0,j]+invA[1,i]*invA[1,j])
        return K

    # degrees of freedom per node
    dof = 1
    # number of elements
    noe = len(elements)
    # nodes per element
    nne = len(elements[0])-1

    KGsize = np.size(nodes,0)*dof
    KG = np.zeros((KGsize,KGsize))

    # list containing all elements' stiffness in global coordinates
    kg = []
    for i in range(noe):
        k = properties[elements[i,3],0]
        coord = np.array([
            [ nodes[elements[i,0],0], nodes[elements[i,0],1], 1 ],
            [ nodes[elements[i,1],0], nodes[elements[i,1],1], 1 ],
            [ nodes[elements[i,2],0], nodes[elements[i,2],1], 1 ]
            ])
        #np.disp(coord)
        kg.append( Ke(k,coord) )

    #for i in range(noe):
        #np.disp(kg[i])

    # assemble global stiffness matrix
    # for every element
    kint=0
    for i in range(noe):
        jint=0
        # for every line
        for j in elements[i][0:-1]:
            # and every column
            for k in elements[i][0:-1]:
                KG[j,k]+=kg[i][jint,kint]
                kint += 1
            jint += 1
            kint=0

    # assemble loads vector
    F = np.zeros((KGsize,1))
    for i in range(len(loads)):
        for j in range(dof):
            F[loads[i,0]*dof + j ] = loads[i, 1+j]

    bignumber = 1000*KG.max()
    # applies fixity conditions to global stiffness matrix
    for i in range(len(fixities)):
        for j in range(dof):
            if fixities[i,j+1] == 0:
                KG[fixities[i,0]*dof+j, fixities[i,0]*dof+j ] = float("inf")
            elif not np.isnan(fixities[i,j+1]):
                KG[fixities[i,0]*dof+j, fixities[i,0]*dof+j ] = \
                        KG[fixities[i,0]*dof+j, fixities[i,0]*dof+j ] * bignumber
                F[fixities[i,0]*dof+j ] = fixities[i,1+j]* \
                        KG[fixities[i,0]*dof+j, fixities[i,0]*dof+j ]

    print("\n== GLOBAL STIFFNESS MATRIX ==\n")
    np.disp(KG)
    #print("\n== GLOBAL LOADS VECTOR ==\n")
    #np.disp(F)
    # solves system
    X = np.linalg.solve(KG,F)

    return X
# end of function solvemembrane()

# vim: ts=4 et sw=4 sts=4 ai
