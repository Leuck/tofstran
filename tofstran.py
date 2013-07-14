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

# solver for 2d frames
def solveframe2d(nodes, elements, properties, loads, fixities, gravityvector ):
    """
    Finds node displacements, element forces and stresses in a \
2-dimensional linear elastic frame.
All arguments are numpy.array arrays, e.g:
nodes = np.array([[x1, y1], [x2, y2]])
elements = np.array([[0,1]])"""

    # element's local stiffness
    def Kle(E, A, I, L, t=0):
        """Returns a beam element's local stiffness matrix."""
        if (t != 1) and (t != 2) :
            K = np.array([
                [ E*A/L, 0, 0, -E*A/L, 0, 0],
                [ 0, 12*E*I/L**3, 6*E*I/L**2, 0, -12*E*I/L**3, 6*E*I/L**2],
                [ 0, 6*E*I/L**2, 4*E*I/L, 0, -6*E*I/L**2, 2*E*I/L],
                [ -E*A/L, 0, 0, E*A/L, 0, 0],
                [ 0, -12*E*I/L**3, -6*E*I/L**2, 0, 12*E*I/L**3, -6*E*I/L**2],
                [ 0, 6*E*I/L**2, 2*E*I/L, 0, -6*E*I/L**2, 4*E*I/L]
                ])
        elif t==1:
            # element articulated at end A
            # remember to write matrix reduction routine before using this
            K = np.array([
                [ E*A/L, 0, 0, -E*A/L, 0, 0 ],
                [ 0, 3*E*I/L**3, 0, 0, -3*E*I/L**3, 3*E*I/L**2 ],
                [ 0, 0, 0, 0, 0, 0 ],
                [ -E*A/L, 0, 0, E*A/L, 0, 0 ],
                [ 0, -3*E*I/L**3, 0, 0, 3*E*I/L**3, -3*E*I/L**2 ],
                [ 0, 3*E*I/L**2, 0, 0, -3*E*I/L**2, 3*E*I/L ]
                ])
        elif t==2:
            # element articulated at both ends
            # remember to write matrix reduction routine before using this
            K = np.array([
                [ E*A/L, 0, 0, -E*A/L, 0, 0 ],
                [ 0, 0, 0, 0, 0, 0 ],
                [ 0, 0, 0, 0, 0, 0 ],
                [ -E*A/L, 0, 0, E*A/L, 0, 0 ],
                [ 0, 0, 0, 0, 0, 0 ],
                [ 0, 0, 0, 0, 0, 0 ]
                ])

        return K

    # element's local mass (Theres something wrong here)
    def Mle(A, rho, L):
        """Returns a beam element's local mass matrix."""
        M = (A*L*rho/2) * np.array([
            [ 1, 0, 0, 0, 0, 0],
            [ 0, 1, 0, 0, 0, 0],
            [ 0, 0, 0, 0, 0, 0],
            [ 0, 0, 0, 1, 0, 0],
            [ 0, 0, 0, 0, 1, 0],
            [ 0, 0, 0, 0, 0, 0]
            ])

        #M = (A*L*rho/420) * np.array([
            #[ 140, 0, 0, 70, 0, 0],
            #[ 0, 156, 22*L, 0, 54, -13*L],
            #[ 0, 22*L, 4*L**2, 0, 13*L, -3*L**2],
            #[ 70, 0, 0, 140, 0, 0],
            #[ 0, 54, 13*L, 0, 156, -22*L],
            #[ 0, -13*L, -3*L**2, 0, -22*L, 4*L**2]
            #])

        return M

    # rotation matrix
    def rotationMatrix(coord1, coord2 ):
        """Returns the rotation matrix for an element based on its\
 coordinates."""

        l = ((coord2[0] -coord1[0])**2 +(coord2[1] -coord1[1])**2)**(1/2)

        # cosine of orientation angle
        co = (coord2[0] - coord1[0] )/l

        # sine of orientation angle
        si = (coord2[1] - coord1[1] )/l

        RM = np.array([
            [co, si, 0, 0, 0, 0 ],
            [-si, co, 0, 0, 0, 0, ],
            [0, 0, 1, 0, 0, 0, ],
            [0, 0, 0, co, si, 0, ],
            [0, 0, 0, -si, co, 0, ],
            [0, 0, 0, 0, 0, 1 ]
            ])
        return RM

    # element's stiffness in global coordinates
    def Kge(K, coord1, coord2 ):
        """Takes an element's local stiffness matrix and node coordinates
        and determines its stiffness in global coordinates."""

        # rotation matrix
        RM = rotationMatrix( coord1, coord2 )

        # element's stiffness in global coordinates
        Kg = np.dot(np.dot(RM.T, K ), RM )
        return Kg

    # degrees of freedom per node
    dof = 3
    # number of elements
    noe = len(elements)

    KGsize = np.size(nodes,0)*3
    KG = np.zeros((KGsize,KGsize))
    MG = np.zeros((KGsize,KGsize))

    # list containing all elements' stiffness in global coordinates
    kg = []
    ke = []
    for i in range(noe):
        E = properties[elements[i,2], 0]
        A = properties[elements[i,2], 2]
        I = properties[elements[i,2], 3]
        t = properties[elements[i,2], 5]
        coord1 = nodes[elements[i,0] ]
        coord2 = nodes[elements[i,1] ]
        L = ((coord2[0] -coord1[0])**2 +(coord2[1] -coord1[1])**2)**(1/2)
        ke.append(Kle(E,A,I,L,t))
        kg.append(Kge(ke[i], coord1, coord2 ))

    ## list containing all elements' mass in global coordinates
    #mg = []
    #for i in range(noe):
        #rho = properties[elements[i,2], 1]
        #A = properties[elements[i,2], 2]
        #coord1 = nodes[elements[i,0] ]
        #coord2 = nodes[elements[i,1] ]
        #L = ((coord2[0] -coord1[0])**2 +(coord2[1] -coord1[1])**2)**(1/2)
        #mg.append(Kge(Mle(A,rho,L), coord1, coord2 ))

    # assemble global stiffness matrix
    for i in range(noe):
        na = elements[i,0]
        nb = elements[i,1]
        for j in range(3):
            for k in range(3):
                KG[na*3+j, na*3+k] += kg[i][j,k]
                KG[nb*3+j, nb*3+k] += kg[i][j+3, k+3]
                KG[na*3+j, nb*3+k] += kg[i][j, k+3]
                KG[nb*3+j, na*3+k] += kg[i][j+3, k]

    ## assemble global mass matrix
    #for i in range(noe):
        #na = elements[i,0]
        #nb = elements[i,1]
        #for j in range(3):
            #for k in range(3):
                #MG[na*3+j, na*3+k] += mg[i][j,k]
                #MG[nb*3+j, nb*3+k] += mg[i][j+3, k+3]
                #MG[na*3+j, nb*3+k] += mg[i][j, k+3]
                #MG[nb*3+j, na*3+k] += mg[i][j+3, k]

    # assemble loads vector
    F = np.zeros((KGsize,1))
    for i in range(len(loads)):
        for j in range(3):
            F[loads[i,0]*3 + j ] = loads[i, 1+j]


    # add gravity loads
    #AG = np.zeros((KGsize,1))
    #for i in range(len(nodes)):
        #AG[i*3] = gravityvector[0]
        #AG[i*3+1] = gravityvector[1]
        #AG[i*3+2] = gravityvector[2]

    #print(AG)
    #F += np.dot(MG, AG)
    #
    #for i in range(noe):
        #na = elements[i,0]
        #nb = elements[i,1]
        #rho = properties[elements[i,2], 1]
        #A = properties[elements[i,2], 2]
        #coord1 = nodes[elements[i,0] ]
        #coord2 = nodes[elements[i,1] ]
        #L = ((coord2[0] -coord1[0])**2 +(coord2[1] -coord1[1])**2)**(1/2)
        #F[na*3+1]+=gravityvector[1]*L*rho*A/2
        #F[nb*3+1]+=gravityvector[1]*L*rho*A/2

    # applies fixity conditions to global stiffness matrix
    for i in range(len(fixities)):
        for j in range(3):
            if not np.isnan(fixities[i,j+1]):
                for k in range(KGsize):
                    KG[ fixities[i,0]*dof+j, k] = 0
                KG[fixities[i,0]*dof+j, fixities[i,0]*dof+j ] = 1
                F[fixities[i,0]*dof+j ] = fixities[i,j+1]

    # solves system
    X = np.linalg.solve(KG,F)

    # forces on elements (in local coordinates)
    Xge = np.split(X, len(nodes) )
    Fe = []
    for i in range(noe):
        na = elements[i,0]
        nb = elements[i,1]
        coord1 = nodes[na]
        coord2 = nodes[nb]
        Xl = np.dot(rotationMatrix(coord1,coord2), np.vstack((Xge[na], Xge[nb] )))
        Fe.append(np.dot(ke[i], Xl ))

    # normal stresses
    S = []
    for i in range(noe):
        N = abs(Fe[i][0,0])
        Ma = abs(Fe[i][2,0])
        Mb = abs(Fe[i][5,0])
        A = properties[elements[i,2], 2 ]
        I = properties[elements[i,2], 3 ]
        y = properties[elements[i,2], 4 ]/2
        Sxa = N/A + Ma*y/I
        Sxb = N/A + Mb*y/I
        S.append(np.array([Sxa, Sxb ]))

    return [X, Fe, S ]
# end of function solveframe2d()


############################
# solver for 2d trusses
############################

def solvetruss2d(nodes, elements, properties, loads, fixities ):
    """Finds node displacements, element forces and stresses in a
     2-dimensional linear elastic frame."""

    # element's local stiffness
    def Kle(E, A, L):
        """Returns a beam element's local stiffness matrix."""
        K = np.array([
            [ E*A/L, 0, -E*A/L, 0 ],
            [ 0, 0, 0, 0 ],
            [ -E*A/L, 0, E*A/L, 0 ],
            [ 0, 0, 0, 0 ]
            ])
        return K

    # rotation matrix
    def rotationMatrix(coord1, coord2 ):
        """Returns the rotation matrix for an element based on its
        coordinates."""

        l = ((coord2[0] -coord1[0])**2 +(coord2[1] -coord1[1])**2)**(1/2)

        # cosine of orientation angle
        co = (coord2[0] - coord1[0] )/l

        # sine of orientation angle
        si = (coord2[1] - coord1[1] )/l

        RM = np.array([
            [co, si, 0, 0 ],
            [-si, co, 0, 0 ],
            [0, 0, co, si ],
            [0, 0, -si, co ]
            ])
        return RM

    # element's stiffness in global coordinates
    def Kge(K, coord1, coord2 ):
        """Takes an element's local stiffness matrix and node coordinates
        and determines its stiffness in global coordinates."""

        # rotation matrix
        RM = rotationMatrix( coord1, coord2 )

        # element's stiffness in global coordinates
        Kg = np.dot(np.dot(RM.T, K ), RM )
        return Kg

    # degrees of freedom per node
    dof = 2
    # number of elements
    noe = len(elements)

    KGsize = np.size(nodes,0)*dof
    KG = np.zeros((KGsize,KGsize))

    # list containing all elements' stiffness in global coordinates
    kg = []
    for i in range(noe):
        E = properties[elements[i,2], 0]
        A = properties[elements[i,2], 1]
        coord1 = nodes[elements[i,0] ]
        coord2 = nodes[elements[i,1] ]
        L = ((coord2[0] -coord1[0])**2 +(coord2[1] -coord1[1])**2)**(1/2)
        kg.append(Kge(Kle(E,A,L), coord1, coord2 ))

    # assemble global stiffness matrix
    for i in range(noe):
        na = elements[i,0]
        nb = elements[i,1]
        for j in range(dof):
            for k in range(dof):
                KG[na*dof+j, na*dof+k] += kg[i][j,k]
                KG[nb*dof+j, nb*dof+k] += kg[i][j+dof, k+dof]
                KG[na*dof+j, nb*dof+k] += kg[i][j, k+dof]
                KG[nb*dof+j, na*dof+k] += kg[i][j+dof, k]

    # assemble loads vector
    F = np.zeros((KGsize,1))
    for i in range(len(loads)):
        for j in range(dof):
            F[loads[i,0]*dof + j ] = loads[i, 1+j]

    # applies fixity conditions to global stiffness matrix
    for i in range(len(fixities)):
        for j in range(dof):
            if not np.isnan(fixities[i,j+1]):
                for k in range(KGsize):
                    KG[ fixities[i,0]*dof+j, k] = 0
                KG[fixities[i,0]*dof+j, fixities[i,0]*dof+j ] = 1
                F[fixities[i,0]*dof+j ] = fixities[i,j+1]

    # solves system
    X = np.linalg.solve(KG,F)

    # forces on elements (in local coordinates)
    Xge = np.split(X, len(nodes) )
    Fe = []
    for i in range(noe):
        na = elements[i,0]
        nb = elements[i,1]
        coord1 = nodes[na]
        coord2 = nodes[nb]
        Fe.append(np.dot(kg[i], np.vstack((Xge[na], Xge[nb] ))))
        Fe[i] = np.dot( rotationMatrix(coord1, coord2 ), Fe[i] )

    # normal stresses
    S = []
    for i in range(noe):
        N = Fe[i][0,0]
        A = properties[elements[i,2], 1 ]
        Sx = N/A
        S.append(np.array([Sx ]))

    return [X, Fe, S ]
# end of function solvetruss2d()

# vim: ts=4 et sw=4 sts=4 ai
