#!/usr/bin/env python

import numpy as np
from scipy.spatial import ConvexHull

eps= 2.2204e-16

def icosahedron(level):

    den = 5
    pi = np.pi

    C = 1 / np.sqrt(1.25)
    t = (2 * pi / den) * np.reshape([x for x in range(den)], (den, 1))

    u1 = C*np.hstack((np.cos(t), np.sin(t), 0.5 * np.ones((den, 1))))
    u2 = C*np.hstack((np.cos(t + 0.2 * pi), np.sin(t + 0.2 * pi), -0.5 * np.ones((den, 1))))
    u = np.vstack((np.array([[0, 0, 1]]), u1, u2, np.array([[0, 0, -1]])))

    if level>0:

        for lev in range(level):

            fcs= ConvexHull(u).simplices
            N= fcs.shape[0]

            U=np.zeros((3*N,3))

            for k in range(N):

                A= u[fcs[k,0],: ]
                B= u[fcs[k,1],: ]
                C= u[fcs[k,2],: ]

                U[3*k:3*(k+1),: ]= 0.5*np.array([A+B, B+C, A+C])

            norm= np.reshape(np.linalg.norm(U, axis= 1), (U.shape[0],1))

            U= U/norm
            u= np.vstack((u, U))

        ind= np.argsort(u[ :,2])[::-1]
        u= u[ind,: ]

        index= np.where(u[ :,2]==0)[0]
        v = u[index, :]

        ind = np.argsort(v[ :,1])[::-1]
        u[index,: ]= v[ind,: ]

    else:
        norm= np.linalg.norm(u)
        u= u/[norm for x in range(3)]
        fcs= ConvexHull(u).simplices

    print(u)
    print(fcs)
    return (u, fcs)


if __name__=='__main__':
    # for debugging purpose
    icosahedron(1)