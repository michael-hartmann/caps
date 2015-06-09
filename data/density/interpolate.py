#!/usr/bin/python

from __future__ import division
from math import pi

hbarc = 3.161526510740123e-26
kb    = 1.3806488e-23

class Interpolate:
    def __init__(self, filename, nearest=False):
        # L/R, R/(L+R), T, F, S, lmax, nmax, time
        d = []
        fh = open(filename, "r")
        LbyR_last = None
        for line in fh:
            line = line.strip()
            if line == "" or line[0] == '#':
                continue
            LbyR, RbyScriptL, T, F, S, FbyF0, lmax,nmax,time = map(float, line.split(","))

            if LbyR_last != LbyR:
                d.append([])

            d[-1].append((LbyR, T, S, F, FbyF0))
            LbyR_last = LbyR
    
        self.nearest = nearest
        self.d = d            


    def interpolateS(self, LbyR, T):
        return self.interpolate(LbyR, T, 2)


    def interpolate(self, LbyR, T, index):
        d = self.d
        index_x = None
        index_y = None

        left = 0
        right = len(d)-1
        while True:
            mid = left+int((right-left)/2)
            if d[left][0][0] <= LbyR and d[mid][0][0] >= LbyR:
                right = mid
            else:
                left = mid

            if (right-left) == 1:
                index_x = left
                break

        left = 0
        right = len(d[index_x])-1
        while True:
            mid = left+int((right-left)/2)
            if d[index_x][left][1] <= T and d[index_x][mid][1] >= T:
                right = mid
            else:
                left = mid

            if (right-left) == 1:
                index_y = left
                break

        if self.nearest:
            Q12 = (d[index_x+0][index_y+1][0]-LbyR)**2 + (d[index_x+0][index_y+1][1]-T)**2
            Q22 = (d[index_x+1][index_y+1][0]-LbyR)**2 + (d[index_x+1][index_y+1][1]-T)**2
            Q11 = (d[index_x+0][index_y+0][0]-LbyR)**2 + (d[index_x+0][index_y+0][1]-T)**2
            Q12 = (d[index_x+1][index_y+0][0]-LbyR)**2 + (d[index_x+1][index_y+0][1]-T)**2

            if Q12 < Q22 and Q12 < Q11 and Q12 < Q11:
                return LbyR, T, d[index_x+0][index_y+1][index]
            elif Q22 < Q11 and Q22 < Q12:
                return LbyR, T, d[index_x+1][index_y+1][index]
            elif Q11 < Q12:
                return LbyR, T, d[index_x+0][index_y+0][index]
            else:
                return LbyR, T, d[index_x+1][index_y+0][index]

        Q12 = d[index_x+0][index_y+1]
        Q22 = d[index_x+1][index_y+1]
        R2 = ( LbyR, Q12[1], Q12[index]+ (Q22[index]-Q12[index])/(Q22[0]-Q12[0])*(LbyR-Q12[0]) )

        Q11 = d[index_x+0][index_y+0]
        Q21 = d[index_x+1][index_y+0]
        R1 = ( LbyR, Q21[1], Q21[index]+ (Q11[index]-Q21[index])/(Q11[0]-Q21[0])*(LbyR-Q21[0]) )

        S = R1[2] + (R2[2]-R1[2])/(R2[1]-R1[1])*(T-R1[1])

        return [LbyR, T, S]


    def interpolate_SI(self, LbyR, R, T_SI, index):
        T = 2*pi*kb*R*(LbyR+1)/hbarc*T_SI
        p = self.interpolate(LbyR, T, index)

        S_SI = 2*pi*kb*p[2]

        return [LbyR, T_SI, S_SI]
