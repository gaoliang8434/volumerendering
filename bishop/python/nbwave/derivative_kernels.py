#!/usr/bin/python

import os
import sys
import numpy as np
import math


def DerivativeCoefficients( N ):
    denomArray = []
    for j in range(1,N+1):
        row = []
        for i in range(0,N):
            row.append( math.pow( j, 2*i+1 ) )
        denomArray.append(row)
    denominator = np.linalg.det(denomArray)
    basereturn = [0.5]
    for i in range(1,N):
        basereturn.append(0)
    coeffs = []
    for i in range(0,N):
        coeffArray = list( denomArray )
        coeffArray[i] = basereturn
        coeff = np.linalg.det(coeffArray)/denominator
        coeffs.append(coeff)
    return coeffs



def TestFunction( x ):
    value = math.sin(2.24345*x)
    return value

def TestFunctionDerivative( x ):
    value = math.cos(2.24345*x)*2.24345
    return value


for n in range(1,17):
    dc = DerivativeCoefficients(n)
    derivative_string = "df(x)/dx = ( \n"
    for i in range(0,n):
        ss = "\t+ ( f(x+%d dx) - f(x-%d dx) ) * ( %25.21e )" % (i+1, i+1, dc[i])
        derivative_string += ss + "\n"
    derivative_string += ")/dx"
    print derivative_string
    print "Demo"
    Nrange = 20
    dx = 4.0 / float(Nrange)
    x0 = -dx*Nrange*0.5
    dxv = 0.01
    for ix in range(0,Nrange):
        x = ix * dx + x0
        derivative = 0
        for i in range(0,n):
            v = TestFunction( float(x + dxv*(i+1)) ) * dc[i]
            v = v - TestFunction( float(x - dxv*(i+1)) ) * dc[i]
            derivative = derivative + v
        derivative = derivative / dxv
        exact = TestFunctionDerivative( float(x) )
        print str(x) + "  " + str(exact) + "  " + str(derivative) + "   " + str(exact-derivative)
    print "\n\n"
