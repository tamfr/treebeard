# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 14:19:08 2014

@author: Mott

Numerical root finder using Newton's method for systems of equations with form Ax+Bx^2+C.
"""
# Imports necessary math modules for linear algebra, sig figs, division, and importing from Excel.

from __future__ import division
from math import floor, log10, pi
from numpy import matrix, zeros, ones, multiply, divide, arange, tile
from pysnatch import xlsnatch
import os
mypath = os.path.dirname(__file__) # Makes available the absolute path to the directory of this file such that if it is imported and gets executed outside its original directory, any subsequently imported files can be found.

def pipe():
    """
    Solves systems of equations to determine mass flow rates and pressures in a fluidics network.
    """ 
    total_unknowns = 222 # Total number of unkown variables.  
    file_name = 'coefficients-20.xlsx'
    variable_names = xlsnatch(os.path.join(mypath,file_name),1,2,2,1,total_unknowns+2,"list") # Import variable names from Excel as list. 
    num_pressure_unknowns = 0 # Initiate variable.
    for name in variable_names[0][:total_unknowns]:
        num_pressure_unknowns += (name[:1].encode('ascii','ignore')=="P") # Determine how many pressure unknowns there are.
    
    num_channels = 0
    for name in variable_names[0][:total_unknowns]:
        num_channels += (name[:3].encode('ascii','ignore')=="P_B") 
       
    Diameters = xlsnatch(os.path.join(mypath,file_name),0,50,2,1,total_unknowns-num_pressure_unknowns+1,"matrix") # Snatch diameters [m] for later calculations.
    P_exit = xlsnatch(os.path.join(mypath,file_name),0,7,2,1,1,"list")[0][0] # Snatch the exit pressure.
    vel_given = xlsnatch(os.path.join(mypath,file_name),0,32,2,1,1,"list")[0][0] # Snatch the given velocity.
    rho = 1000 # Density of water [kg/m^3].
    
    # Set coefficients for equation Ax+Bx^2+C.
    A = matrix(zeros((total_unknowns,total_unknowns))) # Initiate matrix.
    B = matrix(zeros((total_unknowns,total_unknowns))) # Initiate matrix.
    C = matrix(zeros((len(A),1))) # Initiate matrix.
    A = xlsnatch(os.path.join(mypath,file_name),1,3,2,total_unknowns,total_unknowns,"matrix") # Imports data from Excel file by using absolute path to file to propagate matrix of coefficients.
    B = xlsnatch(os.path.join(mypath,file_name),2,3,2,total_unknowns,total_unknowns,"matrix") # Imports data from Excel file by using absolute path to file to propagate matrix of coefficients.
    C = xlsnatch(os.path.join(mypath,file_name),3,3,2,total_unknowns,1,"matrix") # Imports data from Excel file by using absolute path to file to propagate matrix of coefficients.
    
    f=lambda x: A*x+B*multiply(x,x)+C # Define function f to solve for roots.
    
    J = matrix(zeros((A.shape[0],A.shape[1]))) # Initiate Jacobian based on number of equations as defined by matrix A.
    
    # Define parameters to run loop
    n = 0 # Define iteration counter.
    max_iterations = 600 # Define max number of iterations.
    tol = 1e-6 # Define tolerance upon which to interate.
    
    # Initial guesses
    #vel_guess = .003 # Good guest found to be .0000003 to yield all positive roots.
    pressure_guess = 100000
    
    X = matrix(ones((A.shape[0],1))*pressure_guess)
    #X[2:num_pressure_unknowns,0] = multiply(matrix(tile(arange(1,pressure_guess/P_exit+(pressure_guess/P_exit-1)/((num_pressure_unknowns-2)/num_channels),(pressure_guess/P_exit-1)/((num_pressure_unknowns-2)/num_channels-1)),num_channels)).T[::-1],P_exit)
    X[num_pressure_unknowns:X.shape[0],0] = multiply(divide(1,multiply(Diameters[0,:-1].T,Diameters[0,:-1].T)),vel_given*Diameters[0,-1]**2)
    #X[num_pressure_unknowns:X.shape[0],0] = vel_guess

    # Newton's method for finding roots
    while any(abs(v) > tol for v in f(X)):
    
        # propagate jacobian (A+Bx)
        for i in range(0,J.shape[0]):
            for j in range(0,J.shape[1]):
                J[i,j] = A[i,j] + B[i,j]*X[j]
            
        X = X - J.I*f(X)
        
        #X = matrix([(abs(X[x,0]) if X[x,0] < 0 else X[x,0]) for x in range(0,len(X))]).T   
        
        if n < max_iterations:
            n = n + 1
            print n
        else:
            tol = 999999999999999999
    
    pressure_precision = 10 # Define number of sig figs.  
    velocity_precision = 5 # Define number of sig figs. 
    ans = []    
    ans[0:num_pressure_unknowns] = [round(x,-int(floor(log10(abs(x)/10**pressure_precision)))) for x in X[0:num_pressure_unknowns,0]] # Round numbers to certain number of sig figs based on precision      
    ans[num_pressure_unknowns:X.shape[0]] = [round(x,-int(floor(log10(abs(x)/10**velocity_precision)))) for x in X[num_pressure_unknowns:X.shape[0],0]] # Round numbers to certain number of sig figs based on precision
    ans.append(round(vel_given,-int(floor(log10(abs(vel_given)/10**velocity_precision))))) # Append the given velocity
    ans.append(P_exit) # Append the given pressure
    
    m_dot = [round(x,-int(floor(log10(abs(x)/10**velocity_precision)))) for x in multiply(multiply(Diameters[0,:-1].T,Diameters[0,:-1].T),X[num_pressure_unknowns:,0])*pi*rho/4] # Calculate mass flow rates in each branch.
    m_dot.append(round(ans[len(ans)-1]*Diameters[0,len(C)-num_pressure_unknowns]**2*pi*rho/4,-int(floor(log10(abs(ans[len(ans)-1]*Diameters[0,len(C)-num_pressure_unknowns]**2*pi*rho/4)/10**velocity_precision))))) # Append the calculated mass flow rate from given velocity
    
    
    variables = {}
    for i in range(0,len(ans)):
        variables[variable_names[0][i].encode('ascii','ignore')] = ans[i]
        if len(ans)-1 > i >= num_pressure_unknowns:
            variables['m_'+variable_names[0][i].encode('ascii','ignore')[2:]] = m_dot[i-num_pressure_unknowns]

    return variables