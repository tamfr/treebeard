# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 14:55:59 2014

@author: Mott
"""

#import os
import string
from numpy import matrix, zeros, savetxt

class IterResistors(type):
    def __iter__(cls):
        return iter(cls._resistors)

class Resistors(object):
   __metaclass__ = IterResistors
   _resistors = []
   
   def __init__(self, name, part, length, diameter, delta_z, predecessor = None, velocity = None):
        self._resistors.append(self)
        self.name = name
        self.part = part
        self.length = length
        self.diameter = diameter
        self.delta_z = delta_z
        self.predecessor = str(predecessor)
        self.velocity = velocity
        


viscosity = 0.000798 # Pa*s
density = 1000 # kg/m^3
gravity = 9.81 # m/s^2
k_T = 1
k_c = 0
P_atm = 97421.6 # Pa
     

def system(lines):
    # Source
    Resistors('AB','syringe',.1,.267,.1, None,None)
    
    # Line
    for line in range(1,lines+1):
        n = 1
        Resistors(str(line)+string.uppercase[n:n+2],'barb',.001,.00107,.1,string.uppercase[n-1:n+1]); n += 1
        Resistors(str(line)+string.uppercase[n:n+2],'tube',.7985,.0015875,.1,str(line)+string.uppercase[n-1:n+1]); n+=1
        Resistors(str(line)+string.uppercase[n:n+2],'throttle',.0254,.00005,.1,str(line)+string.uppercase[n-1:n+1]); n+=1
        Resistors(str(line)+string.uppercase[n:n+2],'tube',7.1865,.0015875,.1,str(line)+string.uppercase[n-1:n+1]); n+=1
        Resistors(str(line)+string.uppercase[n:n+2],'bioreactor',.07,.0122,.1,str(line)+string.uppercase[n-1:n+1]); n+=1
        Resistors(str(line)+string.uppercase[n:n+2],'tube',3,.0015875,.1,str(line)+string.uppercase[n-1:n+1]); n+=1
    

continuity_equations = {}
n = 1
for resistor in Resistors:
    equation = {}          
    for next_resistor in Resistors:
        if resistor.name == next_resistor.predecessor:
            equation["V_"+ resistor.name] = resistor.diameter**2            
            equation["V_"+ next_resistor.name] = -next_resistor.diameter**2
            continuity_equations[str(n)] = equation
    n += 1

def fillMatrix():
    length = len(Resistors._resistors)
    A = matrix(zeros((2*length,2*length)))
    B = matrix(zeros((2*length,2*length)))
    i = 0; k = 0 
    for resistor in Resistors:                
        j = 0; n = 0        
        for next_resistor in Resistors:
            if resistor.name == next_resistor.predecessor and resistor.velocity == None:            
                A[i,k] = resistor.diameter**2 # Coefficient for the inlet. 
                A[i,j] = -next_resistor.diameter**2 # Coefficient for the outlet(s).
                n = 1
            #if resistor.name == next_resistor.predecessor and resistor.velocity != None:
            #    j -= 1
            j += 1
        if n > 0:
            i += 1 # Iterate i only if n = 1 which signifies coefficients were made for an equation.
        k += 1 # Continue to iterate k regardless of whether an equation was made.
    k = 0
    last_in_line = None
    for resistor in Resistors:                
        j = 0 # Reset j
        n = 0 # Reset n, which gets set to 1 if an equation is made      
        for next_resistor in Resistors:                        
            if resistor.name == next_resistor.predecessor:            
                n += 1                
                last_in_line = next_resistor.predecessor # Tracks the last resistor in a line.                        
                A[i,k] = -32*viscosity*resistor.length/resistor.diameter**2 # Major head loss              
                print i,j                
                A[i+j-(j-1)*(n==1),k+length] = 1 # Function for rows in the A index allows for multiple channels to have same pressure in.
                A[i+j-(j-1)*(n==1),j*(n>1)+(k+1)*(n==1)+length] = -1
                if resistor.diameter < next_resistor.diameter:                 
                    B[i,k] = -density*(1-(resistor.diameter/next_resistor.diameter)**2)**2/2 # Minor headless for sudden expansion.
                elif resistor.diameter > next_resistor.diameter:                 
                    B[i,j] = -density*k_c/2 # Minor headless for sudden contraction.
            j += 1
        if resistor.predecessor == last_in_line: # Exit flows
            A[i,k] = -32*viscosity*resistor.length/resistor.diameter**2 # Assigns the equation coefficient for major head loss for the last resistor in a line.           
            B[i,k] = -density/2 # Minor headless for sudden expansion into large space.                  
            last_in_line = None
            n = 1
        if n > 0:
            i += 1 # Increment i if an equation was made.
        k += 1
    A.resize((i,2*j)); B.resize((i,2*j))
    savetxt("matrix_A.txt", A, fmt='%12.6g', delimiter='\t') # Creates matrix.tex file and enables write mode.
    savetxt("matrix_B.txt", B, fmt='%12.6g', delimiter='\t') # Creates matrix.tex file and enables write mode.
    return A
        