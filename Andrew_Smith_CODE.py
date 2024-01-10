#Written by Andrew Smith, March-26-2022 through March-28-2022. 

#This program uses the Fourth-Order Runge-Kutta method as described below.
#y(n+1) = y(n) + 1/6 [ k1 + 2k2 + 2k3 + k4 ]
#for the following weighted slope approximations k1, k2, k3, and k4:
#k1 = h * f( x(n)      , y(n)       )
#k2 = h * f( x(n) + h/2, y(n) + k1/2)
#k3 = h * f( x(n) + h/2, y(n) + k2/2)
#k4 = h * f( x(n) + h  , y(n) + k3  )
#where h is the step size (arbitrarily small). 

#Methodology:
# 1) Define all differential equations and weighted slope (k1-k4) functions. 
# 2) Loop over all angles in predefined list.
# 3) For each angle:
#     a) Calculate VELOCITY and PHI by Fourth-Order Runge-Kutta method.
#     b) Using VELOCITY and PHI, find instantaneous change in position.
#     c) Store position values in arrays. 
#     c) Plot position at each time step.

#Notes:
#The following updates were made to the code for PROBLEM 4:
#Angles range from 30 to 45 DEGREES with step of 0.01, these were approximated bounds for the maximum range.
#A further narrowing of this maximum range angle was reduced to 39 to 40 DEGREES. 
#To create these angle values, a for loop is created to map values from 0 to 1000 (number of angles to check) to these angle bounds. 
#Range values are printed directly to the console (no plots) and sifted through for associated maximum angles.

#The following updates were made to the code for PROBLEM 5:
#A runner position and velocity are defined using the ratio of intial velocities.
#Similar to PROBLEM 4, the range of angles is gradually narrowed from 60 to 70 DEGREES to 64 to 65 DEGREES.
#Range offsets are printed to the console again and the angle is found when offset=0.

#------------------------------------------------------------------------------------------------------------------


#Import math functions and plotting functions. 
import math
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt



#Allows conversion of input angles to RADIANS. 
def rads(a) :
    '''Converts a DEGREE measurement to RADIANS.'''
    return ( math.pi * a / 180 )



#Set angle values.
p_inp = [10, 30, 45, 40, 60, 70, 80]  #Input PHI angles (DEGREES), constant.
#p_inp = []                           #Array of angle values to calculate, used for PROBLEMS 4&5. 
#p_inp = [64.93955]                   #Ideal launch angle for PROBLEM 5.
p_num = 200 #Number of angles in array. 
p_min = 60  #Minimum angle. 
p_max = 70  #Maximum angle.
#Linearly map specified number of angles (p_num) in defined range (p_min, p_max).
'''
for i in range(p_num):
    ang = p_min + i * (p_max - p_min) / p_num
    p_inp.append( ang ) #Add angle value to array.
'''



#Initial Condition values. 
a = 2.1   #Define alpha constant. PROBLEMS 3&4: a=1.8; PROBLEM 5: a=2.1. 
x = 0     #Initial X position.
y = 0     #Initial Y position.
v = 1     #Initial velocity. 



#Global constants. 
h = 1e-5          #Step size, dT (lower values are more accurate). 
t_max = 3           #Maximum t value, stops computation. 
n = int(t_max // h) #Number of iterations.



#Create a RUNNER with initial conditions given for position and velocity. Used for PROBLEM 5. 
rx = 0          #Runner horizontal position. 
v0 = 26         #Ball initial velocity (m/s).
v1 = 6          #Runner horizontal velocity (m/s). Value varies in PROBLEM 5. 
v_rat = (v1/v0) #Ratio of velocities. The value of "6" is changed to observe effects on the offset value. 



#Define the differential equations. 
def dvdt (v, p) :
    '''DIFF EQ 1: finds dV/dT given velocity (v) and phi (p). Constant a is predefined.'''
    return ( -math.sin(p) - a * v * v )
def dpdt (v, p) :
    '''DIFF EQ 2: finds dPHI/dT given velocity (v) and phi (p).'''
    return ( 1/v * -math.cos(p) )
def dxdt (v, p) :
    '''DIFF EQ 3: finds dX/dT given velocity (v) and phi (p). Used to find x position.'''
    return ( v * math.cos(p) )
def dydt (v, p) :
    '''DIFF EQ 4: finds dY/dT given velocity (v) and phi (p). Used to find y position.'''
    return ( v * math.sin(p) )



#Find the slope approximations k1, k2, k3, and k4 by Runge-Kutta Fourth Order method. 
#Inputs for each function (f) are either dV/dT or dPHI/dT. 
def k1(f, v, p):
    '''Input diff. eq. function (f), velocity (v), and phi (p) values to get slope approximation k1 by Runge-Kutta method.'''
    return f(v,p)
def k2(f, v, p):
    '''Input diff. eq. function (f), velocity (v), and phi (p) values to get slope approximation k2 by Runge-Kutta method.'''
    return f(v,p)
def k3(f, v, p):
    '''Input diff. eq. function (f), velocity (v), and phi (p) values to get slope approximation k3 by Runge-Kutta method.'''
    return f(v,p)
def k4(f, v, p):
    '''Input diff. eq. function (f), velocity (v), and phi (p) values to get slope approximation k4 by Runge-Kutta method.'''
    return f(v,p)



#Define position array to store position values, used to plot results (X,Y) at a given PHI. 
x_vals = []
y_vals = []



#Create plot characteristics before overall computation.
fig = plt.figure() #Reference name for figure. 
ax = fig.add_axes([0.1, 0.1,  0.6, 0.75] )
plt.xlabel('Range (X)')
plt.ylabel('Height (Y)') #Axis labels.
plt.title('Variable Trajectory with Air Resistance') #Title.



#Perform slope calculations using k1, k2, k3, and k4 for VELOCITY and PHI.
#Using VELOCITY and PHI, X and Y positions are defined at each step increment (dT), plot results.
for i in p_inp:
    '''Loop over input angles to track trajectory for each angle.'''
    #Set current PHI value (RADIANS).
    p  = rads( i )
    #Save initial PHI (p0, DEGREES) to reference individual trajectories. 
    p0 = i 
        
    for j in range(0, n):
        '''Perform calculations by evaluating instantaneous VELOCITY and PHI then finding position (X,Y).'''
    
        #Exit case, break out of loop once particle reaches the ground (Y = 0).
        #Y(0) = 0 so values of j > 0 only are reported. 
        if y <= 0 and j > 0: 
            #Print the range as the X value once Y crosses axis. 
            print("ANG: " + str(p0) + ", RANGE: " + str(x) )
            #Print the difference in runner position and ball position (X values).
            #Find the minimum of this offset difference in console, report the associated angle PHI0. 
            #print("ANG: " + str(p0) + ", OFFSET: " + str(x-rx) )
            break
        
        #Find slope values for VELOCITY.
        k1_v = k1(dvdt, v, p)
        k2_v = k2(dvdt, v, p)
        k3_v = k3(dvdt, v, p)
        k4_v = k4(dvdt, v, p)
        #Calculate VELOCITY weighted slope approximation.  
        k_v = (1/6) * (k1_v + 2*k2_v + 2*k3_v + k4_v) * h
        #Update VELOCITY.
        v += k_v

        #Find slope values for PHI.
        k1_p = k1(dpdt, v, p)
        k2_p = k1(dpdt, v, p)
        k3_p = k1(dpdt, v, p)
        k4_p = k1(dpdt, v, p)
        #Calculate PHI weighted slope approximation. 
        k_p = (1/6) * (k1_p + 2*k2_p + 2*k3_p + k4_p) * h   
        #Update PHI.
        p += k_p
        
        #Convert polar velocity to update Cartesian coordinates, add to list of data points.
        x += dxdt(v, p) * h
        x_vals.append(x)
        y += dydt(v, p) * h
        y_vals.append(y)
        
        #Update RUNNER position.
        rx += v_rat * h
               
    #Plot trajectory from X and Y values calculated.
    #Label using angle, loop for each angle to stack trajectories on the same axes. 
    plt.plot( x_vals, y_vals , label="Trajectory of PHI = " + str(p0) + "deg" )
    
    #Reset position and coordinate arrays before next angle computation. 
    x = 0
    y = 0
    x_vals = []
    y_vals = []
    #Reset VELOCITY value.
    v = 1
    #Reset RUNNER position.
    rx = 0


    
#Display and save the final overlapped plotted trajectories.
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) #Legend. 
plt.show() #Display final plot.

#plt.savefig('Problem_3_Paths.png') #Save image as PNG to file folder.


    
