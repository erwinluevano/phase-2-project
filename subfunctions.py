"""
# By submitting this assignment, I agree to the following:
#  “Aggies do not lie, cheat, or steal, or tolerate those who do”
#  “I have not given or received any unauthorized aid on this assignment”
# 
# Name: 		Erwin Luevano
# UIN:          930002548
# Section:		MEEN 357-505
# Assignment:	
# Date:		    

"""
import numpy as np
from math import *
from scipy.interpolate import interp1d

wheel = {'wheel radius': 0.30, 'mass': 1.0} 

speed_reducer = {'type': 'reverted', 'pinion diameter': 0.04, 'gear diameter': 0.07, 'mass': 1.5}

motor = {'stall torque': 170, 'motor no-load torque': 0, 'motor no-load speed': 3.80, 'mass': 5.0}

chassis = {'mass': 659.0}

science_payload = {'mass': 75.0}

power_subsys = {'mass': 90.0}

omega = np.array([0, pi/3, pi/2])

v = np.array([10, 20, 30])

terrain_angle = np.array([-20, -30, 0, 30, 60,75])

planet = {'Mars': 3.72}

wheel_assembly = {1: wheel, 2: speed_reducer, 3: motor}

rover = {1: chassis, 2: wheel_assembly, 3: science_payload, 4: power_subsys} # dict holds rover parts dict titles

def is_array(x):
    if type(x) is not np.ndarray:
        x = np.array([x])
        return x
    else:
        return x



def get_mass(rover):
    if type(rover) is dict:
        m = rover[1]['mass'] + 6*(rover[2][1]['mass'] + rover[2][2]['mass'] + rover[2][3]['mass']) + rover[3]['mass'] + rover[4]['mass']
        return m
        
    else:
        raise Exception(TypeError)




def F_gravity(terrain_angle, rover, planet):
    terrain_angle = is_array(terrain_angle)
    m = get_mass(rover)
    Fgt = []
    for i in terrain_angle:
        Fgt.append((planet['Mars'])*(sin(i*pi/180))*m)
    Fgt = np.array(Fgt)
    return Fgt
    


def tau_dcmotor(omega, motor):
    tau = []
    omega = is_array(omega)

    if type(omega) is np.ndarray and type(motor) is dict:
        for i in omega:
            if i <= 0:
                tau.append(motor['stall torque'])
            elif i > 0 and i < motor['motor no-load speed']:
                torque = motor['stall torque'] - ((motor['stall torque'] - motor['motor no-load torque'])/(motor['motor no-load speed']))*i
                tau.append(torque)
            elif i >= motor['motor no-load speed']:
                tau.append(0)
        #     tau_dcmotor = motor['stall torque'] - ((motor['stall torque'] - motor['motor no-load torque'])/(motor['motor no-load speed']))*i
        #     tau.append(tau_dcmotor)
        # tau = np.array(tau)
        if len(tau)<2:
            tau = tau[0]
        return tau
    else:
        raise Exception(TypeError)
        


def get_gear_ratio(speed_reducer):
    if type(speed_reducer) is dict:
        if speed_reducer['type'] == 'reverted' or speed_reducer['type'] == 'REVERTED':
            Ng = (speed_reducer['gear diameter']/speed_reducer['pinion diameter'])**2
    else:
        raise Exception(TypeError)
    return Ng

def F_drive(omega, rover):
    """

    Parameters
    ----------
    omega : NUMPY ARRAY 
        array of motor shaft speeds
    rover : DICTIONARY
        Dictionary of what the rover consist of and there constants
    Returns
    -------
    Fd : Numpy array 
        force created by tires 

    """
    omega = is_array(omega)
    if type(omega) is not (dict or str) and type(rover) is dict:
        Ng =get_gear_ratio(rover[2][2]) # returns gear ratio 
        tamu = tau_dcmotor(omega, rover[2][3]) # returns torque at shaft
        radius = rover[2][1]['wheel radius'] # radius of wheels 
        torque_output = Ng* np.array([tamu]) # output torque = gear ratio * input torque 
        Fd = 6*torque_output/radius # force = 6 tires * torque / radius 
    return Fd



def F_rolling(omega, terrain_angle, rover, planet, Crr):
    terrain_angle = is_array(terrain_angle)
    omega = is_array(omega)
    Crr = is_array(Crr)
    # checks to make sure all inputs are correct
    if type(omega) and type(terrain_angle) is (dict or str):
        raise Exception('Did not enter a int or float') 
    if len(omega) != len(terrain_angle):
        raise Exception('list do not match')
    if min(terrain_angle) < (-75) or max(terrain_angle) > 75 :
        raise Exception('Angle is not -75<x<75')
    if type(rover) and type(planet) is not dict:
        raise Exception('Did not enter dictionaries')
    if type(Crr) is not (dict or str) and Crr < 0:
        raise Exception('CRR is not a int or float ')
        
    # calculate rolling resistance 
    m = get_mass(rover) # scalar
    Ng = get_gear_ratio(rover[2][2]) #scalar 
    radius = rover[2][1]['wheel radius'] # scalar
    V = radius * omega/ Ng # list # velocity of rover
    normalF = [m* planet['Mars'] * cos(a*pi/180) for a in terrain_angle] # list # normal force acting on each wheel # changed to cos rahter than sin
    Frrsimple = [normalF[x] * Crr[x] for x in range(len(normalF))] # list # return simple rolling resistant 
    Frr = np.array([erf(40*V[i])* Frrsimple[i] for i in range(len(V))]) # equation for rolling resistance # list
    
    return Frr




def F_net(omega, terrain_angle, rover, planet, Crr):
    
    if type(omega) is not np.ndarray:
        omega = np.array([omega])
    if type(terrain_angle) is not np.ndarray:
        terrain_angle = np.array([terrain_angle])
    if type(Crr) is not np.ndarray:
        Crr = np.array([Crr])
    if type(omega) and type(terrain_angle) is (dict or str):
        raise Exception('Did not enter a int or float') 
    if len(omega) != len(terrain_angle):
        raise Exception('list do not match')
    if min(terrain_angle) < (-75) or max(terrain_angle) > 75 :
        raise Exception('Angle is not -75<x<75')
    if type(rover) and type(planet) is not dict:
        raise Exception('Did not enter dictionaries')
    if type(Crr) is not (dict or str) and Crr < 0:
        raise Exception('CRR is not a int or float ')

  
    Fd = F_drive(omega, rover)
    Fg = F_gravity(terrain_angle, rover, planet)
    Fr = F_rolling(omega, terrain_angle, rover, planet, Crr)
    
    Fn = Fd - Fr - Fg

    if len(Fn) < 2:
        Fn = Fn[0]
    return Fn

def motorW(v, rover):
    
    if type(v) is float or int: # this checks to see if the data type of v is a float or integer
        w_sr = v/rover[2][1]['wheel radius'] # determines the angular velocity of the wheel
        w = w_sr*get_gear_ratio(speed_reducer) # uses the gear ratio to determine the angular velocity of the motor
        return(w)
        
    elif type(v) is np.array(): # this checks to see if the data type of v is an array
        for i in v: # loops through the array of v values
            w_sr = i/rover[2][1]['wheel radius'] # determines the angular velocity of the wheel
            w = w_sr*get_gear_ratio(speed_reducer) # uses the gear ratio to determine the angular velocity of the motor
        return(w)

#print(motorW(v, rover))

def mechpower(v, rover):
    
    if type(v) is float or int: # this checks to see if the data type of v is a float or integer
        w_sr = v/rover[2][1]['wheel radius'] # determines the angular velocity of the wheel
        w = w_sr*get_gear_ratio(speed_reducer) # uses the gear ratio to determine the angular velocity of the motor
        torque = tau_dcmotor(w, motor)
        P = w*torque
        return(P)
        
    elif type(v) is np.array(): # this checks to see if the data type of v is an array
        for i in v: # loops through the array of v values
            w_sr = i/rover[2][1]['wheel radius'] # determines the angular velocity of the wheel
            w = w_sr*get_gear_ratio(speed_reducer) # uses the gear ratio to determine the angular velocity of the motor
            torque = tau_dcmotor(w, motor)
            P = w*torque
        return(P)

print(mechpower(v, rover))


































