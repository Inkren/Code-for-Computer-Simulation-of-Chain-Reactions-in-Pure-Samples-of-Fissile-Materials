# -*- coding: utf-8 -*-
"Importing modules"
import numpy as np
import math
import matplotlib.pyplot as plt


"The constants needed"
#for future reference critical radius for uranium 235 is 8.37, plutonium 239 is 6.346
material= input("Which material: Uranium-235 plutonium-239 [ur/pu]:")
if material == "ur":
    pabs=0.089/7.694 #possibility of neutron absorption pg 193 physics of manhattan project
    pfis=1.235/7.694 #possibility of neutron fission
    psca=6.37/7.694  #possibility of neutron scattering
    d=3.596 #mfp of uranium 235 in cm from page 58 of physics of manhattan project
        #Probabilities from https://journals.aps.org/pr/abstract/10.1103/PhysRev.101.1012 page 4 multiplied by 1000
    p0=27
    p1=158
    p2=339
    p3=305
    p4=133
    p5=38
    p6=0
if material == "pu":
    pabs=0.052/7.707
    pfis=1.8/7.707
    psca=5.854/7.707
    d=4.108
    
    p0=0
    p1=110
    p2=130
    p3=560
    p4=110
    p5=60
    p6=50

lim= float(input("What is the radius of the sample in cm:")) #

#Starting values
simtime=0

volume=(4/3)*np.pi*(lim)**3 #need volume to get total nuclei

"Neutron values"
class neu:
    def __init__(self, lif, starttime, pos):
        self.starttime=starttime #Saves creation time
        self.pos=(pos) #saves neutron position
        self.lif=(lif) #saves neutron life or death values

        
"lists"
nel=[]#stores neutrons, time and position
time=[] #Stores total time steps
liveneutrons=[] #Stores number of live neutrons at any given 
fission=[]#Stores number of fissions
    



"Functions"

#vector randomiser
def vector(pos):
    x, y, z =pos 
    thetadegrees=np.random.uniform(0, 360) #Randomiser for azimuthal angle
    phidegrees=np.random.uniform(0, 180) #Randomiser for polar angle
    theta=math.radians(thetadegrees) #Converts azimuthal angles degrees to radians for theta
    phi=math.radians(phidegrees) #Converts polar angles degrees to radians for phi
    
    #Cartesian coordinate conversions below
    dx=d*(np.sin(phi)*np.cos(theta)) 
    dy=d*(np.sin(phi)*np.sin(theta))
    dz=d*(np.cos(phi))
    
    return (x+dx, y+dy, z+dz) #updates current position to old position plus vector

def radius(pos):
    x, y, z = pos
    r=math.sqrt(x*x+y*y+z*z) #Use this to make sure the neutron is in bounds
    return r

def neutroncount(nel):
    return sum(1 for neu in nel if neu.lif=="l") #Sums up current live neutrons in list

def interaction(current):
    p=np.random.uniform(0, 1) #Randomises a number and sees which interaction happens
    if p<=pfis:
        fis(current)
    elif p<=pfis+pabs:
        orp(current)
    else:
        sca(current)



def fis(current):
    current.lif="d" #sets current neutron to dead
    checker=0 
    while checker < len(nel) and nel[checker].starttime <= simtime: #This loop searches through the list until it finds a neutron with a higher start time
        checker +=1
    rando=np.random.uniform(1,1000) #Randomiser to see how many neutrons are added.
    if rando<=p0:
        for _ in range(0):
            newneutron= neu("l", simtime, current.pos)
            nel.insert(checker, newneutron)
            
    if p0<rando<=p0+p1:
        for _ in range(1):
            newneutron= neu("l", simtime, current.pos)
            nel.insert(checker, newneutron)
            
    if p0+p1<rando<=p0+p1+p2:
        for _ in range(2):
            newneutron= neu("l", simtime, current.pos)
            nel.insert(checker, newneutron)
            
    if p0+p1+p2<rando<=p0+p1+p2+p3:
        for _ in range(3):
            newneutron= neu("l", simtime, current.pos)
            nel.insert(checker, newneutron)
            
    if p0+p1+p2+p3<rando<=p0+p1+p2+p3+p4:
        for _ in range(4):
            newneutron= neu("l", simtime, current.pos)
            nel.insert(checker, newneutron)
            
    if p0+p1+p2+p3+p4<rando<=p0+p1+p2+p3+p4+p5:
        for _ in range(5):
            newneutron= neu("l", simtime, current.pos)
            nel.insert(checker, newneutron)
            
    if p0+p1+p2+p3+p4+p5<rando<=p0+p1+p2+p3+p4+p5+p6:
        for _ in range(6):
            newneutron= neu("l", simtime, current.pos)
            nel.insert(checker, newneutron)
                
def orp(current):
    current.lif="d" #Sets the neutron value to dead after it is absorbed

        
def sca(current):
    return #This might be scattering things twice since there is already movement in the main while loop
#fix this to actually do it 10 times
for i in range(10):
    nel.append(neu("l", 0, (0, 0, 0)))#Adds 10 neutrons to the list at the start 

counter=0
while len(nel) > 0 and len(nel) < 500:
    counter+=1
    current=nel.pop(0)
    if current.lif=="d":
        continue #removes dead neutron. Maybe store them for later at some point
    while current.lif=="l":
        simtime=current.starttime+1 #Increases the timestep
        time.append(simtime) #Adds current time to time list
        liveneutrons.append(neutroncount(nel)) #updates number of live neutrons at this time step
        current.pos=vector(current.pos) #Changes current neutrons position
        r=radius(current.pos) #Checks whether the neutron is in the limit of the simulation
        if r>lim:
            orp(current)
        else:
            interaction(current) #Sends the neutrons through the interaction function
            
plt.figure()
plt.plot(time, liveneutrons) #plots a graph of live neutrons over time
plt.xlabel("Time")
plt.ylabel("live nuetrons")
plt.title(f"Live neutrons over time for radius = {lim} cm")
plt.grid()
plt.show()


#look at dunlop again
#record a few different results to find crit radius and find crit mass
#record percentage of neutrons lost to escape vs absorption
#later track nuclei lost and effects on neutrons
#Adjust distances and times between interactions
#Add tamper optionsur
#Adjust mean free paths and time look at pg452 and onward of "An Introduction to Computer Simulation Methods Third Edition" equation 11 63

#Tamper stats on pg 71 of physics of manhattan project






