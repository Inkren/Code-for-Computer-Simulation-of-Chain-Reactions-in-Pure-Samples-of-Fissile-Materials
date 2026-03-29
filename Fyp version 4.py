# -*- coding: utf-8 -*-
"Importing modules"
import numpy as np
import math
import matplotlib.pyplot as plt


"The constants needed"
#for future reference critical radius for uranium 235 is 8.37, plutonium 239 is 6.346
material= input("Which material: Uranium-235, plutonium-239 [ur/pu]:")
if material == "ur":
    mat= "Uranium-235"
    pabs=0.089/7.694 #possibility of neutron absorption pg 193 physics of manhattan project
    pfis=1.235/7.694 #possibility of neutron fission
    psca=6.37/7.694  #possibility of neutron scattering
    mfp=3.596 #mfp of uranium 235 in cm from page 58 of physics of manhattan project
    density=18.71 #grams/cm^3
    #Probabilities of how many neutrons will be released from https://journals.aps.org/pr/abstract/10.1103/PhysRev.101.1012 page 4 multiplied by 1000
    p0=27
    p1=158
    p2=339
    p3=305
    p4=133
    p5=38
    p6=0
if material == "pu":
    mat= "Plutonium-239"
    pabs=0.052/7.706
    pfis=1.8/7.706
    psca=5.854/7.706
    mfp=4.108
    density=15.6 #grams/cm^3
    p0=0
    p1=110
    p2=130
    p3=560
    p4=110
    p5=60
    p6=50

tamper= input("What material do you want for the tamper if any?: No tamper, Aluminium, Depleted Uranium, Berilium Oxide,Tungsten Carbide[no/al/ur/bo/tc]")
if tamper=="al":
    tampermfp=5.595 #tamper mean free path(in cm) found on pg 71 of physics of manhatten project
if tamper=="ur":
    tampermfp=4.342
if tamper=="bo":
    tampermfp= 2.541
if tamper=="tc":
    tampermfp=3.159

lim= float(input("What is the radius of the sample in cm:")) 

#Starting values
simtime=0

velocity=1500000000 #Velocity of a fast neutron past 1MeV from https://www.nuclear-power.com/nuclear-power/reactor-physics/atomic-nuclear-physics/fundamental-particles/neutron/fast-neutrons-high-energy-neutrons/
volume=(4/3)*np.pi*(lim)**3 #need volume to get total nuclei
mass= volume*density

"Neutron values"
class neu:
    def __init__(self, lif, starttime, pos):
        self.starttime=starttime #Saves creation time
        self.pos=(pos) #saves neutron position
        self.lif=(lif) #saves neutron life or death values
        
"Dictionary sorting system"#Gotten from here https://www.geeksforgeeks.org/python/python-sort-python-dictionaries-by-key-or-value/        
#check if a given key exists in the dictionary
def key_exists(dictionary, value:float, Tolerance=0.1):
    return any(abs(key - value) < Tolerance for key in dictionary)

# adds a liveneutron for keys where: startime <= key <= endtime.
# appends both these times if they are not found in the liveneutron
def AddNeutronToSimTime(LiveNeutrons, starttime:int, endtime:int):
    maxNeutron = 0
    timelim = starttime
    #increase the live neutroncount for keys where: startime <= key <= endtime. 
    for key in sorted(LiveNeutrons):
        if starttime <= key <= endtime:
            LiveNeutrons[key] += 1
            maxNeutron = max(LiveNeutrons[key],maxNeutron)
            timelim = key
    
    #check if the starttime exists
    if not key_exists(LiveNeutrons, starttime):
        smaller_keys = [k for k in LiveNeutrons if k < starttime]
        if smaller_keys:
            closest_smaller = max(smaller_keys)
            LiveNeutrons[starttime] = LiveNeutrons[closest_smaller]
        else:
            LiveNeutrons[starttime] = 1
    
    #check if the end time exists
    if not key_exists(LiveNeutrons, endtime):
        smaller_keys = [k for k in LiveNeutrons if k > endtime]
        if smaller_keys:
            closest_bigger = min(smaller_keys)
            LiveNeutrons[endtime] = LiveNeutrons[closest_bigger]
        else:
            LiveNeutrons[endtime] = 1
    
    return maxNeutron >= 500, timelim

#sort from smallest to largest by keyes
def sortLibrary(library, cutoff_time):
    sorted_keys = sorted(library.keys())
    return {k: library[k] for k in sorted_keys if k < cutoff_time}
        
"lists"
nel=[]#stores neutrons, time and position
time={} #Stores total time steps
liveneutrons=[] #Stores number of live neutrons at any given 
fission=[]#Stores number of fissions
    
"Functions"

#calculates distance using formula from pg431 and onward of "An Introduction to Computer Simulation Methods Third Edition" equation 11 63
def distance(mfp):
    random=np.random.uniform(0,1)
    d=-mfp*(np.log(random))
    return d

#vector randomiser
def vector(pos, d):
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
    checker = 0
    while checker < len(nel) and nel[checker].starttime <= simtime:
        checker += 1
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
timelim = 0
counter=0
simdone = False
while len(nel) > 0 and not simdone:
    counter+=1
    current=nel.pop(0)
    if current.lif=="d":
        continue #removes dead neutron. Maybe store them for later at some point
    while current.lif=="l":
        d=distance(mfp)
        current.pos=vector(current.pos, d) #Changes current neutrons position
        sampleradius=radius(current.pos) #Checks whether the neutron is in the limit of the simulation
        if sampleradius>lim: #+ (mfp*0.71):#This is the effective radius. Double check source library.institutpendidikan.ac.id/wp-content/uploads/2024/08/Stacey-Nuclear_Reactor_Physics-2nd.pdf page 47 eqn 3.15
            orp(current)
        else:
            interaction(current) #Sends the neutrons through the interaction function
        oldsimtime=simtime
        simtime=current.starttime+(d/velocity)*10**(9) #Increases the time (Unit is nanoseconds)
        simdone, timelim = AddNeutronToSimTime(time,oldsimtime,simtime)

time = sortLibrary(time,timelim)      
plt.figure()
plt.plot(time.keys(), time.values()) #plots a graph of live neutrons over time
plt.xlim(0, timelim)#cuts off graph after the peak
plt.xlabel("Time(nanoseconds)")
plt.ylabel("live nuetrons")
plt.title(f"Live neutrons over time for radius = {lim} cm in {mat}")
plt.grid()
plt.show()

print("Mass =", {mass},  "grams")