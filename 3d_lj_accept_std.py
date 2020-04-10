import matplotlib.pyplot as plt
import numpy as np
from numpy import random
import copy
from matplotlib.animation import FuncAnimation
import matplotlib
import os
from mpl_toolkits.mplot3d import Axes3D
import copy
import math
import pandas as pd
import datetime

tstart = datetime.datetime.now()

class Atom:
    x = y = z = 0.0
    mass = 0.0 
    radius = 0.0
    energy = 0.0
class Configuration:
    def __init__(self,na):
        self.natoms = na
        self.COM = [0.0,0.0,0.0]
        self.atom = []
        for i in range(na):
            self.atom.append(Atom())
    
    def CalcCOM(self):
        M = 0.0
        sumx = sumy = sumz = 0.0
        for i in range(0, self.natoms):
            m = self.atom[i].mass
            sumx += self.atom[i].x * m
            sumy += self.atom[i].y * m
            sumz += self.atom[i].z * m
            M += m

        Cx = sumx/M
        Cy = sumy/M
        Cz = sumz/M
        self.COM[0]=Cx
        self.COM[1]=Cy
        self.COM[2]=Cz
    
    def RadGyr(self):
        sumgyr = 0.0
        M = 0.0
        self.CalcCOM()
        comx = self.COM[0]
        comy = self.COM[1]
        comz = self.COM[2]
        for i in range(0, self.natoms):
            M += self.atom[i].mass
            mc = self.atom[i].mass
            rx = self.atom[i].x
            ry = self.atom[i].y
            rz = self.atom[i].z
            sgx = (rx - comx)**2.0
            sgy = (ry - comy)**2.0
            sgz = (rz - comz)**2.0

            sumgyr += mc*(sgx+sgy+sgz)

        Rsq = sumgyr / M
        R = math.sqrt(Rsq)
        return R
N = 10  #number of partical
L = 15 #length of box
d = L/(N)  #diameter of hard ball
steps = 50000  #total steps
sigma=1.5*d   #sigma in lj potential
epsilon=2*d  #epsilon in lj potential
temp=300 #temperature 
kb=1.38e-23*6.022e23/(1000)  #bolzmann constant
cutoff=L/2-1  #cutoff of pair-wise interaction
delta=0.5  #maximum move of x
mass = 12.0 #mass of partical
center = [L/2,L/2,L/2] #center of box
#check if inbox
def inbox(x,y,z):     
    if x > L-d/2:
        x -= L-d
    elif x < d/2:
        x += L - d

    if y > L-d/2:
        y -= L-d
    elif y < d/2:
        y += L - d
  
    if z > L-d/2:
        z -= L-d
    elif z < d/2:
        z += L - d   
    return x,y,z
#lj potential
def lj_potential(r):
    energy = 4*epsilon*((sigma/r)**12-(sigma/r)**6)
    return energy
#pair-wise interaction calculation
def energy(r,x,y,z,n,compute_all):
    if compute_all:
        r[n][0],r[n][1],r[n][2]=x,y,z
    x_left = list(r[:,0]-L)
    x_right = list(r[:,0]+L)
    y_left = list(r[:,1]-L)
    y_right = list(r[:,1]+L)
    z_left = list(r[:,2]-L)
    z_right = list(r[:,2]+L)
    all_box = np.array([(x_left*3+list(r[:,0])*3+x_right*3)*3, \
                        ((y_left+list(r[:,1])+y_right)*3)*3, \
                        z_right*9+list(r[:,2])*9+z_left*9])
    all_box = all_box.transpose((1,0))
    #print(np.shape(all_box))
    i = 13*N
    distance = []
    if compute_all: 
        while i < 14*N:
            distance += [np.sqrt((all_box[i][0]-all_box[n][0])**2+ \
                                 (all_box[i][1]-all_box[n][1])**2+ \
                                 (all_box[i][2]-all_box[n][2])**2) for n in \
                         list(np.linspace(0,13*N-1,13*N,dtype=int))+ \
                         list(np.linspace(i+1,27*N-1,27*N-i-1,dtype=int))]
            i += 1
      #  print(distance)
        energy = [lj_potential(d) for d in distance if d < cutoff and d > 0]
        #print(energy)
        energy_sum = sum(energy)
    else:
        #print(all_box[13*N**2+n])       
        r0 = np.array([[x for d in range(27*N)],[y for d in range(27*N)],\
                       [z for d in range(27*N)]]) 
        a = r0.transpose((1,0)) - all_box  #compare the updated particle coordinate with others
        distance = [np.sqrt(a[i][0]**2+a[i][1]**2+a[i][2]**2) for i in range(27*N)]  #calculate distance    
        #print([d for d in distance if d < cutoff and d > 0])    
        Energy = [lj_potential(d) for d in distance if d < cutoff and d > 0]
        #print(Energy)
        energy_sum = sum(Energy)
        #print(energy_sum)
    return energy_sum
#metropolis move criterion
def move(r, x,y,z,n,e, diff_e, temp,accepted_step):                         
    '''x,y are new move coordinate; n represents all coordinate, 
       i is which particle, N is how many particles in total. this
       function will compute the total distance after moving one particle '''
   # n[i] = [x,y,z]      #update coordinate             
   # n0 = np.array([[x for d in range(N)],[y for d in range(N)],[z for d in range(N)]]) 
   # a = n0.transpose((1,0)) - n  #compare the updated particle coordinate with others
   # distance = [np.sqrt(a[i][0]**2+a[i][1]**2+a[i][2]**2) for i in range(N)]  #calculate distance

    if diff_e < 0.0:
    #    print(diff_e)
        e = e_new
        #pre_e += diff_e
        r[n][0],r[n][1],r[n][2] = inbox(x_new,y_new,z_new)
        accepted_step+=1
    else:
        rand = random.uniform(0,1)
       # print(diff_e)
        #print(rand)
        #print(math.exp(-diff_e/(kb*temp)))
        if math.exp(-diff_e/(kb*temp)) > rand:
            e = e_new
            accepted_step+=1
           # pre_e += e
            r[n][0],r[n][1],r[n][2] = inbox(x_new,y_new,z_new)
       # else:
            #r = origi
            #e = pre_e
    return r,e,accepted_step
#initialize system
c = Configuration(N) 
#x_init = sorted([i*d for i in range(1,N+1)]*N)
#y_init = [i*d for i in range(1,N+1)]*N
#z_init = [d/2 for i in range(N)]*N
for index in range(N):
    c.atom[index].x = center[0]+ (L-d) * (random.random()-0.5)
    c.atom[index].y = center[1]+ (L-d) * (random.random()-0.5)
    c.atom[index].z = center[2]+ (L-d) * (random.random()-0.5)
    c.atom[index].mass = mass
    c.atom[index].radius = d/2
point = np.array([[c.atom[i].x for i in range(N)],[c.atom[i].y for i in range(N)],
                  [c.atom[i].z for i in range(N)]])
point = point.transpose((1,0))
e = energy(point,0,0,0,0,compute_all=True)
#print(e)
all_energy = np.zeros(steps+1)
all_energy[0]=e
all_frame = np.zeros((steps+1,N,3))
all_frame[0] = point
all_accepted_ratio=[]
accepted_step=0
average_energy= np.zeros(steps+1)
average_energy[0] = e
RG = np.zeros(steps+1)
RG[0] = c.RadGyr()
all_particle_energy = np.zeros((steps,N))
for i in range(steps):
    #pre_e = e
  #  origi = copy.deepcopy(point)
    for m in range(N): 
        x,y,z = point[m][0],point[m][1],point[m][2]
       # e_old = energy(point,x,y,z,m,compute_all=True)
       # print(e_old)
        x_new = x+delta*random.uniform(-1,1)
        y_new = y+delta*random.uniform(-1,1)
        z_new = z+delta*random.uniform(-1,1)
        e_new = energy(point,x_new,y_new,z_new,m,compute_all=True)
       # x_new,y_new,z_new = inbox(x_new,y_new,z_new)
        #e_new = energy(point,x_new,y_new,z_new,m,compute_all=False)
        #print(e_new)
        diff_e = e_new - e
        point, e, accepted_step = move(point,x_new,y_new,z_new,m,e,diff_e,temp,accepted_step)
        c.atom[m].x,c.atom[m].y,c.atom[m].z = point[m][0],point[m][1],point[m][2]
        all_particle_energy[i][m] = energy(point,point[m][0],point[m][1],point[m][2],m,compute_all=False)
        #if point[m][0] == x and point[m][1] == y and point[m][2] == z:
           # accepted_step += 1
         #   all_particle_energy[i][m] = e_new
        #else:
       #     accepted_step += 1
          #  all_particle_energy[i][m] = e 
    RG[i+1] = c.RadGyr()
    all_energy[i+1] = e
    all_frame[i+1] = point
    average_energy[i+1] = all_energy[0:i+2].sum()/(i+2)
tend=datetime.datetime.now()
tdiff=tend-tstart
print(tdiff)
all_accepted_ratio.append(accepted_step/(steps*N))
data={}
RG = pd.DataFrame(RG[::100],index=range(0,steps+1,100))
RG.to_csv('/ufrc/alberto.perezant/liweichang/Tutorial/HW4/noexchg_2/std_rg_40.csv',index=False)
all_frame_rearrange=np.zeros((len(range(0,steps+1,100)),N*3))
for i in range(len(range(0,steps+1,100))):
    for n in range(N):
        for a in range(3):
            all_frame_rearrange[i][n*3+a]=all_frame[i][n][a]
all_frame = pd.DataFrame(all_frame_rearrange,index=range(0,steps+1,100),columns=['x','y','z']*N)
all_frame.to_csv('/ufrc/alberto.perezant/liweichang/Tutorial/HW4/noexchg_2/std_all_frame_40.csv',index=False)
E_tot = pd.DataFrame(all_energy[::100],index=range(0,steps+1,100)) 
E_tot.to_csv('/ufrc/alberto.perezant/liweichang/Tutorial/HW4/noexchg_2/std_E_tot_40.csv',index=False)
E_ave = pd.DataFrame(average_energy[::100],index=range(0,steps+1,100))
E_ave.to_csv('/ufrc/alberto.perezant/liweichang/Tutorial/HW4/noexchg_2/std_E_ave_40.csv',index=False) 
atomic_ave_e = pd.DataFrame(all_particle_energy[::100],index=range(0,steps,100))
atomic_ave_e.to_csv('/ufrc/alberto.perezant/liweichang/Tutorial/HW4/noexchg_2/std_atomic_E_ave_40.csv',index=False) 
print(all_accepted_ratio)
#E_tot = pd.DataFrame(all_energy,index=range(0,steps+1))
#print(E_tot)
#print(atomic_ave_e)
