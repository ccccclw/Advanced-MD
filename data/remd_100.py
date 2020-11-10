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
import xarray as xr
import netCDF4 as ncf
import datetime
tstart = datetime.datetime.now()
class Atom:
    x = y = z = 0.0
    mass = 0.0 
    radius = 0.0
    energy = []
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
            
        
N = 10
L = 15
d = L/N
steps = 200
numexchg = 250
sigma=1.5*d  #Angstrom
epsilon=2*d #kJ/mol
#temp=[100.00,246.00,300.00,362.00,436.00,520.00]
temp=[300,400,500,600,700,800,900,1000]
#temp=[300.00, 353.49, 414.11, 482.85, 560.86, 600.00]
#temp=[300.00, 353.49, 414.11, 682.85, 860.86, 10000.00]
kb=1.38e-23*6.022e23/(1000)  #kJ/mol*K
#kb=1.38e-23
cutoff=L/2 - 1
max_move=0.5
mass = 12.0 #mass of partical
center = [L/2,L/2,L/2] #center of box
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
def lj_potential(r):
    energy = 4*epsilon*(((sigma/r)**12)-((sigma/r)**6))
    return energy
def energy(r,x,y,z,n,compute_all=True):
    r_new = r
    if compute_all:
        r_new[n][0],r_new[n][1],r_new[n][2]=x,y,z
    x_left = list(r_new[:,0]-L)
    x_right = list(r_new[:,0]+L)
    y_left = list(r_new[:,1]-L)
    y_right = list(r_new[:,1]+L)
    z_left = list(r_new[:,2]-L)
    z_right = list(r_new[:,2]+L)
    all_box = np.array([(x_left*3+list(r_new[:,0])*3+x_right*3)*3, \
                        ((y_left+list(r_new[:,1])+y_right)*3)*3, \
                        z_right*9+list(r_new[:,2])*9+z_left*9])
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
        all_box[13*N+n] = [x,y,z]      #update coordinate      
        #print(all_box[13*N**2+n])       
        r0 = np.array([[x for d in range(27*N)],[y for d in range(27*N)],\
                       [z for d in range(27*N)]]) 
        a = r0.transpose((1,0)) - all_box  #compare the updated particle coordinate with others
        distance = [np.sqrt(a[i][0]**2+a[i][1]**2+a[i][2]**2) for i in range(27*N)]  #calculate distance    
        #print(sorted(distance))    
        energy = [lj_potential(d) for d in distance if d < cutoff and d > 0]
        #print(energy)
        energy_sum = sum(energy)
        #print(energy_sum)
    return energy_sum
def move(r,x,y,z,n,e,diff_e,temp):                         
    if diff_e < 0:
        e += diff_e 
        r[n][0],r[n][1],r[n][2] = inbox(x_new,y_new,z_new)
        
    else:
        rand = random.uniform(0,1)
        if math.exp(-diff_e/(kb*temp)) > rand:
            e += diff_e
            r[n][0],r[n][1],r[n][2] = inbox(x_new,y_new,z_new)
    return r,e
def delta(temp1,temp2,energy1,energy2):
    ddelta = (1/(kb*temp1) - 1/(kb*temp2))*(energy2-energy1)
    return ddelta
def exchange(frame,rg,energy,numexchg):
    if numexchg % 2 == 0:
        for i in [[0,1],[2,3],[4,5],[6,7]]:
            dd = delta(temp[i[0]],temp[i[1]],energy[i[0]][(numexchg+1)*(steps+1)-1],energy[i[1]][(numexchg+1)*(steps+1)-1])
            if dd <= 0:
                frame[i[0]][(numexchg+1)*(steps+1)] = frame[i[1]][(numexchg+1)*(steps+1)-1]
                frame[i[1]][(numexchg+1)*(steps+1)] = frame[i[0]][(numexchg+1)*(steps+1)-1]
                energy[i[0]][(numexchg+1)*(steps+1)] = energy[i[1]][(numexchg+1)*(steps+1)-1]
                energy[i[1]][(numexchg+1)*(steps+1)] = energy[i[0]][(numexchg+1)*(steps+1)-1]
                rg[i[0]][(numexchg+1)*(steps+1)] = rg[i[1]][(numexchg+1)*(steps+1)-1]
                rg[i[1]][(numexchg+1)*(steps+1)] = rg[i[0]][(numexchg+1)*(steps+1)-1]
                temp_exchange[numexchg+1][i[0]] = temp_exchange[numexchg][i[1]]
                temp_exchange[numexchg+1][i[1]] = temp_exchange[numexchg][i[0]]
                accepted_exchange[numexchg][i[0]] = 1
                accepted_exchange[numexchg][i[1]] = 1
                accepted_exchange_ratio[numexchg][i[0]] = accepted_exchange[0:numexchg+1,i[0]].sum()/(numexchg+1)
                accepted_exchange_ratio[numexchg][i[1]] = accepted_exchange[0:numexchg+1,i[1]].sum()/(numexchg+1)
            else:
                rand = random.uniform(0,1)
                if math.exp(-dd) > rand:
                    frame[i[0]][(numexchg+1)*(steps+1)] = frame[i[1]][(numexchg+1)*(steps+1)-1]
                    frame[i[1]][(numexchg+1)*(steps+1)] = frame[i[0]][(numexchg+1)*(steps+1)-1]
                    energy[i[0]][(numexchg+1)*(steps+1)] = energy[i[1]][(numexchg+1)*(steps+1)-1]
                    energy[i[1]][(numexchg+1)*(steps+1)] = energy[i[0]][(numexchg+1)*(steps+1)-1]
                    rg[i[0]][(numexchg+1)*(steps+1)] = rg[i[1]][(numexchg+1)*(steps+1)-1]
                    rg[i[1]][(numexchg+1)*(steps+1)] = rg[i[0]][(numexchg+1)*(steps+1)-1]
                    temp_exchange[numexchg+1][i[0]] = temp_exchange[numexchg][i[1]]
                    temp_exchange[numexchg+1][i[1]] = temp_exchange[numexchg][i[0]]
                    accepted_exchange[numexchg][i[0]] = 1
                    accepted_exchange[numexchg][i[1]] = 1
                    accepted_exchange_ratio[numexchg][i[0]] = accepted_exchange[0:numexchg+1,i[0]].sum()/(numexchg+1)
                    accepted_exchange_ratio[numexchg][i[1]] = accepted_exchange[0:numexchg+1,i[1]].sum()/(numexchg+1)
                else:
                    frame[i[0]][(numexchg+1)*(steps+1)] = frame[i[0]][(numexchg+1)*(steps+1)-1]
                    frame[i[1]][(numexchg+1)*(steps+1)] = frame[i[1]][(numexchg+1)*(steps+1)-1]
                    energy[i[0]][(numexchg+1)*(steps+1)] = energy[i[0]][(numexchg+1)*(steps+1)-1]
                    energy[i[1]][(numexchg+1)*(steps+1)] = energy[i[1]][(numexchg+1)*(steps+1)-1]
                    rg[i[0]][(numexchg+1)*(steps+1)] = rg[i[0]][(numexchg+1)*(steps+1)-1]
                    rg[i[1]][(numexchg+1)*(steps+1)] = rg[i[1]][(numexchg+1)*(steps+1)-1]
                    temp_exchange[numexchg+1][i[0]] = temp_exchange[numexchg][i[0]]
                    temp_exchange[numexchg+1][i[1]] = temp_exchange[numexchg][i[1]]
                    accepted_exchange[numexchg][i[0]] = 0
                    accepted_exchange[numexchg][i[1]] = 0
                    accepted_exchange_ratio[numexchg][i[0]] = accepted_exchange[0:numexchg+1,i[0]].sum()/(numexchg+1)
                    accepted_exchange_ratio[numexchg][i[1]] = accepted_exchange[0:numexchg+1,i[1]].sum()/(numexchg+1)
                
    else:
        for i in [[0,7],[1,2],[3,4],[5,6]]:
            dd = delta(temp[i[0]],temp[i[1]],energy[i[0]][(numexchg+1)*(steps+1)-1],energy[i[1]][(numexchg+1)*(steps+1)-1])
            if dd <= 0:
                frame[i[0]][(numexchg+1)*(steps+1)] = frame[i[1]][(numexchg+1)*(steps+1)-1]
                frame[i[1]][(numexchg+1)*(steps+1)] = frame[i[0]][(numexchg+1)*(steps+1)-1]
                energy[i[0]][(numexchg+1)*(steps+1)] = energy[i[1]][(numexchg+1)*(steps+1)-1]
                energy[i[1]][(numexchg+1)*(steps+1)] = energy[i[0]][(numexchg+1)*(steps+1)-1]
                rg[i[0]][(numexchg+1)*(steps+1)] = rg[i[1]][(numexchg+1)*(steps+1)-1]
                rg[i[1]][(numexchg+1)*(steps+1)] = rg[i[0]][(numexchg+1)*(steps+1)-1]
                temp_exchange[numexchg+1][i[0]] = temp_exchange[numexchg][i[1]]
                temp_exchange[numexchg+1][i[1]] = temp_exchange[numexchg][i[0]]
                accepted_exchange[numexchg][i[0]] = 1
                accepted_exchange[numexchg][i[1]] = 1
                
                accepted_exchange_ratio[numexchg][i[0]] = accepted_exchange[0:numexchg+1,i[0]].sum()/(numexchg+1)
                accepted_exchange_ratio[numexchg][i[1]] = accepted_exchange[0:numexchg+1,i[1]].sum()/(numexchg+1)
            else:
                rand = random.uniform(0,1)
                if math.exp(-dd) > rand:
                    frame[i[0]][(numexchg+1)*(steps+1)] = frame[i[1]][(numexchg+1)*(steps+1)-1]
                    frame[i[1]][(numexchg+1)*(steps+1)] = frame[i[0]][(numexchg+1)*(steps+1)-1]
                    energy[i[0]][(numexchg+1)*(steps+1)] = energy[i[1]][(numexchg+1)*(steps+1)-1]
                    energy[i[1]][(numexchg+1)*(steps+1)] = energy[i[0]][(numexchg+1)*(steps+1)-1]
                    rg[i[0]][(numexchg+1)*(steps+1)] = rg[i[1]][(numexchg+1)*(steps+1)-1]
                    rg[i[1]][(numexchg+1)*(steps+1)] = rg[i[0]][(numexchg+1)*(steps+1)-1]
                    temp_exchange[numexchg+1][i[0]] = temp_exchange[numexchg][i[1]]
                    temp_exchange[numexchg+1][i[1]] = temp_exchange[numexchg][i[0]]
                    accepted_exchange[numexchg][i[0]] = 1
                    accepted_exchange[numexchg][i[1]] = 1
                    accepted_exchange_ratio[numexchg][i[0]] = accepted_exchange[0:numexchg+1,i[0]].sum()/(numexchg+1)
                    accepted_exchange_ratio[numexchg][i[1]] = accepted_exchange[0:numexchg+1,i[1]].sum()/(numexchg+1)
                else:
                    frame[i[0]][(numexchg+1)*(steps+1)] = frame[i[0]][(numexchg+1)*(steps+1)-1]
                    frame[i[1]][(numexchg+1)*(steps+1)] = frame[i[1]][(numexchg+1)*(steps+1)-1]
                    energy[i[0]][(numexchg+1)*(steps+1)] = energy[i[0]][(numexchg+1)*(steps+1)-1]
                    energy[i[1]][(numexchg+1)*(steps+1)] = energy[i[1]][(numexchg+1)*(steps+1)-1]
                    rg[i[0]][(numexchg+1)*(steps+1)] = rg[i[0]][(numexchg+1)*(steps+1)-1]
                    rg[i[1]][(numexchg+1)*(steps+1)] = rg[i[1]][(numexchg+1)*(steps+1)-1] 
                    temp_exchange[numexchg+1][i[0]] = temp_exchange[numexchg][i[0]]
                    temp_exchange[numexchg+1][i[1]] = temp_exchange[numexchg][i[1]]
                    accepted_exchange[numexchg][i[0]] = 0
                    accepted_exchange[numexchg][i[1]] = 0
                    accepted_exchange_ratio[numexchg][i[0]] = accepted_exchange[0:numexchg+1,i[0]].sum()/(numexchg+1)
                    accepted_exchange_ratio[numexchg][i[1]] = accepted_exchange[0:numexchg+1,i[1]].sum()/(numexchg+1)

c = Configuration(N) 
#x_init = sorted([i*d for i in range(1,N+1)]*N)
#y_init = [i*d for i in range(1,N+1)]*N
#z_init = [d/2 for i in range(N)]*N
for index in range(N):
    c.atom[index].x = center[0]+ (L-d) * (random.uniform(0,1)-0.5)
    c.atom[index].y = center[1]+ (L-d) * (random.uniform(0,1)-0.5)
    c.atom[index].z = center[2]+ (L-d) * (random.uniform(0,1)-0.5)
    c.atom[index].mass = mass
    c.atom[index].radius = d/2
point = np.array([[c.atom[i].x for i in range(N)],[c.atom[i].y for i in range(N)],
                  [c.atom[i].z for i in range(N)]])
point = point.transpose((1,0))
e = energy(point,0,0,0,0,compute_all=True)
#print(e)
all_energy = np.zeros((len(temp),(numexchg)*(steps+1)+1))
all_energy[:,0] = e
#print(all_energy[0])
all_frame = np.zeros((len(temp),(numexchg)*(steps+1)+1,N,3))
all_frame[:,0] = point
accepted_exchange = np.zeros((numexchg,len(temp)))
accepted_exchange_ratio = np.zeros((numexchg,len(temp)))
temp_exchange = np.zeros((numexchg+1,len(temp)))
temp_exchange[0] = temp
RG = np.zeros((len(temp),(numexchg)*(steps+1)+1))
RG[:,0] = c.RadGyr()
average_energy=np.zeros((len(temp),(numexchg)*(steps+1)+1))
average_energy[:,0] = e
for n in range(numexchg):
    
    for t in range(len(temp)):
        #final_frame=np.zeros((len(temp),N*N,3))
        x_init = all_frame[t][n*(steps+1)][:,0]
        y_init = all_frame[t][n*(steps+1)][:,1]
        z_init = all_frame[t][n*(steps+1)][:,2]
        e = all_energy[t][n*(steps+1)]
        #print(e)
        point = np.array([x_init,y_init,z_init])
        point = point.transpose((1,0))
        for i in range(steps):
            #pre_e = e
            #origi = copy.deepcopy(point)
            for m in range(N): 
                x,y,z = point[m][0],point[m][1],point[m][2]
                #e_old = energy(point,x,y,z,m,compute_all=False)
                #print(e_old)
                x_new = x+max_move*random.uniform(-1,1)
                y_new = y+max_move*random.uniform(-1,1)
                z_new = z+max_move*random.uniform(-1,1)
              # x_new,y_new,z_new = inbox(x_new,y_new,z_new)
                e_new = energy(point,x_new,y_new,z_new,m,compute_all=True)
                diff_e = e_new - e
                #print(diff_e)
                point, e = move(point,x_new,y_new,z_new,m,e,diff_e,temp[t])
                c.atom[m].x,c.atom[m].y,c.atom[m].z = point[m][0],point[m][1],point[m][2]
            
            RG[t][i+1+(steps+1)*n] = c.RadGyr() 
            all_energy[t][i+1+(steps+1)*n] = e
            all_frame[t][i+1+(steps+1)*n] = point      
           # average_energy[t][i+1+(steps+1)*n] = all_energy[t][0:i+2+(steps+1)*n].sum()/(i+2+(steps+1)*n)
      #  final_frame[t] = all_frame[t][-1]
   # print(n)
    exchange(all_frame,RG,all_energy,n)
for t in range(len(temp)): 
    for i in range((numexchg)*(steps+1)+1): 
        average_energy[t][i] = all_energy[t][0:i+1].sum()/(i+1)
tend=datetime.datetime.now()
tdiff=tend-tstart
print(tdiff)
data={}
for i in range(len(temp)):
    data['temp{}'.format(i)] = all_energy[i]
all_energy = pd.DataFrame(data=data,index=range((numexchg)*(steps+1)+1),columns=['temp{}'.format(i) for i in range(len(temp))])
data={}
for i in range(len(temp)):
    data['temp{}'.format(i)] = average_energy[i]
average_energy = pd.DataFrame(data=data,index=range((numexchg)*(steps+1)+1),columns=['temp{}'.format(i) for i in range(len(temp))])
data={}
for i in range(len(temp)):
    data['temp{}'.format(i)] = RG[i]
RG = pd.DataFrame(data=data,index=range((numexchg)*(steps+1)+1))
all_frame_rearrange=np.zeros(((numexchg)*(steps+1)+1,len(temp)*N*3))  
for i in range((numexchg)*(steps+1)+1):
    for t in range(len(temp)):
        for n in range(N):
            for c in range(3):
                
                all_frame_rearrange[i][t*N*3+n*3+c]=all_frame[t][i][n][c]
all_frame = pd.DataFrame(all_frame_rearrange,index=range((numexchg)*(steps+1)+1),columns=['x','y','z']*N*len(temp))
data={}
for i in range(len(temp)): 
    data['walker{}'.format(i)] = temp_exchange[:,i]
all_temp_exchange = pd.DataFrame(data=data,index=range(numexchg+1))
data={}
for i in range(len(temp)):
    data['walker{}'.format(i)] = accepted_exchange[:,i]
accepted_exchange = pd.DataFrame(data=data,index=range(numexchg))
data={}
for i in range(len(temp)):
    data['walker{}'.format(i)] = accepted_exchange_ratio[:,i] 
accepted_exchange_ratio = pd.DataFrame(data=data,index=range(numexchg))
print(average_energy)

#print(temp_exchange)   
#print(accepted_exchange)
#print(accepted_exchange_ratio)
#print(all_energy) 
#print(all_frame)
#print(RG)
#print(all_energy)
RG.to_csv(r'/ufrc/alberto.perezant/liweichang/Tutorial/HW4/remd_4/remd_rg_100.csv',index=False)
all_energy.to_csv(r'/ufrc/alberto.perezant/liweichang/Tutorial/HW4/remd_4/energy_dataframe_100.csv',index=False)
average_energy.to_csv(r'/ufrc/alberto.perezant/liweichang/Tutorial/HW4/remd_4/ave_energy_dataframe_100.csv',index=False)
all_frame.to_csv(r'/ufrc/alberto.perezant/liweichang/Tutorial/HW4/remd_4/frame_dataframe_100.csv',index=False)
accepted_exchange.to_csv(r'/ufrc/alberto.perezant/liweichang/Tutorial/HW4/remd_4/accepted_exchange_100.csv',index=False)
accepted_exchange_ratio.to_csv(r'/ufrc/alberto.perezant/liweichang/Tutorial/HW4/remd_4/accepted_exchange_ratio_100.csv',index=False)
all_temp_exchange.to_csv(r'/ufrc/alberto.perezant/liweichang/Tutorial/HW4/remd_4/temp_ex_dataframe.csv_100',index=False)
