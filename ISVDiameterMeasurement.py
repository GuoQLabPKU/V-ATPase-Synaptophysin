#fit a cycle accoring to three points on the edge of circle 
from numpy import *
from scipy import optimize
import functools
from matplotlib import pyplot as p, cm, colors
import pandas as pd
import re
import os
import seaborn as sns
import numpy as np


#measure the vesicle diameter of syp wild type vesicles 
wildTypeDiaList = [ ]

if 1:
    #read the file 
    data_dir = 'C:/Users/DELL/Desktop/Guolab/projects/VATPase/versicleSizeMeasure/versiclesWithV-ATPase/WT/'
    txtDir = os.listdir(data_dir)
    for sgFile in txtDir:
        if '.txt' in sgFile:
            
            sgMicro = pd.read_csv(data_dir + 'ISV_20201210_378_00014_X+1Y+0-2_goldless_addV.txt', names = ['x'])
            #process into readable dataframe 
            xpool = [ ]
            ypool = [ ]
            for sgRow in range(sgMicro.shape[0]):
                splitpool = [float(i) for i in re.findall(r"\d+\.?\d*",sgMicro['x'].values[sgRow])]
                xpool.append(splitpool[0])
                ypool.append(splitpool[1])
                
            sgMicro['x'] = xpool
            sgMicro['y'] = ypool
            #smaple the point each three points 
            assert sgMicro.shape[0]%3 == 0
            

            for i in arange(0, sgMicro.shape[0], 3):
                x = [sgMicro['x'].values[i],  sgMicro['x'].values[i+1], sgMicro['x'].values[i+2]]
                y = [sgMicro['y'].values[i],  sgMicro['y'].values[i+1], sgMicro['y'].values[i+2]]
                
                x_m = mean(x)
                y_m = mean(y)
                def countcalls(fn):
                    "decorator function count function calls "
                
                    @functools.wraps(fn)
                    def wrapped(*args):
                        wrapped.ncalls +=1
                        return fn(*args)
                
                    wrapped.ncalls = 0
                    return wrapped
                
                def calc_R(xc, yc):
                
                    return sqrt((x - xc) ** 2 + (y - yc) ** 2)
                
                @countcalls
                def f_2(c):
                    Ri = calc_R(*c)
                    return Ri - Ri.mean()
                
                def plot_all(xc_2,yc_2, R_2):   
                    theta_fit = linspace(-pi, pi, 180)
                
                    x_fit2 = xc_2 + R_2 * cos(theta_fit)
                    y_fit2 = yc_2 + R_2 * sin(theta_fit)
                    
                    return x_fit2, y_fit2
            
                center_estimate = x_m, y_m
                center_2, _ = optimize.leastsq(f_2, center_estimate)
            
                xc_2, yc_2 = center_2
                Ri_2       = calc_R(xc_2, yc_2)
                R_2        = Ri_2.mean()
                
                
                if R_2 > 1000:
                    print(R_2, x,y )
                    break
                 
                          
                    
                wildTypeDiaList.append(R_2)
            
                

KOTypeDiaList = [ ]
if 0:
    #read the file 
    data_dir = 'C:/Users/DELL/Desktop/Guolab/projects/VATPase/versicleSizeMeasure/versiclesWithV-ATPase/ko/'
    txtDir = os.listdir(data_dir)
    for sgFile in txtDir:
        if '.txt' in sgFile:
            
            sgMicro = pd.read_csv(data_dir + sgFile, names = ['x'])
            #process into readable dataframe 
            xpool = [ ]
            ypool = [ ]
            for sgRow in range(sgMicro.shape[0]):
                splitpool = [float(i) for i in re.findall(r"\d+\.?\d*",sgMicro['x'].values[sgRow])]
                xpool.append(splitpool[0])
                ypool.append(splitpool[1])
                
            sgMicro['x'] = xpool
            sgMicro['y'] = ypool
            #smaple the point each three points 
            assert sgMicro.shape[0]%3 == 0
            
            for i in arange(0, sgMicro.shape[0], 3):
                x = [sgMicro['x'].values[i],  sgMicro['x'].values[i+1], sgMicro['x'].values[i+2]]
                y = [sgMicro['y'].values[i],  sgMicro['y'].values[i+1], sgMicro['y'].values[i+2]]
                
                x_m = mean(x)
                y_m = mean(y)
                def countcalls(fn):
                    "decorator function count function calls "
                
                    @functools.wraps(fn)
                    def wrapped(*args):
                        wrapped.ncalls +=1
                        return fn(*args)
                
                    wrapped.ncalls = 0
                    return wrapped
                
                def calc_R(xc, yc):
                
                    return sqrt((x - xc) ** 2 + (y - yc) ** 2)
                
                @countcalls
                def f_2(c):
                    Ri = calc_R(*c)
                    return Ri - Ri.mean()
                
                def plot_all(xc_2,yc_2, R_2):   
                    theta_fit = linspace(-pi, pi, 180)
                
                    x_fit2 = xc_2 + R_2 * cos(theta_fit)
                    y_fit2 = yc_2 + R_2 * sin(theta_fit)
                    
                    return x_fit2, y_fit2
            
                center_estimate = x_m, y_m
                center_2, _ = optimize.leastsq(f_2, center_estimate)
            
                xc_2, yc_2 = center_2
                Ri_2       = calc_R(xc_2, yc_2)
                R_2        = Ri_2.mean()
                KOTypeDiaList.append(R_2)    

#show and plot the distribution 
wildTypeDiaListnm =  np.array([i*2/10 for i in wildTypeDiaList])
KOTypeDiaListnm = np.array([i*2/10 for i in KOTypeDiaList])
KOTypeDiaListnm = KOTypeDiaListnm[KOTypeDiaListnm>30]


