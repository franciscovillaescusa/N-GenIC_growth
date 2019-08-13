import numpy as np
import sys


########################## create Hz file ###################################  
f_in  = '0.0eV_hubble.txt'  #file provided by Matteo  
f_out = 'Hz.txt'            #file to be used by Gadget  
bins  = 1000                #number of bins in the new file  
z_max = 99.0                #maximum redshift of the new file  
z_min = 0.0                 #minimum redshift of the new file  

#read original file; H(z) in km/s/Mpc   
z,Hz = np.loadtxt(f_in,unpack=True)

#create a new (1+z) array    
z_new = np.logspace(np.log10(1.0+z_min),np.log10(1.0+z_max),bins)

#interpolate from the original table. The H(z) has to be in Gadget units:  
#km/s/(kpc/h), thus we need to multiply the Hz of Matteo by 1.0/(1000.0*h)    
Hz_new = np.interp(z_new-1,z,Hz)/(Hz[0]/100.0)/1000.0

#create the w=-1 array     
w = np.ones(bins)*(-1.0)

#the format of the H(z) file is 1+z,w,H(z), where 1+z should be decreasing   
np.savetxt(f_out,np.transpose([z_new[::-1],w,Hz_new[::-1]]))
############################################################################# 
