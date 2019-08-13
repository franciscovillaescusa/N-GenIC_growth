import numpy as np
import sys

################################## INPUT ######################################
z1=108
z2=90

prefixes = ['CDM', 'B', 'CDMB', 'NU', 'm']

###############################################################################
z = (1+z1 + 1+z2)/2 - 1
print 'Computing growth factors at z =',z

for prefix in prefixes:

    #set name of the files
    Pk_z1 = 'Pk_CAMB/'+prefix+'_Pk_CAMB_z='+str(z1)+'.dat'
    Pk_z2 = 'Pk_CAMB/'+prefix+'_Pk_CAMB_z='+str(z2)+'.dat'
    f_out = 'growth_'+prefix+'_z='+str(z)+'.dat'

    #read power spectra
    k1,Pk1=np.loadtxt(Pk_z1,unpack=True)
    k2,Pk2=np.loadtxt(Pk_z2,unpack=True)

    #check that k-arrays are exactly the same
    if np.any(k1!=k2):
        print 'k-arrays are different!!!!'; sys.exit()

    #compute the growth: f = d log(delta)/d log(a) = - d log(delta)/d log(1+z)
    f = - np.log(np.sqrt(Pk2/Pk1)) / np.log((1.0+z2)/(1.0+z1))

    #save growth to file
    np.savetxt(f_out,np.transpose([k1,f]))

