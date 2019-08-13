#This code reads the matterpower spectra and the transfer functions from CAMB
#and computes the different power spectra

import numpy as np
import sys,os


################################## INPUT ######################################
root_CAMB     = '../CAMB_TABLES/ics_matterpow_'
root_transfer = '../CAMB_TABLES/ics_transfer_'

Omega_CDM = 0.2685
Omega_B   = 0.049

redshifts = [108, 99, 90, 80, 70, 60, 50, 40, 30, 20, 10,
             9, 8, 7, 6, 5, 4, 3, 2, 1, 0.5, 0]
###############################################################################

print '\n\n####################### WARNING #########################'
print 'Values of Omega_CDM and Omega_B used to compute Pk_CDMB:'
print 'Omega_CDM =',Omega_CDM
print 'Omega_B   =',Omega_B
print '#########################################################\n\n'


#do a loop over all redshifts and compute the different power spectra
for z in redshifts:

    #obtain the name of the matter power spectrum and transfer function files
    f_Pk       = root_CAMB     + str(z) + '.dat' 
    f_transfer = root_transfer + str(z) + '.dat' 

    #copy the matter power spectrum files to this folder
    os.system('cp '+f_Pk+' m_Pk_CAMB_z='+str(z)+'.dat')

    # read CAMB matter power spectrum file
    k_m,Pk_m = np.loadtxt(f_Pk,unpack=True)

    # read CAMB transfer function file
    k,Tcdm,Tb,dumb,dumb,Tnu,Tm = np.loadtxt(f_transfer,unpack=True)

    #compute the CDM+B transfer function
    Tcdmb = (Omega_CDM*Tcdm+Omega_B*Tb)/(Omega_CDM+Omega_B)

    #Interpolate to find P(k)_matter in the same ks as the transfer functions
    Pk_m = 10**(np.interp(np.log10(k),np.log10(k_m),np.log10(Pk_m)))

    #set the name of the output files
    f_c          = 'CDM_Pk_CAMB_z='           + str(z) + '.dat'
    f_b          = 'B_Pk_CAMB_z='             + str(z) + '.dat'
    f_n          = 'NU_Pk_CAMB_z='            + str(z) + '.dat'
    f_cb         = 'CDMB_Pk_CAMB_z='          + str(z) + '.dat'
    f_cross_c_b  = 'CDM_B_cross_Pk_CAMB_z='   + str(z) + '.dat'
    f_cross_c_n  = 'CDM_NU_cross_Pk_CAMB_z='  + str(z) + '.dat'
    f_cross_b_n  = 'B_NU_cross_Pk_CAMB_z='    + str(z) + '.dat'
    f_cross_cb_n = 'CDMB_NU_cross_Pk_CAMB_z=' + str(z) + '.dat'

    #compute the different power spectra and save them
    Pk_c  = Pk_m*(Tcdm/Tm)**2;   np.savetxt(f_c ,np.transpose([k,Pk_c]))
    Pk_b  = Pk_m*(Tb/Tm)**2;     np.savetxt(f_b ,np.transpose([k,Pk_b]))
    Pk_n  = Pk_m*(Tnu/Tm)**2;    np.savetxt(f_n ,np.transpose([k,Pk_n]))
    Pk_cb = Pk_m*(Tcdmb/Tm)**2;  np.savetxt(f_cb,np.transpose([k,Pk_cb]))

    Pk_x_c_b  = Pk_m*Tcdm*Tb/Tm**2;   np.savetxt(f_cross_c_b,
                                                 np.transpose([k,Pk_x_c_b]))
    Pk_x_c_n  = Pk_m*Tcdm*Tnu/Tm**2;  np.savetxt(f_cross_c_n,
                                                 np.transpose([k,Pk_x_c_n]))
    Pk_x_b_n  = Pk_m*Tb*Tnu/Tm**2;    np.savetxt(f_cross_b_n,
                                                 np.transpose([k,Pk_x_b_n]))
    Pk_x_cb_n = Pk_m*Tcdmb*Tnu/Tm**2; np.savetxt(f_cross_cb_n,
                                                 np.transpose([k,Pk_x_cb_n]))
    




