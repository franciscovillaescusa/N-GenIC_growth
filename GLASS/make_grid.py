# So far Im assuming that N_cdm, N_gas and N_nu will be the same or 0. Need to update this for the
# general case
import numpy as np
import sys,os

################################ INPUT #######################################
fout = 'grid_file_GAS_1_CDM_1.dat'
N_cdm    = 1 #total number of CDM particles in the file will be N_cdm^3
N_gas    = 1 #total number of gas particles in the file will be N_gas^3
N_nu     = 0  #total number of NU  particles in the file will be N_nu^3
##############################################################################

N = N_cdm

# declare different arrays
npart         = np.zeros(6,  dtype=np.int32)
massarr       = np.zeros(6,  dtype=np.float64)
time          = np.zeros(1,  dtype=np.float64)
redshift      = np.zeros(1,  dtype=np.float64)
flag_sfr      = np.zeros(1,  dtype=np.int32)
flag_feedback = np.zeros(1,  dtype=np.int32)
npartall      = np.zeros(6,  dtype=np.int32)
flag_cooling  = np.zeros(1,  dtype=np.int32)
num_files     = np.ones(1,   dtype=np.int32)
BoxSize       = np.zeros(1,  dtype=np.float64);  BoxSize[0] = 100.0 
extra_array   = np.zeros(30, dtype=np.float32) #120 bytes

# fill the number of particles arrays
npart[1], npartall[1] = N_gas**3, N_gas**3 #GAS
npart[2], npartall[2] = N_cdm**3, N_cdm**3 #CDM
npart[3], npartall[3] = N_nu**3,  N_nu**3  #NU

# declare the arrays hosting the particle positions
pos_gas = np.zeros((N_gas**3, 3), dtype=np.float32)
pos_cdm = np.zeros((N_cdm**3, 3), dtype=np.float32)
pos_nu  = np.zeros((N_nu**3,  3), dtype=np.float32)

# N-body
if N_gas==0:
    
    # find the particle positions
    for i in xrange(0,N):
        for j in xrange(0,N):
            for k in xrange(0,N_cdm):
                pos_cdm[(i*N+j)*N+k, 0] = (i+0.0)*BoxSize/N
                pos_cdm[(i*N+j)*N+k, 1] = (j+0.0)*BoxSize/N
                pos_cdm[(i*N+j)*N+k, 2] = (k+0.0)*BoxSize/N

    if N_nu>0:
        # find the NU particle positions
        pos_nu[:] = pos_cdm[:] + 0.5*BoxSize/N
                
                
# hydro without neutrinos
elif N_nu==0 and N_gas>0:

    # find the gas particle positions
    for i in xrange(0,N):
        for j in xrange(0,N):
            for k in xrange(0,N):
                pos_gas[(i*N+j)*N+k, 0] = (i+0.0)*BoxSize/N
                pos_gas[(i*N+j)*N+k, 1] = (j+0.0)*BoxSize/N
                pos_gas[(i*N+j)*N+k, 2] = (k+0.0)*BoxSize/N

    # find the CDM particle positions
    pos_cdm[:] = pos_gas[:] + 0.5*BoxSize/N


# hydro with neutrinos
elif N_gas>0 and N_nu>0:
    
    # find the gas,cdm and nu particle positions
    for i in xrange(0,N):
        for j in xrange(0,N):
            for k in xrange(0,N):
                pos_gas[(i*N+j)*N+k, 0] = (i+0.5)*BoxSize/N
                pos_gas[(i*N+j)*N+k, 1] = (j+0.5)*BoxSize/N
                pos_gas[(i*N+j)*N+k, 2] = (k+0.0)*BoxSize/N

                pos_cdm[(i*N+j)*N+k, 0] = (i+0.5)*BoxSize/N
                pos_cdm[(i*N+j)*N+k, 1] = (j+0.0)*BoxSize/N
                pos_cdm[(i*N+j)*N+k, 2] = (k+0.5)*BoxSize/N

                pos_nu[(i*N+j)*N+k, 0]  = (i+0.0)*BoxSize/N
                pos_nu[(i*N+j)*N+k, 1]  = (j+0.5)*BoxSize/N
                pos_nu[(i*N+j)*N+k, 2]  = (k+0.5)*BoxSize/N

                
# arrays with the number of bytes of each block
blocksize1 = np.array([256],                                       dtype=np.int32)
blocksize2 = np.array([4*3*N_gas**3 + 4*3*N_cdm**3 + 4*3*N_nu**3], dtype=np.int32)

# write the binary file
f = open(fout, 'wb')
for array in [blocksize1, npart, massarr, time, redshift, flag_sfr,
              flag_feedback, npartall, flag_cooling, num_files,
              BoxSize, extra_array, blocksize1,
              blocksize2, pos_gas, pos_cdm, pos_nu, blocksize2]:
    array.tofile(f)
f.close()
