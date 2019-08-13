# We use this script to generate N-GenIC parameter files
# for standard and paired simulations. It only needs a N-GenIC
# fiducial parameter file where Seed is set to 1340, Phase_flip
# to 9998 and RayleighSampling to 4447
import numpy as np
import sys,os

############################### INPUT ###################################
fiducial_file = 'N-GenIC.param'

realizations = 100

fiducial_seed = 1340
fiducial_flip = 9998
fiducial_RS   = 4447
#########################################################################


for i in xrange(realizations):

    seed = fiducial_seed + 5*i

    # create the folder if it does not exists
    folder1 = str(i)
    folder2 = 'NCV_0_'+str(i)
    folder3 = 'NCV_1_'+str(i)
    for folder in [folder1, folder2, folder3]:
        if not(os.path.exists(folder)):  os.system('mkdir %s'%folder)

        # open input and output files
        fin  = open(fiducial_file, 'r')
        fout = open(folder+'/%s'%fiducial_file, 'w')

        if folder==folder1:
            for line in fin:
                if str(fiducial_seed) in line.split():
                    fout.write(line.replace(str(fiducial_seed), 
                                            str(seed)))
                elif str(fiducial_RS) in line.split():
                    fout.write(line.replace(str(fiducial_RS),
                                            str(1)))
                elif str(fiducial_flip) in line.split():
                    fout.write(line.replace(str(fiducial_flip), str(0)))
                else:
                    fout.write(line)

        if folder==folder2:
            for line in fin:
                if str(fiducial_seed) in line.split():
                    fout.write(line.replace(str(fiducial_seed), 
                                            str(seed)))
                elif str(fiducial_RS) in line.split():
                    fout.write(line.replace(str(fiducial_RS),
                                            str(0)))
                elif str(fiducial_flip) in line.split():
                    fout.write(line.replace(str(fiducial_flip), str(0)))
                else:
                    fout.write(line)

        if folder==folder3:
            for line in fin:
                if str(fiducial_seed) in line.split():
                    fout.write(line.replace(str(fiducial_seed), 
                                            str(seed)))
                elif str(fiducial_RS) in line.split():
                    fout.write(line.replace(str(fiducial_RS),
                                            str(0)))
                elif str(fiducial_flip) in line.split():
                    fout.write(line.replace(str(fiducial_flip), str(1)))
                else:
                    fout.write(line)
        
        fin.close(); fout.close()


