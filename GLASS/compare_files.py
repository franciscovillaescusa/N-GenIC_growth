# This script can be used to compare two grid files. I used it to validate that
# my python script produces exactly the same result as the IDL script of Matteo
import numpy as np
import sys,os

f1 = 'grid_file_2comp_64.dat'
f2 = 'dummy_glass_CDM_B_64_64.dat'

f = open(f1, 'rb')
g = open(f2, 'rb')

#blocksize
d1 = np.fromfile(f, dtype=np.int32, count=1)
d2 = np.fromfile(g, dtype=np.int32, count=1)
print d1;  print d2; print  ' '

#npart
d1 = np.fromfile(f, dtype=np.int32, count=6)
d2 = np.fromfile(g, dtype=np.int32, count=6)
print d1;  print d2; print  ' '

#massarr
d1 = np.fromfile(f, dtype=np.float64, count=6)
d2 = np.fromfile(g, dtype=np.float64, count=6)
print d1;  print d2; print  ' '

#time
d1 = np.fromfile(f, dtype=np.float64, count=1)
d2 = np.fromfile(g, dtype=np.float64, count=1)
print d1;  print d2; print  ' '

#redshift
d1 = np.fromfile(f, dtype=np.float64, count=1)
d2 = np.fromfile(g, dtype=np.float64, count=1)
print d1;  print d2; print  ' '

#flag_sfr
d1 = np.fromfile(f, dtype=np.int32, count=1)
d2 = np.fromfile(g, dtype=np.int32, count=1)
print d1;  print d2; print  ' '

#flag_feedback
d1 = np.fromfile(f, dtype=np.int32, count=1)
d2 = np.fromfile(g, dtype=np.int32, count=1)
print d1;  print d2; print  ' '

#npartall
d1 = np.fromfile(f, dtype=np.int32, count=6)
d2 = np.fromfile(g, dtype=np.int32, count=6)
print d1;  print d2; print  ' '

# flag cooling
d1 = np.fromfile(f, dtype=np.int32, count=1)
d2 = np.fromfile(g, dtype=np.int32, count=1)
print d1;  print d2; print  ' '

#num_files
d1 = np.fromfile(f, dtype=np.int32, count=1)
d2 = np.fromfile(g, dtype=np.int32, count=1)
print d1;  print d2; print  ' '

#BoxSize
d1 = np.fromfile(f, dtype=np.float64, count=1)
d2 = np.fromfile(g, dtype=np.float64, count=1)
print d1;  print d2; print  ' '

#extra array
d1 = np.fromfile(f, dtype=np.int32, count=30)
d2 = np.fromfile(g, dtype=np.int32, count=30)
print d1;  print d2; print  ' '

#blocksize
d1 = np.fromfile(f, dtype=np.int32, count=1)
d2 = np.fromfile(g, dtype=np.int32, count=1)
print d1;  print d2; print  ' '

#blocksize
d1 = np.fromfile(f, dtype=np.int32, count=1)
d2 = np.fromfile(g, dtype=np.int32, count=1)
print d1;  print d2; print  ' '

#pos1
d1 = np.fromfile(f, dtype=np.float32, count=64**3*3)
d2 = np.fromfile(g, dtype=np.float32, count=64**3*3)
print d1;  print d2; print  ' '

#pos2
d1 = np.fromfile(f, dtype=np.float32, count=64**3*3)
d2 = np.fromfile(g, dtype=np.float32, count=64**3*3)
print d1;  print d2; print  ' '

#blocksize
d1 = np.fromfile(f, dtype=np.int32, count=1)
d2 = np.fromfile(g, dtype=np.int32, count=1)
print d1;  print d2; print  ' '

f.close();  g.close()
