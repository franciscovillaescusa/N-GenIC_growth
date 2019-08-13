

N  = 64L

fout = "dummy_glass_2comp_64.dat"


Ntot = N * N * N


npart=lonarr(6)	
massarr=dblarr(6)
time=0.0D
redshift=0.0D
flag_sfr=0L
flag_feedback=0L
npartall=lonarr(6)	
flag_cooling= 0L
num_files= 1L  
BoxSize = 0.0D

bytesleft=120
la=intarr(bytesleft/2)



BoxSize = 100.0D


npart(1) = Ntot
npartall(1) = Ntot

npart(2) = Ntot
npartall(2) = Ntot


pos= fltarr(3, Ntot)

for i=0L, N-1 do begin
  for j=0L, N-1 do begin
    for k=0L, N-1 do begin
      pos(0, (i*N+j)*N+k) = (i+0.0)/N * BoxSize
      pos(1, (i*N+j)*N+k) = (j+0.0)/N * BoxSize
      pos(2, (i*N+j)*N+k) = (k+0.0)/N * BoxSize
    endfor
  endfor
endfor


pos2 = pos + float((0.5)/N * BoxSize)

openw,1,fout,/f77_unformatted
writeu,1, npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall,flag_cooling,num_files,BoxSize,la
writeu,1, pos, pos2
close,1

end


