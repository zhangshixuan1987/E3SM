program main 

integer i,j,k,N
real lambda(6)
real pi,sum
integer mltbc_istep
real mltbc_wgtstep
pi = 3.14159265
dtime = 1800.0
dtwin = 3.0 * 3600.0

mltbc_nstep = int(dtwin/dtime)
M = int(mltbc_nstep / 2)
sum = 0.0 
do mltbc_istep  = 1, mltbc_nstep 
  k = 0
  do i = 1,M
     k = k + 2*i
  end do
  if (mltbc_istep <= M ) then
     mltbc_wgtstep = mltbc_istep * 1.0 / k
  else
     mltbc_wgtstep = (mltbc_nstep - mltbc_istep + 1.0) / k
  end if
  print*,'mltbc_wgtstep=',mltbc_wgtstep
  sum = sum + mltbc_wgtstep
end do 
print*,sum
end 
