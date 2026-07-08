!reference:https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.chebwin.html
program main 
integer i,j,k,m,n
integer mltbc_istep
real,parameter   :: pi = 3.14159265
real,parameter   :: dtime = 1800.0 !15.0 
real,allocatable :: mltbc_wgtstep(:)
real sum,sum1,dtwin,dtcut 

dtwin = 3.0 * 3600.0
dtcut = 3.0 * 3600.0 

mltbc_nstep = int(dtwin/dtime) + 1
m = int(mltbc_nstep / 2)

allocate(mltbc_wgtstep(mltbc_nstep))

fcs = 2.0*pi*dtime/dtwin
fcb = pi*dtwin/m/dtcut

do i = 1, mltbc_nstep
  n = i - M - 1
  if ( n == 0  ) then
    mltbc_wgtstep(i) = fcb/pi
  else if ( abs(n) > M ) then
    mltbc_wgtstep(i) = 0.0
  else 
    mltbc_wgtstep(i) = sin(n*pi/(M+1))* (M+1)/(n*pi) * sin(n*fcb)/(n*pi) 
  end if 
  print*,m,mltbc_nstep,n,mltbc_wgtstep(i),sum(mltbc_wgtstep)
end do
mltbc_wgtstep = mltbc_wgtstep / sum(mltbc_wgtstep)
do i = 1,mltbc_nstep
  print*,m,mltbc_nstep,i,mltbc_wgtstep(i),sum(mltbc_wgtstep)
end do 
end 
