subroutine distribute
use var
use type
use fix
use parallel
implicit none
include 'mpif.h'
integer :: i, j, k, ii, jj, kk, inside, ncpu, zp(0:size), nodesPcpu
real :: RR, Rdroplet


! the box is diveded in the z direction to the number of processors
allocate(Zproc(0:size-1))
allocate(TZproc(0:size-1))
allocate(deltaZ(0:size-1))
allocate(SumZproc(0:size))
allocate(Nproc(0:size-1))
allocate(TNproc(0:size-1))
allocate(deltaN(0:size-1))
allocate(SumNproc(0:size))

If (channel.or.PBC3) then
!check the size and number of processors to optimize distribute
If (int(dime(3)/size)+2>size) then
  Do i=0,size-2
    Zproc(i) = int(dime(3)/size) + 1
  enddo
  Zproc(size-1) = dime(3) - (size-1) * (int(dime(3)/size) + 1) 
else
  Do i=0,size-2
    Zproc(i) = int(dime(3)/(size-1))
  enddo
  Zproc(size-1) = dime(3) - (size-1) * (int(dime(3)/(size-1))) 
endif

Nproc = Zproc * dime(1) * dime(2)
SumNproc = 0
Do i=1,size
   SumNproc(i) = SumNproc(i-1) + Nproc(i-1)
enddo
SumZproc = 0
Do i=1,size
   SumZproc(i) = SumZproc(i-1) + Zproc(i-1)
enddo

MaxZ = maxval(Zproc)
MaxN = maxval(Nproc)
Do i=0,size-1
   TZproc(i) = i*MaxZ
   TNproc(i) = i*MaxN
enddo

Do i=0,size-1
   deltaZ(i) = MaxZ - Zproc(i) 
   deltaN(i) = MaxN - Nproc(i)
enddo

if (ID==master) PRINT *, "Zproc",Zproc
if (ID==master) PRINT *, "TZproc",TZproc
if (ID==master) PRINT *, "sumZproc",sumZproc
if (ID==master) PRINT *, "deltaZ",deltaZ
if (ID==master) PRINT *, "Nproc",Nproc
if (ID==master) PRINT *, "TNproc",TNproc
if (ID==master) PRINT *, "sumNproc",SumNproc

if (ID==master) PRINT *, "Number of nodes per processor =", MaxN  

elseif (droplet) then

if (ID==master) then
center = int(dime/2) + 1
Rdroplet = center(1) - 1
inside = 0
Do k=1,dime(3)
   Do j=1,dime(2)
      Do i=1,dime(1)
         ii = i - center(1)
         jj = j - center(2)
         kk = k - center(3)
         RR = sqrt(real(ii*ii+jj*jj+kk*kk))
         if (RR<=Rdroplet) then
            inside = inside + 1
         endif
      enddo
   enddo
enddo
Print *, inside

nodesPcpu = int(inside/size)
ncpu = 1
inside = 0
zp=0
Do k=1,dime(3)
   Do j=1,dime(2)
      Do i=1,dime(1)
         ii = i - center(1)
         jj = j - center(2)
         kk = k - center(3)
         RR = sqrt(real(ii*ii+jj*jj+kk*kk))
         if (RR<=Rdroplet) then
            inside = inside + 1
         endif
         if (inside>=ncpu*nodesPcpu) then
            Zp(ncpu) = k
            if (j<center(2)) zp(ncpu) = k - 1
            ncpu = ncpu + 1
         endif
      enddo
   enddo
enddo
zp(size) = dime(3)

ncpu = 1
inside = 0
Do k=1,dime(3)
   Do j=1,dime(2)
      Do i=1,dime(1)
         ii = i - center(1)
         jj = j - center(2)
         kk = k - center(3)
         RR = sqrt(real(ii*ii+jj*jj+kk*kk))
         if (RR<=Rdroplet) then
            inside = inside + 1
         endif
         if (k>Zp(ncpu)) then
            inside = 0
            ncpu = ncpu + 1
         endif
      enddo
   enddo
enddo

Do i=0,size-1
   Zproc(i) = Zp(i+1)-Zp(i)
   sumZproc(i) = Zp(i)
enddo

if (ID==master) PRINT *, "Zproc",Zproc
if (ID==master) PRINT *, "TZproc",TZproc
if (ID==master) PRINT *, "sumZproc",sumZproc
if (ID==master) PRINT *, "deltaZ",deltaZ
if (ID==master) PRINT *, "Nproc",Nproc
if (ID==master) PRINT *, "TNproc",TNproc
if (ID==master) PRINT *, "sumNproc",SumNproc

if (ID==master) PRINT *, "Number of nodes per processor =", MaxN  


endif
endif
End subroutine distribute
