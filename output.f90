Subroutine output_files(step)
use type
use var
use parallel
implicit none
include 'mpif.h'
integer :: i, j, k, indx, outlQ, Nsend, Nrec(size), displas(size+1)
integer :: N, NN
integer, intent(in) :: step
real(dp), dimension(5) ::   q_tmp
real(dp),dimension(3) :: temp_scalar,temp_director
real(dp),dimension(:,:), allocatable :: Tscalar, Tdirector
real(dp),dimension(:,:), allocatable :: scalar, director
real(sp),dimension(:,:), allocatable :: qsingle

If (output>=1.and.ID==master.and.step/=0) then
   outlQ = 400 + ID
   open ( outlQ , FILE = "lQ" ,form='UNFORMATTED', access='SEQUENTIAL')
   allocate (qsingle(6,Mnodes))

   do i=1,Mnodes
      qsingle(:,i) = real(q(:,i))
   enddo

   write(outlQ)"traj frame"
   write(outlQ)frame,step
   write(outlQ)qsingle

   deallocate(qsingle)
endif

if (output>=2) then

if (ID==master) allocate(Tscalar(3,length*size), Tdirector(3,length*size))
allocate(scalar(3,length), director(3,length))
scalar = 0.d0
director = 0.d0

do indx=1,length
   q_tmp = q(1:5,indx)
   call eigen ( q_tmp(:), temp_scalar, temp_director(:))
   scalar(:,indx) = temp_scalar
   if (temp_director(3)<0.0) temp_director = -1.d0 * temp_director
   director(:,indx) = temp_director
enddo

CALL MPI_BARRIER(MPI_COMM_NEW, ierr)
CALL MPI_GATHER (scalar, 3*length, MPI_DOUBLE_PRECISION, Tscalar, 3*length, MPI_DOUBLE_PRECISION, master, MPI_COMM_NEW, IERR)
CALL MPI_GATHER (director, 3*length, MPI_DOUBLE_PRECISION, Tdirector, 3*length, MPI_DOUBLE_PRECISION, master, MPI_COMM_NEW, IERR)

if (ID==master) call paraview(step)
if (ID==master .and. step /= 0) call CrossPol
CALL MPI_BARRIER(MPI_COMM_NEW, ierr)
if (ID==master) deallocate(Tscalar, Tdirector)
deallocate(scalar, director)
endif

contains
subroutine paraview(step)
use type
use var
implicit none
integer, intent(in) :: step
integer (i4b) :: i,iout=50,j,k,label
real(dp) :: r2,lbox,cx,x,y,z,eta,splaybend_tmp, temp_scalar(3), temp_director(3)
!********************************************************************
Character(len=10) :: namefile1, namefile2
Character(len=20), Parameter :: sfile_3= 'surface',end='vtk'  

if (step==0) then
   open(iout,FILE=trim(adjustr(sfile_3))//"0."//trim(adjustl(end)),STATUS='replace')
else
   open(iout,FILE=trim(adjustr(sfile_3))//"1."//trim(adjustl(end)),STATUS='replace')
endif

!*************************************************************************************

write(iout,100)'# vtk DataFile Version 3.0'
write(iout,100)'vtk output'
write(iout,100)'ASCII'
write(iout,100)'DATASET STRUCTURED_GRID'
write(iout,110)'DIMENSIONS',dime(1),dime(2),dime(3)
write(iout,120)'POINTS',dime(1)*dime(2)*dime(3) ,'  float'

do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
          write(iout,151)  real(i),real(j),real(k),'  '
     enddo
   enddo
enddo


write(iout,*)
write(iout,120)'POINT_DATA', dime(1)*dime(2)*dime(3)
write(iout,100)'SCALARS scalars float'
write(iout,100)'LOOKUP_TABLE default'


do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
            label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
            if (lb2indx(label)==0) then
               temp_scalar = 0.0
            else
               temp_scalar = Tscalar(:,lb2indx(label))
            endif
            write(iout,152 )  maxval(temp_scalar)*1.5d0,'  '
      enddo
   enddo
enddo

write(iout,*)
write(iout,100)'SCALARS biaxial float'
write(iout,100)'LOOKUP_TABLE default'


do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
            label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
            if (lb2indx(label)==0) then
               temp_scalar = 0.0
               temp_director = 0.0
            else
               temp_scalar = Tscalar(:,lb2indx(label))
               temp_director = Tdirector(:,lb2indx(label))
            endif
            write(iout,152 ) temp_scalar(2)*3.d0 + maxval(temp_scalar)*1.5d0,'  ' 
      enddo
   enddo
enddo


write(iout,*)
write(iout,100)'VECTORS directorsz float'

do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
            label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
            if (lb2indx(label)==0) then
               temp_director = 0.0
            else
               temp_director = Tdirector(:,lb2indx(label))
            endif
!            write(iout,153 )  abs(temp_director(1)),abs(temp_director(2)),abs(temp_director(3)),'  '
            write(iout,153 )  (temp_director(1)),(temp_director(2)),(temp_director(3)),'  '
      enddo
   enddo
enddo

  write(iout,*)
  write(iout,100)'SCALARS SplayBend float'
  write(iout,100)'LOOKUP_TABLE default'

  do k=1,dime(3)
     do j=1,dime(2)
        do i=1,dime(1)
              label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
               if (lb2indx(label)==0) then
               splaybend_tmp = -1.0
               else
              !call Sbend (lb2indx(label),splaybend_tmp)
              call Sbend (label,splaybend_tmp)
              ! splaybend_tmp=splaybend(lb2indx(label))
               endif
              write(iout,152 ) splaybend_tmp,'  '
        enddo
     enddo
  enddo




close (iout)
!*****************************************************************************************
!*****************************************************************************************
100 format('',A)
120 format('',A,I12,A)
110 format('',A, I6, 1X, I6, 1X, I6)
151 format (3f12.6 ,A)
152 format (f12.6 ,A)
153 format (3f12.6 ,A)
end subroutine paraview

subroutine CrossPol
use type
use var
implicit none
integer :: i,j,k,l,label,step,iout11=50,iout21=50,iout31=50,iout41=50
integer :: iout12=50,iout22=50,iout32=50,iout42=50
integer :: iout13=50,iout23=50,iout33=50,iout43=50
integer :: iout14=50,iout24=50,iout34=50,iout44=50
integer :: iout15=50,iout25=50,iout35=50,iout45=50
integer :: iout16=50,iout26=50,iout36=50,iout46=50
real(dp) :: LdG = 0.d0, Elastic = 0.d0, surface0 = 0.d0, surface1 = 0.d0, Chiral=0.d0
!real(dp), intent(out) :: energy(6)
real(dp), dimension(5) :: q_tmp
real(dp), dimension(2) :: S11,S12,S21,S22,S11new,S12new,S21new,S22new,S11old,S12old,S21old,S22old
real(dp) :: lambda(6), no, ne, dmesh,ggama,phio,phie,bbeta,nebeta,scalar(3),Pvalue,gama1,gama0
  real(dp),dimension(:,:), allocatable :: Tscalar, Tdirector,Tspbend
  real(dp),dimension(3) :: dx,director,xs
real(dp) :: thetax, x,y,z
real(dp)  xxr(3), rot(3,3),xx(3)

!********************************************************************
Character(len=10) :: namefile1, namefile2
Character(len=20), Parameter :: sfile_3= 'view',end='dat'

!*************************************************************************************

no=1.5d0; ne=1.7d0; dmesh=1.d0

do l=1,6

lambda(1)=10.d0 
lambda(2)=20.d0 
lambda(3)=30.d0 
lambda(4)=40.d0 
lambda(5)=50.d0 
lambda(6)=60.d0 

if(l.eq.1) open(unit=iout11,FILE=trim(adjustr(sfile_3))//"zBlue."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.2) open(unit=iout12,FILE=trim(adjustr(sfile_3))//"zGreen."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.3) open(unit=iout13,FILE=trim(adjustr(sfile_3))//"zRed."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.4) open(unit=iout14,FILE=trim(adjustr(sfile_3))//"zBlue2."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.5) open(unit=iout15,FILE=trim(adjustr(sfile_3))//"zGreen2."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.6) open(unit=iout16,FILE=trim(adjustr(sfile_3))//"zRed2."//trim(adjustl(end)),action='write',STATUS='replace')

S11old=1.d0; S12old=0.d0; S21old=0.d0; S22old=1.d0

Pvalue=0.d0
gama0=0.d0
  do i=1,dime(1)
           gama0=0.d0
     do j=1,dime(2)
           gama0=0.d0

        do k=1,dime(3)
!            write(iout,151)  real(i),real(j),real(k),'  '

            if (boundary_master(i,j,k).or.bulk_master(i,j,k)) then
!            print *, i,j,k, boundary_master(i,j,k), bulk_master(i,j,k)
              label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
              q_tmp = q(1:5,lb2indx(label))
              call eigen ( q_tmp(:), scalar, director(:))
              Tscalar(:,lb2indx(label)) = scalar
              Tdirector(:,lb2indx(label)) = director
!  bbeta=acos(abs(director(3))); ggama=atan2(director(2),director(1));
  bbeta=acos(abs(director(3))); gama1=atan2(director(2),director(1)); ggama=gama1-gama0; gama0=gama1; 
phio=2.d0*pi*no*dmesh/lambda(l)
  nebeta=(no*ne/sqrt(no**2*sin(bbeta)**2+ne**2*cos(bbeta)**2));
phie=2.d0*pi*nebeta*dmesh/lambda(l)

!  S11 = cos(amma)**2* exp(ii* phie) + sin(gama)**2*exp(ii* phio)
!  S12 = cos(gama)* sin(gama) *(exp(ii* phie) - exp(ii* phio));
!  S21 = cos(gama)* sin(gama) *(exp(ii* phie) - exp(ii* phio));
!  S22 = sin(gmma)**2 exp(ii* phie) + cos(gama)**2* exp(ii* phio);

S11(1)= cos(ggama)**2 * cos(phie) + sin(ggama)**2* cos(phio)
S11(2)= cos(ggama)**2 * sin(phie) + sin(ggama)**2* sin(phio)
S12(1)= cos(ggama)* sin(ggama)*(cos(phie)-cos(phio) )
S12(2)= cos(ggama)* sin(ggama)*(sin(phie)-sin(phie) )
S21(1)= cos(ggama)* sin(ggama)*(cos(phie)-cos(phio) )
S21(2)= cos(ggama)* sin(ggama)*(sin(phie)-sin(phie) )
S22(1)= sin(ggama)**2 * cos(phie) + cos(ggama)**2* cos(phio)
S22(2)= sin(ggama)**2 * sin(phie) + cos(ggama)**2* sin(phio)

!!!!!Product of matrix S^2 !!!!!!!!!!

S11new(1)=(S11(1)*S11old(1)-S11(2)*S11old(2))+(S12(1)*S12old(1)-S12(2)*S12old(2))
S11new(2)=(S11(1)*S11old(2)+S11old(1)*S11(2))+(S12(1)*S12old(2)+S12old(1)*S12(2))

S12new(1)=(S11old(1)*S21(1)-S11old(2)*S21(2))+(S21old(1)*S22(1)-S21old(2)*S22(2))
S12new(2)=(S11old(1)*S21(2)+S11old(2)*S21(1))+(S21old(1)*S22(2)+S21old(2)*S22(1))
S21new(1)=(S11(1)*S21old(1)-S11(2)*S21old(2))+(S21(1)*S22old(1)-S21(2)*S22old(2))
S21new(2)=(S11(1)*S21old(2)+S11(2)*S21old(1))+(S21(1)*S22old(2)+S21(2)*S22old(1))

S22new(1)=(S22(1)*S22old(1)-S22(2)*S22old(2))+(S12(1)*S12old(1)-S12(2)*S12old(2))
S22new(2)=(S22(1)*S22old(2)+S22old(1)*S22(2))+(S12(1)*S12old(2)+S12old(1)*S12(2))


S11old=S11new
S12old=S12new
S21old=S21new
S22old=S22new

!            print *, i,j,k, boundary_master(i,j,k), bulk_master(i,j,k)
           else 
           gama0=0.d0
           endif
                  if(k .eq. dime(3)) then
                  Pvalue=S12new(1)**2+S12new(2)**2
                 ! Pvalue=(S12new(1)+S22new(1))**2+(S12new(2)+S22new(1))**2
                 ! write(24,*) i,j,Pvalue
if(l.eq.1) write(iout11,*) i,j,Pvalue
if(l.eq.2) write(iout12,*) i,j,Pvalue
if(l.eq.3) write(iout13,*) i,j,Pvalue
if(l.eq.4) write(iout14,*) i,j,Pvalue
if(l.eq.5) write(iout15,*) i,j,Pvalue
if(l.eq.6) write(iout16,*) i,j,Pvalue
                 ! write(iout1,*) i,j,Pvalue
                  S11old=1.d0; S12old=0.d0; S21old=0.d0; S22old=1.d0
Pvalue=0.d0
                  gama0=0.d0
                  endif
       enddo
     enddo
  enddo

close (iout11)
close (iout12)
close (iout13)
close (iout14)
close (iout15)
close (iout16)

!  open(unit=iout2,FILE=trim(adjustr(sfile_3))//"x."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.1) open(unit=iout21,FILE=trim(adjustr(sfile_3))//"xBlue."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.2) open(unit=iout22,FILE=trim(adjustr(sfile_3))//"xGreen."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.3) open(unit=iout23,FILE=trim(adjustr(sfile_3))//"xRed."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.4) open(unit=iout24,FILE=trim(adjustr(sfile_3))//"xBlue2."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.5) open(unit=iout25,FILE=trim(adjustr(sfile_3))//"xGreen2."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.6) open(unit=iout26,FILE=trim(adjustr(sfile_3))//"xRed2."//trim(adjustl(end)),action='write',STATUS='replace')
S11old=1.d0; S12old=0.d0; S21old=0.d0; S22old=1.d0
Pvalue=0.d0
  do k=1,dime(3)
gama0=0.d0
     do j=1,dime(2)
gama0=0.d0
        do i=1,dime(1)
!            write(iout,151)  real(i),real(j),real(k),'  '

Pvalue=0.d0
!            if (boundary_master(i,j,k).or.drop(i,j,k)) then
            if (boundary_master(i,j,k).or.bulk_master(i,j,k)) then
              label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
              q_tmp = q(1:5,lb2indx(label))
              call eigen ( q_tmp(:), scalar, director(:))
              Tscalar(:,lb2indx(label)) = scalar
              Tdirector(:,lb2indx(label)) = director
!  bbeta=acos(abs(director(1))); ggama=atan2(director(3),director(2));
  bbeta=acos(abs(director(1))); gama1=atan2(director(3),director(2)); ggama=gama1-gama0; gama0=gama1; 
phio=2.d0*pi*no*dmesh/lambda(l)
  nebeta=(no*ne/sqrt(no**2*sin(bbeta)**2+ne**2*cos(bbeta)**2));
phie=2.d0*pi*nebeta*dmesh/lambda(l)

!  S11 = cos(amma)**2* exp(ii* phie) + sin(gama)**2*exp(ii* phio)
!  S12 = cos(gama)* sin(gama) *(exp(ii* phie) - exp(ii* phio));
!  S21 = cos(gama)* sin(gama) *(exp(ii* phie) - exp(ii* phio));
!  S22 = sin(gmma)**2 exp(ii* phie) + cos(gama)**2* exp(ii* phio);

S11(1)= cos(ggama)**2 * cos(phie) + sin(ggama)**2* cos(phio)
S11(2)= cos(ggama)**2 * sin(phie) + sin(ggama)**2* sin(phio)
S12(1)= cos(ggama)* sin(ggama)*(cos(phie)-cos(phio) )
S12(2)= cos(ggama)* sin(ggama)*(sin(phie)-sin(phie) )
S21(1)= cos(ggama)* sin(ggama)*(cos(phie)-cos(phio) )
S21(2)= cos(ggama)* sin(ggama)*(sin(phie)-sin(phie) )
S22(1)= sin(ggama)**2 * cos(phie) + cos(ggama)**2* cos(phio)
S22(2)= sin(ggama)**2 * sin(phie) + cos(ggama)**2* sin(phio)

!!!!!Product of matrix S^2 !!!!!!!!!!

S11new(1)=(S11(1)*S11old(1)-S11(2)*S11old(2))+(S12(1)*S12old(1)-S12(2)*S12old(2))
S11new(2)=(S11(1)*S11old(2)+S11old(1)*S11(2))+(S12(1)*S12old(2)+S12old(1)*S12(2))
S12new(1)=(S11old(1)*S21(1)-S11old(2)*S21(2))+(S21old(1)*S22(1)-S21old(2)*S22(2))
!S12new(2)=(S11old(1)*S21(2)+S11old(2)*S21(1))+(S21old(1)*S22(2)-S21old(2)*S22(1))
S12new(2)=(S11old(1)*S21(2)+S11old(2)*S21(1))+(S21old(1)*S22(2)+S21old(2)*S22(1))

S21new(1)=(S11(1)*S21old(1)-S11(2)*S21old(2))+(S21(1)*S22old(1)-S21(2)*S22old(2))
!S21new(2)=(S11(1)*S21old(2)+S11(2)*S21old(1))+(S21(1)*S22old(2)-S21(2)*S22old(1))
S21new(2)=(S11(1)*S21old(2)+S11(2)*S21old(1))+(S21(1)*S22old(2)+S21(2)*S22old(1))

S22new(1)=(S22(1)*S22old(1)-S22(2)*S22old(2))+(S12(1)*S12old(1)-S12(2)*S12old(2))
S22new(2)=(S22(1)*S22old(2)+S22old(1)*S22(2))+(S12(1)*S12old(2)+S12old(1)*S12(2))


S11old=S11new
S12old=S12new
S21old=S21new
S22old=S22new
           else
gama0=0.d0
           endif
                  if(i .eq. dime(1)) then
                  Pvalue=S12new(1)**2+S12new(2)**2
                  !Pvalue=(S12new(1)+S22new(1))**2+(S12new(2)+S22new(1))**2
if(l.eq.1) write(iout21,*) j,k,Pvalue
if(l.eq.2) write(iout22,*) j,k,Pvalue
if(l.eq.3) write(iout23,*) j,k,Pvalue
if(l.eq.4) write(iout24,*) j,k,Pvalue
if(l.eq.5) write(iout25,*) j,k,Pvalue
if(l.eq.6) write(iout26,*) j,k,Pvalue
                  !write(25,*) j,k,Pvalue
                  !write(iout2,*) j,k,Pvalue
                  S11old=1.d0; S12old=0.d0; S21old=0.d0; S22old=1.d0
Pvalue=0.d0
gama0=0.d0
                  endif
       enddo
     enddo
  enddo
close (iout21)
close (iout22)
close (iout23)
close (iout24)
close (iout25)
close (iout26)
!  open(unit=iout3,FILE=trim(adjustr(sfile_3))//"lambda"lambda(l)"y."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.1) open(unit=iout31,FILE=trim(adjustr(sfile_3))//"yBlue."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.2) open(unit=iout32,FILE=trim(adjustr(sfile_3))//"yGreen."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.3) open(unit=iout33,FILE=trim(adjustr(sfile_3))//"yRed."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.4) open(unit=iout34,FILE=trim(adjustr(sfile_3))//"yBlue2."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.5) open(unit=iout35,FILE=trim(adjustr(sfile_3))//"yGreen2."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.6) open(unit=iout36,FILE=trim(adjustr(sfile_3))//"yRed2."//trim(adjustl(end)),action='write',STATUS='replace')
S11old=1.d0; S12old=0.d0; S21old=0.d0; S22old=1.d0
Pvalue=0.d0
  do k=1,dime(3)
gama0=0.d0
     do i=1,dime(1)
gama0=0.d0
        do j=1,dime(2)
!            write(iout,151)  real(i),real(j),real(k),'  '

     !       if(boundary_master(i,j,k).or.drop(i,j,k)) then
            if (boundary_master(i,j,k).or.bulk_master(i,j,k)) then
              label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
              q_tmp = q(1:5,lb2indx(label))
              call eigen ( q_tmp(:), scalar, director(:))
              Tscalar(:,lb2indx(label)) = scalar
              Tdirector(:,lb2indx(label)) = director
 ! bbeta=acos(abs(director(2))); ggama=atan2(director(1),director(3));
  bbeta=acos(abs(director(2))); gama1=atan2(director(1),director(3)); ggama=gama1-gama0; gama0=gama1; 
phio=2.d0*pi*no*dmesh/lambda(l)
  nebeta=(no*ne/sqrt(no**2*sin(bbeta)**2+ne**2*cos(bbeta)**2));
phie=2.d0*pi*nebeta*dmesh/lambda(l)

!  S11 = cos(amma)**2* exp(ii* phie) + sin(gama)**2*exp(ii* phio)
!  S12 = cos(gama)* sin(gama) *(exp(ii* phie) - exp(ii* phio));
!  S21 = cos(gama)* sin(gama) *(exp(ii* phie) - exp(ii* phio));
!  S22 = sin(gmma)**2 exp(ii* phie) + cos(gama)**2* exp(ii* phio);

S11(1)= cos(ggama)**2 * cos(phie) + sin(ggama)**2* cos(phio)
S11(2)= cos(ggama)**2 * sin(phie) + sin(ggama)**2* sin(phio)
S12(1)= cos(ggama)* sin(ggama)*(cos(phie)-cos(phio) )
S12(2)= cos(ggama)* sin(ggama)*(sin(phie)-sin(phie) )
S21(1)= cos(ggama)* sin(ggama)*(cos(phie)-cos(phio) )
S21(2)= cos(ggama)* sin(ggama)*(sin(phie)-sin(phie) )
S22(1)= sin(ggama)**2 * cos(phie) + cos(ggama)**2* cos(phio)
S22(2)= sin(ggama)**2 * sin(phie) + cos(ggama)**2* sin(phio)

!!!!!Product of matrix S^2 !!!!!!!!!!

S11new(1)=(S11(1)*S11old(1)-S11(2)*S11old(2))+(S12(1)*S12old(1)-S12(2)*S12old(2))
S11new(2)=(S11(1)*S11old(2)+S11old(1)*S11(2))+(S12(1)*S12old(2)+S12old(1)*S12(2))
S12new(1)=(S11old(1)*S21(1)-S11old(2)*S21(2))+(S21old(1)*S22(1)-S21old(2)*S22(2))
S12new(2)=(S11old(1)*S21(2)+S11old(2)*S21(1))+(S21old(1)*S22(2)+S21old(2)*S22(1))

S21new(1)=(S11(1)*S21old(1)-S11(2)*S21old(2))+(S21(1)*S22old(1)-S21(2)*S22old(2))
S21new(2)=(S11(1)*S21old(2)+S11(2)*S21old(1))+(S21(1)*S22old(2)+S21(2)*S22old(1))

S22new(1)=(S22(1)*S22old(1)-S22(2)*S22old(2))+(S12(1)*S12old(1)-S12(2)*S12old(2))
S22new(2)=(S22(1)*S22old(2)+S22old(1)*S22(2))+(S12(1)*S12old(2)+S12old(1)*S12(2))


S11old=S11new
S12old=S12new
S21old=S21new
S22old=S22new
           else
          
gama0=0.d0
           endif
                  if(j .eq. dime(2)) then
                  Pvalue=S12new(1)**2+S12new(2)**2
                  !Pvalue=(S12new(1)+S22new(1))**2+(S12new(2)+S22new(1))**2
                  !write(26,*) k,i,Pvalue
if(l.eq.1) write(iout31,*) k,i,Pvalue
if(l.eq.2) write(iout32,*) k,i,Pvalue
if(l.eq.3) write(iout33,*) k,i,Pvalue
if(l.eq.4) write(iout34,*) k,i,Pvalue
if(l.eq.5) write(iout35,*) k,i,Pvalue
if(l.eq.6) write(iout36,*) k,i,Pvalue
                  !write(iout3,*) k,i,Pvalue
                  S11old=1.d0; S12old=0.d0; S21old=0.d0; S22old=1.d0
Pvalue=0.d0
                  gama0=0.d0
                  endif
       enddo
     enddo
  enddo

!close (iout1)
!close (iout2)
close (iout31)
close (iout32)
close (iout33)
close (iout34)
close (iout35)
close (iout36)
!!!Rotated system!!!!!

!  open(unit=iout4,FILE=trim(adjustr(sfile_3))//"lambda"lambda(l)"zR."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.1) open(unit=iout41,FILE=trim(adjustr(sfile_3))//"zRBlue."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.2) open(unit=iout42,FILE=trim(adjustr(sfile_3))//"zRGreen."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.3) open(unit=iout43,FILE=trim(adjustr(sfile_3))//"zRRed."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.4) open(unit=iout44,FILE=trim(adjustr(sfile_3))//"zRBlue2."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.5) open(unit=iout45,FILE=trim(adjustr(sfile_3))//"zRGreen2."//trim(adjustl(end)),action='write',STATUS='replace')
if(l.eq.6) open(unit=iout46,FILE=trim(adjustr(sfile_3))//"zRRed2."//trim(adjustl(end)),action='write',STATUS='replace')
S11old=1.d0; S12old=0.d0; S21old=0.d0; S22old=1.d0
Pvalue=0.d0
  do i=1,dime(1)
     do j=1,dime(2)
        do k=1,dime(3)

       !     if(boundary_master(i,j,k).or.drop(i,j,k)) then
            if (boundary_master(i,j,k).or.bulk_master(i,j,k)) then
              label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)

                  x = real(i); y = real(j); z = real(k)

                    thetax=pi/4.d0

                    rot(1,1) =  cos(thetax)
                    rot(1,2) = -sin(thetax)
                    rot(1,3) = 0.0
                    rot(2,1) = sin(thetax)
                    rot(2,2) = cos(thetax)
                    rot(2,3) = 0.0
                    rot(3,1) = 0.0
                    rot(3,2) = 0.0
                    rot(3,3) = 1.0

                    xx(1) = x
                    xx(2) = y
                    xx(3) = z

                    xxr=matmul(rot,xx)

                    x = xxr(1)
                    y = xxr(2)
                    z = xxr(3)


              q_tmp = q(1:5,lb2indx(label))
              call eigen ( q_tmp(:), scalar, director(:))
              Tscalar(:,lb2indx(label)) = scalar
              Tdirector(:,lb2indx(label)) = director


                   director=matmul(rot,director) 
                      


!bbeta=acos(abs(director(3))); ggama=atan2(director(2),director(1));
  bbeta=acos(abs(director(3))); gama1=atan2(director(2),director(1)); ggama=gama1-gama0; gama0=gama1; 
phio=2.d0*pi*no*dmesh/lambda(l)
nebeta=(no*ne/sqrt(no**2*sin(bbeta)**2+ne**2*cos(bbeta)**2));
phie=2.d0*pi*nebeta*dmesh/lambda(l)


S11(1)= cos(ggama)**2 * cos(phie) + sin(ggama)**2* cos(phio)
S11(2)= cos(ggama)**2 * sin(phie) + sin(ggama)**2* sin(phio)
S12(1)= cos(ggama)* sin(ggama)*(cos(phie)-cos(phio) )
S12(2)= cos(ggama)* sin(ggama)*(sin(phie)-sin(phie) )
S21(1)= cos(ggama)* sin(ggama)*(cos(phie)-cos(phio) )
S21(2)= cos(ggama)* sin(ggama)*(sin(phie)-sin(phie) )
S22(1)= sin(ggama)**2 * cos(phie) + cos(ggama)**2* cos(phio)
S22(2)= sin(ggama)**2 * sin(phie) + cos(ggama)**2* sin(phio)

!!!!!Product of matrix S^2 !!!!!!!!!!

S11new(1)=(S11(1)*S11old(1)-S11(2)*S11old(2))+(S12(1)*S12old(1)-S12(2)*S12old(2))
S11new(2)=(S11(1)*S11old(2)+S11old(1)*S11(2))+(S12(1)*S12old(2)+S12old(1)*S12(2))

S12new(1)=(S11old(1)*S21(1)-S11old(2)*S21(2))+(S21old(1)*S22(1)-S21old(2)*S22(2))
S12new(2)=(S11old(1)*S21(2)+S11old(2)*S21(1))+(S21old(1)*S22(2)+S21old(2)*S22(1))
S21new(1)=(S11(1)*S21old(1)-S11(2)*S21old(2))+(S21(1)*S22old(1)-S21(2)*S22old(2))
S21new(2)=(S11(1)*S21old(2)+S11(2)*S21old(1))+(S21(1)*S22old(2)+S21(2)*S22old(1))

S22new(1)=(S22(1)*S22old(1)-S22(2)*S22old(2))+(S12(1)*S12old(1)-S12(2)*S12old(2))
S22new(2)=(S22(1)*S22old(2)+S22old(1)*S22(2))+(S12(1)*S12old(2)+S12old(1)*S12(2))


S11old=S11new
S12old=S12new
S21old=S21new
S22old=S22new
           endif
                  if(k .eq. dime(3)) then
                 Pvalue=S12new(1)**2+S12new(2)**2
                  !Pvalue=(S12new(1)+S22new(1))**2+(S12new(2)+S22new(1))**2
                  !write(28,*) x,y,Pvalue
if(l.eq.1) write(iout41,*) x,y,Pvalue
if(l.eq.2) write(iout42,*) x,y,Pvalue
if(l.eq.3) write(iout43,*) x,y,Pvalue
if(l.eq.4) write(iout44,*) x,y,Pvalue
if(l.eq.5) write(iout45,*) x,y,Pvalue
if(l.eq.6) write(iout46,*) x,y,Pvalue
                  !write(iout4,*) x,y,Pvalue
                  S11old=1.d0; S12old=0.d0; S21old=0.d0; S22old=1.d0
Pvalue=0.d0
                  endif
       enddo
     enddo
  enddo

!close (iout1)
!close (iout2)
!close (iout3)
close (iout41)
close (iout42)
close (iout43)
close (iout44)
close (iout45)
close (iout46)

enddo !!!do for lambda

151 format (3f12.6 ,A)

End subroutine  CrossPol








end subroutine output_files
