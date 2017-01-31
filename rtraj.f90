subroutine rtraj(step)
use var
use type
use parallel
implicit none
integer,intent(out)  :: step
integer :: i, k, j, N, NN, label, ioerr, outQ
character*20 :: string
real, allocatable :: qsingle(:,:)

allocate (qsingle(6,Mnodes))
write(*,*)"READ THE TRAJ FILE"
write(777,*)"READ THE TRAJ FILE"
outQ = 200 + ID

do
   read(outQ,IOSTAT=ioerr)string
   if (ioerr == 0) then
      read(outQ)frame,step
      read(outQ)qsingle
   else
      exit
   endif
enddo

do i=1,Mnodes
   q_reduced(:,i) = qsingle(:,i) 
enddo

deallocate(qsingle)
write(*,*)"FRAME",frame,"STEP",step
write(777,*)"FRAME",frame,"STEP",step
end subroutine

