subroutine wtraj(step)
use var
use type
use parallel
implicit none
integer, intent(in) :: step
integer :: i, j, k, N, NN, label , outQ
character*20 :: string = "traj frame"
real, allocatable :: qsingle(:,:)

allocate (qsingle(6,Mnodes))
outQ = 200 + ID

do i=1,Mnodes
   qsingle(:,i) = real(q(:,i))
enddo

write(outQ)string
write(outQ)frame,step
write(outQ)qsingle
frame = frame + 1
deallocate (qsingle)
end subroutine

