subroutine GL(step)
use var
use type
use parallel
implicit none
include 'mpif.h'
integer, intent(inout) :: step
integer :: i, outE, outlQ
real(dp) :: delta, dt_cur
real(dp) :: energy(6), E1, E2, DE, ac
real, allocatable :: qsingle(:,:)
logical :: flag = .true. , restart_flag=.true.

allocate(q_new(6,length))
outE = 300 + ID

Call MPI_BCAST  (step,1,MPI_INT,master,MPI_COMM_NEW,ierr)
delta = (dt(2)-dt(1))/dt(3)
dt_cur = dt(1) 

Call output_files(step)

if (restart==1) restart_flag=.false.

flag = .true.
Do while (flag)
   if (mod(step,traj)==0.and.ID==master.and.restart_flag) then
      write(777,'(I8,2X,5(f15.6),3(f20.6))')step,energy(1:6),Escale*energy(6)/kbT,ac
      call wtraj(step)
   endif
   if (mod(step,sample)==0) then
   call free_energy(energy)
      ! calculate accuracy; energy between two step
      E2 = energy(6)
      DE = (E1-E2)
      ac = abs(DE/E1)*100
      E1 = E2

      if (ID==master) then
        open (22, file = 'timestep')
        write(22,'(I8,2X,5(f15.6),2(f20.6),(f15.10))')step,energy(1:6),Escale*energy(6)/kbT, ac
        write(outE,'(I8,2X,5(f15.6),2(f20.6))')step,energy(1:6),Escale*energy(6)/kbT
        write(*,'(I8,2X,5(f15.6),2(f20.6),(f15.10))')step,energy(1:6),Escale*energy(6)/kbT, ac
        close(22)
      endif

      if (ac < accuracy) flag = .false.

      dt_cur = dt(1) + delta*step
      if (step>dt(3)) dt_cur = dt(2)

      restart_flag=.true.
   endif
   
   call MPI_BARRIER(MPI_COMM_NEW,ierr)
   call  ginzburg_landau(dt_cur)
   step = step + 1
enddo

deallocate(q_new) !to save the memory
if (ID==master) write(777,'(I8,2X,5(f15.6),2(f20.6),(f15.10))')step,energy(1:6),Escale*energy(6)/kbT,ac

call MPI_BARRIER(MPI_COMM_NEW,ierr)
Call output_files(step)

endsubroutine GL
