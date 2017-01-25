Program main
use type
use var
use parallel
implicit none
include 'mpif.h'
integer :: i ,step, outQ, outE, sum
real (dp) :: temp_E(7)
real (dp) :: clock_start,clock_finish,c_ini,c_fin 
character (len=8) :: fname='Energy', fname2='Q', fname3='lQ'

call MPI_INIT(ierr)
CALL MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, MPI_COMM_NEW, ierr)
call MPI_COMM_SIZE(MPI_COMM_NEW, size, ierr)
call MPI_COMM_RANK(MPI_COMM_NEW, ID, ierr)

call input()

! master generate initial configuration
if (ID==master) then
  call cpu_time(clock_start)
  Print *, "Number of processors = ", size

  open(777,file='data.out')
  outQ = 200 + ID
  outE = 300 + ID
  if (restart == 0) then
     frame = 0
     step = 0
     open ( outQ , FILE = trim(adjustr(fname2)),form='UNFORMATTED', access='SEQUENTIAL')
     open ( outE , FILE = trim(adjustr(fname)),STATUS='replace')
     write( outE,*) "@ step LdG Elastic Energy"
  else
     open ( outQ , FILE = trim(adjustr(fname2)),form='UNFORMATTED', access='SEQUENTIAL',status='old')
     open ( outE , FILE = trim(adjustr(fname)),STATUS='old')
     read ( outE,*)
  endif

  call allocate_master()
  call initial()
  call define_parameter()
  call convert()
 !if (number_NP > 0) then
!   allocate (NP_position(3,number_NP))
!   call position_NP()
!!   if (droplet) call make_np_droplet()
!   if (channel) call make_np_channel()
!   call MPI_BARRIER(MPI_COMM_NEW,ierr)
!endif
!
!
if (restart == 1) then
   call rtraj(step) ! keep the old style
   Print *, "Simulation continue from step ", step
   Do while (i/=step)
      read(outE,'(I8,2X,5(f15.6),2(f20.6))') i
   enddo
   frame = frame + 1
endif
endif

! Distribute nodes between all processors
call MPI_BARRIER(MPI_COMM_NEW,ierr)
call scater()

!!!Gibzburg - Landau
CALL GL(step)

call MPI_BARRIER(MPI_COMM_NEW,ierr)
if (ID==master) then
  call cpu_time (clock_finish)
  call timer(clock_finish,clock_start)
endif
call MPI_BARRIER(MPI_COMM_NEW,ierr)
call MPI_FINALIZE(ierr)
end program main
