subroutine input
use var
use type
use fix
implicit none

!--- read input file ---
open(111,file='control',status='old')
read(111,nml = indata)
close(111)

channel=.false.;droplet=.false.;PBC3=.false.
if (system=="channel") then
   channel=.true.
elseif (system=="droplet") then
   droplet=.true.
elseif (system=="bulk") then
   PBC3=.true.
else 
   Print *, "ERORR System type"
   stop
endif

infinite=.false.;nondegenerate=.false.;degenerate=.false.
if (anchor=="inf") then
   infinite=.true.
elseif (anchor=="nond") then
   nondegenerate=.true.
elseif (anchor=="deg") then
   degenerate=.true.
else
   Print *, "ERORR"
   stop
endif

NP_infinite=.false.;NP_nondegenerate=.false.;NP_degenerate=.false.
if (NPanchor=="inf") then
   NP_infinite=.true.
elseif (NPanchor=="nond") then
   NP_nondegenerate=.true.
elseif (NPanchor=="deg") then
   NP_degenerate=.true.
else
   Print *, "ERORR"
   stop
endif


! Initial parameter
Sbulk    = S( ULdG )            ! Scalar Order parameter  
Sinit    = S(2.8d0)
! Allocate variabel
  
tnodes = dime(1)*dime(2)*dime(3)

!if(droplet) qchiral = Pi*chiral/(dble(dime(1)-3))  ! chiral
if(droplet) qchiral = 2.d0*Pi/pitch  ! chiral
!if(channel .or. PBC3) qchiral = Pi*chiral/(dble(box(3)))  ! chiral
if(channel .or. PBC3) qchiral = 2.d0*Pi/pitch  ! chiral
!qchiral = Pi*chiral/(dble(dime(1)-4))  ! chiral Ye

! RANDOM SEED
call RANDOM_SEED
iseed(1:2) = iseed(1:2) 
call RANDOM_SEED(put = iseed(1:2))
  
contains
   Function S(u)
	  use type
	  implicit none
	  real(dp), intent(in):: u
	  real(dp) S
	  S = 0.25d0 + 0.75*sqrt( 1.d0 - 8.d0/(3.d0*U))
   end function S
end subroutine input 
