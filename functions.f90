subroutine position_NP()
use type
use var
use parallel 
implicit none
integer :: coun = 0, r, i, nn
real :: temp_pos(3)
real(dp) :: dr(3), rand(3)
logical :: flag, file_exists
Character(len=10) :: tt

R = radius + 1
!Check size of the NP
if ((2*R>=dime(1).or.2*R>=dime(2).or.2*R>=dime(3)).and.number_NP>0) then
   if (ID==master) Print *, "NP is bigger than box"
   stop
endif

! Check for the NP_position file
INQUIRE(FILE="NP_position", EXIST=file_exists)
if (file_exists) then
   if (ID==master) Print *,"read the position of NPs from NP_position file"
   open(123,file='NP_position', status = 'old' )
   read(123,*) nn
   if (nn/=number_NP) then
      if (ID==master) Print *, "Number of NPs in the file is not equal to the number NPs in the control"
      stop
   endif
   read(123,*)
   Do i=1,number_NP
      read(123,*) tt,NP_position(:,i)
   enddo
else

! Locate 1 NP in the center of the box
if (number_NP==1) then
   NP_position(:,1) = int(dime/2) + 1
elseif (number_NP > 1) then !for more than one NP
Do while (coun <= number_NP)
   call random_number(rand)
   temp_pos = int(rand*dime) + 1
   temp_pos(3) = int(dime(3)/2) + 1 ! Fix NP in Z direction
   if (coun == 0) then
      coun = coun + 1
      NP_position(:,coun) = temp_pos
   else
      flag = .true.
      do i=1 , coun
          dr = (temp_pos - NP_position(:,i))
          dr = dr - dime * anint( dr / dime )
          dr = dr * dr
          if (sqrt(sum(dr))<2*R) flag = .false.
      enddo
      if (flag) then
         coun = coun + 1
         NP_position(:,coun) = temp_pos
      endif
   endif
enddo
endif
endif

end subroutine position_NP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine director_tensor( s, n, qq )
  use type
  implicit none
  real(dp), intent(in) :: s
  real(dp), intent(in),  dimension(3) :: n
  real(dp), intent(out), dimension(5) :: qq
  ! this definition uses the SCALAR ORDER PARAMETER
  !   not the maximum eigenvalue 
  qq(1) = s*( n(1)*n(1) - 1.d0/3.d0 )
  qq(2) = s*( n(1)*n(2) )
  qq(3) = s*( n(1)*n(3) )
  qq(4) = s*( n(2)*n(2) - 1.d0/3.d0 )
  qq(5) = s*( n(2)*n(3) )
end subroutine director_tensor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transformation( q, a )
  use type
  implicit none
  real(dp), dimension(5), intent(in) :: q
  real(dp), dimension(5), intent(out) :: a
  
  a = 0.d0
  a(1) = - sqrt(3.d0/2.d0)*( q(1) + q(4) )
  a(2) = sqrt(2.d0)*q(2)
  a(3) = sqrt(2.d0)*q(3)
  a(4) = ( q(1) - q(4) )/sqrt(2.d0)
  a(5) = sqrt(2.d0)*q(5) 

end subroutine transformation
!----------------------------------------------------------
subroutine eigen( qq, eigenvalues, vector )
  use type
  implicit none
  real(dp), intent(in),  dimension(5) :: qq
  real(dp) :: dumpy
  real(dp), intent(out), dimension(3) :: vector
  integer(i4b) lwork, info,i,j
  real(dp), dimension(3)   :: eigenvalues
  real(dp), dimension(10)  :: work
  real(dp), dimension(3,3) :: a

  lwork  = 10
  ! Matrix form
  a(1,1) = qq(1);   a(1,2) = qq(2);   a(1,3) = qq(3)
  a(2,1) = a(1,2); a(2,2) = qq(4);   a(2,3) = qq(5)
  a(3,1) = a(1,3); a(3,2) = a(2,3); a(3,3) = - a(1,1) - a(2,2)
  call dsyev( 'V', 'U', 3, a, 3, eigenvalues, work, lwork, info )
!  value  = maxval( eigenvalues )
  vector = a(:,3)
  
  dumpy =  eigenvalues(1)
  eigenvalues(1) =  eigenvalues(3)
  eigenvalues(3) = dumpy
 
end subroutine eigen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inv_transformation( a, q )
  use type
  implicit none
  real(dp), dimension(5), intent(in) :: a
  real(dp), dimension(5), intent(out) :: q
  real(dp) q33, aux1, aux2

  aux1 = 1.d0/sqrt(6.d0) ; aux2 = 1.d0/sqrt(2.d0)
  q = 0.d0
  q(1) = -a(1)*aux1 + a(4)*aux2
  q(2) =  a(2)*aux2
  q(3) =  a(3)*aux2
  q(4) = -a(1)*aux1 - a(4)*aux2
  q(5) =  a(5)*aux2
  q33 = sqrt(2.d0/3.d0)*a(1)
  q33 = q33 + q(1) + q(4) 
!  if ( abs(q33) > 1.d-10 )stop('Trace is non-zero')

end subroutine inv_transformation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_bound( q_tem, flag )
  use type
  use var

  implicit none
  logical, intent(out) :: flag
  real(dp) valor, tr
  real(dp), dimension(3) :: vector, dump
  real(dp), dimension(3,3) :: dumi
  real(dp), dimension(6) :: q_tem

  flag = .true.

  ! EigenValues
  call eigen( q_tem(1:5), vector, dump )
  vector(1) =     vector(1)*1.5d0
  vector(2) =     vector(2) + vector(1)/3.d0
  vector(3) = - ( vector(3) + vector(1)/3.d0 )

  valor = ( 1.d0 - vector(1) )/3.d0
  
  ! Boundness for the tensor and its eigenvalues
  if ( vector(1) < -0.5d0 .or. vector(1) > 1.d0 )   flag = .false.
  if ( vector(2) < - valor .or. vector(2) > valor ) flag = .false.

  ! Boundness for the trace --> Included Correct limits from Orlando's arguments
  tr =  2*q_tem(1)**2 + 2*q_tem(2)**2 + 2*q_tem(3)**2 + 2*q_tem(4)**2 + 2*q_tem(5)**2 + 2*q_tem(1)*q_tem(4)

  valor = 2.d0*( 3.d0*vector(2)*vector(2)+ vector(1)*vector(1))/3.d0
 
  if ( abs( tr - valor ) > 1.d-5 ) then
     print *,' Trace Q2 is not the correct value '
     print *, tr, valor, abs( tr - valor )
     stop
  end if
  if ( tr > 2.d0/3.d0 ) flag = .false. 

end subroutine check_bound
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine timer(clock_finish,clock_start)
  double precision :: clock_finish,clock_start
  if((clock_finish-clock_start)/60 < 60.d0)then
     write (*,999)(clock_finish-clock_start)/60.d0
     write (777,999)(clock_finish-clock_start)/60.d0
  else
     write (*,998)(clock_finish-clock_start)/3600.d0
     write (777,998)(clock_finish-clock_start)/3600.d0
  endif
999 format (' Total time of CPU  :',f26.10,'  minutes')
998 format (' Total time of CPU  :',f26.10,'  hours')
end subroutine timer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_bound2( a_tem, flag )
  use type
  use fix
  implicit none
  logical, intent(out) :: flag
  real(dp),dimension(5), intent(in) :: a_tem
  real(dp) valor, tr
  real(dp), dimension(3) :: vector
  real(dp), dimension(3,3) :: dumi
!

real(dp) , dimension(3) :: eigenvalue
real(dp)  ::  tr2

real(dp)  tr3, aux1, aux2, u, v, theta,sqrtu

  flag = .true.

  ! EigenValues

  tr = trq(a_tem) ;

  u =   tr  /6.d0  ;
  v = - trq3(a_tem) /6.d0 ;

  theta = v/sqrt( u**3 ) ;
  sqrtu = -2.d0*sqrt(u)

  if ( theta > 1.d0 ) theta = 1.d0
  if ( theta < -1.d0 ) theta = -1.d0
  theta = acos( theta )

  vector(1) = sqrtu*cos( ( theta + 2*pi )/3.d0 )
  vector(2) = sqrtu*cos( ( theta - 2*pi )/3.d0 )
  vector(3) = sqrtu*cos( theta/3.d0)
  vector(1) =     vector(1)*1.5d0
  vector(2) =     vector(2) + vector(1)/3.d0
  vector(3) = - ( vector(3) + vector(1)/3.d0 )

  valor = ( 1.d0 - vector(1) )/3.d0

  ! Boundness for the tensor and its eigenvalues
  if ( vector(1) < -0.5d0 .or. vector(1) > 1.d0 )   flag = .false.
  if ( vector(2) < - valor .or. vector(2) > valor ) flag = .false.

  ! Boundness for the trace --> Included Correct limits from Orlando's arguments

  valor = 2.d0*( 3.d0*vector(2)*vector(2)+ vector(1)*vector(1))/3.d0

  if ( abs( tr - valor ) > 1.d-5 ) then
     print *,' Trace Q2 is not the correct value '
     print *, tr, valor, abs( tr - valor )
     stop
  end if
  if ( tr > 2.d0/3.d0 ) flag = .false.
contains
function trq (a)
  use type
  implicit none
  real(dp), dimension(5) :: a
  real(dp)  :: trq
  trq = a(1)**2 + a(2)**2 + a(3)**2 + a(4)**2 + a(5)**2
  return
end function trq
!----------------------------------------------------------
function trqt (a,typ)
  use type
  implicit none
  real(dp), dimension(5) :: a
  real(dp)  :: trqt
  integer :: typ
  trqt = a(typ)**2
  return
end function trqt
!----------------------------------------------------------
!----------------------------------------------------------
function trq3 (a)
  use type
  implicit none
  real(dp), dimension(5) :: a
  real(dp)  :: trq3,aux1,aux2,aux3

  aux1 = 1.d0/sqrt(6.d0)
  aux2 = sqrt(1.5d0)
  aux3 = sqrt(3.d0)

  trq3 = a(1)**3*aux1 - aux2*a(1)*a(2)**2 + 0.5d0*aux2*a(1)*a(3)**2 + 3.d0*a(3)**2*a(4)/2.d0/sqrt(2.d0) &
       - aux2*a(1)*a(4)**2 + 3.d0*a(2)*a(3)*a(5)/sqrt(2.d0) + 0.5d0*aux2*a(1)*a(5)**2 - 3.d0*a(4)*a(5)**2/2.d0/sqrt(2.d0)

  return
end function trq3

end subroutine check_bound2

