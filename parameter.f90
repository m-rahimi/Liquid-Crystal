subroutine define_parameter
use var
use type
implicit none
real(dp) :: dr, NP_volume
integer :: i,j,k,label,NPin

df = dble(box)/dble(dime-1) ! distance between two nodes
dr = df(1) ! Correct if NP density equal in all direction
df = 1/df ! Invers of df
df2 = df*0.5 ! Invers of 2*df
ddf = df*df

nodes = 0
bulk_master = .false.
do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
         if (channel.or.PBC3) then
            if (boundary_master(i,j,k)==.false.) then
               bulk_master(i,j,k) = .true.
               nodes = nodes + 1
            endif
         elseif (drop(i,j,k)) then
            bulk_master(i,j,k) = .true.
            nodes = nodes + 1
         endif
      enddo
   enddo
enddo

! differential 
if (channel.or.PBC3) then
!if(channel)  Volume = (box(1)*box(2)*box(3)) - number_NP*((4.0/3.0)*Pi*((radius)*dr)**3) 
if(channel)  Volume = (dime(1)*dime(2)*(dime(3)-2.0)) - number_NP*((4.0/3.0)*Pi*((radius)*dr)**3) 
if(PBC3)  Volume = (dime(1)*dime(2)*dime(3)) - number_NP*((4.0/3.0)*Pi*((radius)*dr)**3) 
   ad = 2.d0*(dime(1)*dime(2)) / bnodes 
   if (number_NP>0) then
      ad_NP = number_NP*(4.d0*Pi*((radius)*dr)**2) / tNPbnodes
   else
      ad_NP = 1 !to avoid division by zero 
   endif 
elseif (droplet) then
   Volume = (4.d0/3.d0)*Pi*((center(1)-2)*dr)**3
   ad = 4.d0*Pi*((center(1)-2)*dr)**2 / bnodes
!   Volume = (4.d0/3.d0)*Pi*(( dime(1)/2.d0-2.d0)*dr)**3
!   ad = 4.d0*Pi*(( dime(1)/2.d0-2.d0)*dr)**2 / bnodes
endif
vd = Volume / nodes
print *, "vd is =", vd
print *, "q =", qchiral
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Eneargy Scale

Escale = Elastic_L1 * nematic_l * 1.0d-9 ! (J/m)*(nm)*(10^-9)

Print *, "Total nodes             = ", tnodes 
Print *, "Total bulk nodes        = ", nodes 
Print *, 'Total boundary nodes    = ', bnodes
Print *, 'Total NP boundary nodes = ', tNPbnodes
Print *, 'Total nodes inside NP   = ', NPin
end subroutine  define_parameter


