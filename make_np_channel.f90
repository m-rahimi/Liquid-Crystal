Subroutine make_np_channel
use type
use var
use parallel
implicit none
include 'mpif.h'
integer :: i,label, indx!, deltax, deltay, deltaz
integer, dimension(:,:), allocatable :: temp_po 
integer, dimension(:), allocatable :: temp_NP

allocate (temp_po(4,int(4*radius**3)))

NPnodes = 0
NPbnodes = 0
NP = .false.
boundary_NP = .false.

call grid_NP()
call find_normal()

call MPI_Win_fence(0, win, ierr)
call anchoring_q()
call MPI_Win_fence(0, win, ierr)
CALL MPI_BARRIER(MPI_COMM_NEW, ierr)

call find_boundary_neighbors()

!if (ID==0) call qmga(0,0)
deallocate(temp_po)
!deallocate(share,bshare)
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine grid_NP()
use var
use parallel
implicit none
integer(i4b) :: i,j,k, ii, kk 
integer(i4b) :: cont_b = 0 , cont_np = 0, cont_sh = 0 ,R
real(dp) :: checkx, checky, checkz, check1, check2, check3
real(dp) :: Lx, Ly, Lz
real, dimension(3) :: xs

allocate(global_NP(dime(1),dime(2),dime(3)),global_boundary_NP(dime(1),dime(2),dime(3)))
R = radius

Do ii=1,number_NP !loop for number of NPs
xs = NP_position(:,ii) 

! select all nodes located in the NPs
do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)

           Lx = i - xs(1)
           Ly = j - xs(2)
           Lz = k - xs(3)

           checkx = (Lx)/real((R-0.5));
           checkx = checkx*checkx;
           checky = (Ly)/real((R-0.5));
           checky = checky*checky;
           checkz = (Lz)/real((R-0.5));
           checkz = checkz*checkz;
           check2 = checkx+checky+checkz;
           
           if (check2<= 1.d0) then
              global_NP(i,j,k) = .true.
              cont_np = cont_np + 1
           endif
      enddo
   enddo
enddo
enddo ! number of NP

! select all nodes which are neighbor of NP and choose them as boundary nodes
do k=2,dime(3)-1
   do j=2,dime(2)-1
      do i=2,dime(1)-1
         if (global_NP(i,j,k)==.false.) then
           if (global_NP(i+1,j,k).or.global_NP(i-1,j,k).or.global_NP(i,j+1,k).or.global_NP(i,j-1,k).or.global_NP(i,j,k+1).or.global_NP(i,j,k-1)) then
              global_boundary_NP(i,j,k) = .true.
              cont_b = cont_b + 1
           endif
         endif
      enddo
   enddo
enddo

! convert global to local
cont_b = 0
cont_np = 0

! select boundary_NP node for each NP 
do k=1,MaxZ
   do j=1,dime(2)
      do i=1,dime(1)
         kk = k + SumZproc(ID)
         if (kk<=dime(3)) then
            if (global_boundary_NP(i,j,kk).or.global_NP(i,j,kk)) then
               Do ii=1,number_NP !loop for number of NPs
                  xs = NP_position(:,ii) 

                  Lx = i-xs(1)
                  Ly = j-xs(2)
                  Lz = kk-xs(3)
       
                  checkx = (Lx)/real((R+0.5));
                  checkx = checkx*checkx;
                  checky = (Ly)/real((R+0.5));
                  checky = checky*checky;
                  checkz = (Lz)/real((R+0.5));
                  checkz = checkz*checkz;
                  check2 = checkx+checky+checkz;
           
                  if (check2<= 1.d0) then
                     if (global_boundary_NP(i,j,kk)) then
                        cont_b = cont_b + 1
                        temp_po(1,cont_b) = i
                        temp_po(2,cont_b) = j
                        temp_po(3,cont_b) = kk
                        temp_po(4,cont_b) = ii !save the center of NP
                        boundary_NP(i,j,k) = .true.
                     elseif (global_NP(i,j,kk)) then
                        cont_np = cont_np + 1
                        NP(i,j,k) = .true.
                     endif
                  endif
               enddo
            endif   
         endif
      enddo
   enddo
enddo

NPnodes = cont_np 
NPbnodes = cont_b
End subroutine grid_NP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine find_normal()
use type
use var
use parallel
implicit none
integer :: i, j, k, ii, jj,indx, label
real(dp) :: theta, phi, r, x,y,z
if (allocated(normal_NP)) deallocate(normal_NP)
if (allocated(tangent_NP)) deallocate(tangent_NP)
allocate(normal_NP(3,NPbnodes),tangent_NP(3,NPbnodes))

indx = 1
Do jj=1, NPbnodes
   i = temp_po(1,jj); j = temp_po(2,jj); k = temp_po(3,jj); ii = temp_po(4,jj)
   label = 1 + (i-1) + (j-1)*dime(1) + (k - SumZproc(ID) -1)*dime(1)*dime(2)
   x = real(i - NP_position(1,ii)); y = real(j - NP_position(2,ii)); z = real(k - NP_position(3,ii))
   ! appling periodic boundary
   if (x>dime(1)/2) then
      x=-dime(1)+x
   elseif (x<-dime(1)/2) then
      x=dime(1)+x
   endif
   if (y>dime(2)/2) then
      y=-dime(2)+y
   elseif (y<-dime(2)/2) then
      y=dime(2)+y
   endif
   if (z>dime(3)/2) then
      z=-dime(3)+z
   elseif (z<-dime(3)/2) then
      z=dime(3)+z
   endif

   r = sqrt(x*x + y*y + z*z)
   theta = acos(z/r)
   if (x/=0) then
      phi = atan(y/x)
   else
      phi = Pi /2.d0
   endif
   if (x<0.d0.and.y<0.d0) then
      phi = phi + Pi
   elseif (x<0.d0.and.y>0.d0) then 
      phi = phi + Pi
   elseif (x>0.d0.and.y<0.d0) then
      phi = phi + 2*Pi
   endif
   if (phi==0.d0.and.x<0) then
      theta = Pi - theta
   endif
   if (phi==Pi/2.d0.and.y<0) then
      theta = Pi - theta
   endif
   normal_NP(1,indx) = sin(theta)*cos(phi)
   normal_NP(2,indx) = sin(theta)*sin(phi)
   normal_NP(3,indx) = cos(theta)
   tangent_NP(1,indx) = cos(theta)*cos(phi)
   tangent_NP(2,indx) = cos(theta)*sin(phi)
   tangent_NP(3,indx) = -sin(theta)
   if (phi==0.d0.and.x<0) then
      normal_NP(:,indx) = -normal_NP(:,indx)
   endif
   if (phi==Pi/2.d0.and.y<0) then
      normal_NP(:,indx) = -normal_NP(:,indx)
   endif
   indx = indx + 1
enddo
endsubroutine find_normal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine anchoring_q()
use type
use var
use parallel
implicit none
include 'mpif.h'
integer :: i, j, k, x, ii, jj, label, indx
real(dp), dimension(3) :: dr, normal, eigenvalues, r
real(dp) :: s, rand(3)
if (allocated(qboundary_NP)) deallocate(qboundary_NP)
allocate(qboundary_NP(6,NPbnodes))
indx = 1

Do jj=1, NPbnodes
   
   if (NPtype=="xx") normal = tangent_NP(:,indx)
   if (NPtype=="zz") normal = normal_NP(:,indx)
   
   if (itype==3) then
   call RANDOM_NUMBER( rand )
   rand = rand - 0.5d0
   dr = sqrt( dot_product( rand, rand ) )
   rand = rand/dr
   normal = rand
   endif

   call director_tensor( Sbulk, normal, qboundary_NP(1:5,indx) )
   call eigen( qboundary_NP(:,indx), eigenvalues, r )

   s = maxval( eigenvalues )*3.d0/2.d0
   if ( abs(s - Sbulk) > 1.d-10 ) then 
   print *, '--> Problem with anchoring 1', i, r
   stop
   end if
   s = sqrt( dot_product( r,r ) )
   if ( abs(s - 1.d0) > 1.d-10 ) then 
   print *, '--> Problem with anchoring 2', i, r
   stop
   end if
   
   qboundary_NP(6,indx) = - qboundary_NP(1,indx) - qboundary_NP(4,indx)
   i = temp_po(1,jj); j = temp_po(2,jj); k = temp_po(3,jj)
   label = 1 + (i-1) + (j-1)*dime(1) + (k - SumZproc(ID) - 1)*dime(1)*dime(2)
   q(:,label) = qboundary_NP(:,indx)

   if (itype==3) then
      normal = normal_NP(:,indx)
      call director_tensor( Sinit, normal, qboundary_NP(1:5,indx) )
      qboundary_NP(6,indx) = - qboundary_NP(1,indx) - qboundary_NP(4,indx)
   endif


   indx = indx+1
enddo

End Subroutine anchoring_q 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine find_boundary_neighbors()
use var
use type
implicit none
integer :: i,j,k, indx, label
integer :: deltax,deltay,deltaz

if (allocated(local_boundary_neighbors_NP)) deallocate(local_boundary_neighbors_NP)
allocate(local_boundary_neighbors_NP(6,NPbnodes))
local_boundary_neighbors_NP = 0
indx = 1

do k=1,MaxZ
   do j=1,dime(2)
      do i=1,dime(1)
         if (boundary_NP(i,j,k)==.true.) then
            if (normal_NP(1,indx)>=0) deltax=1
            if (normal_NP(1,indx)< 0) deltax=-1
            if (normal_NP(2,indx)>=0) deltay=1
            if (normal_NP(2,indx)< 0) deltay=-1
            if (normal_NP(3,indx)>=0) deltaz=1
            if (normal_NP(3,indx)< 0) deltaz=-1
            local_boundary_neighbors_NP(1,indx) = 1 + (i-1+deltax) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
            local_boundary_neighbors_NP(2,indx) = 1 + (i-1+2*deltax) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
            local_boundary_neighbors_NP(3,indx) = 1 + (i-1) + (j-1+deltay)*dime(1) + (k-1)*dime(1)*dime(2)
            local_boundary_neighbors_NP(4,indx) = 1 + (i-1) + (j-1+2*deltay)*dime(1) + (k-1)*dime(1)*dime(2)
            local_boundary_neighbors_NP(5,indx) = 1 + (i-1) + (j-1)*dime(1) + (k-1+deltaz)*dime(1)*dime(2)
            local_boundary_neighbors_NP(6,indx) = 1 + (i-1) + (j-1)*dime(1) + (k-1+2*deltaz)*dime(1)*dime(2)
            indx = indx + 1
         endif
     enddo
   enddo
enddo
end subroutine find_boundary_neighbors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
End Subroutine make_np_channel

