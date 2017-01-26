Subroutine make_np_droplet
use type
use var
use parallel
implicit none
integer :: label, indx, deltax, deltay, deltaz
integer, dimension(:,:,:), allocatable :: temp_po 

allocate(share(dime(1),dime(2),dime(3)))
allocate(bshare(dime(1),dime(2),dime(3)))
share = .false.
bshare =.false.
allocate (temp_po(3,int(4*radius**3),number_NP))

NPnodes = 0
NPbnodes = 0
temp_po = 0
NP = .false.
boundary_NP = .false.

call grid_NP()
call find_normal1()
if (restart == 0) call anchoring_q()
!call qmga(0,0)
call find_boundary_neighbors1()
call correct_neighbors()
deallocate(temp_po)
deallocate(share,bshare)
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine grid_NP()
use var
implicit none
integer(i4b) :: i,j,k, ii 
integer(i4b) :: cont_b = 0 , cont_np = 0, cont_sh = 0 ,R
real(dp) :: checkx, checky, checkz, check1, check2, check3
real(dp) :: Lx, Ly, Lz
real, dimension(3) :: xs

R = radius

Do ii=1,number_NP !loop for number of NPs
xs = NP_position(:,ii) 

! select all nodes located in the NPs
do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
         if(drop(i,j,k)) then

           Lx = i-xs(1)
           Ly = j-xs(2)
           Lz = k-xs(3)

           checkx = (Lx)/real((R-0.5));
           checkx = checkx*checkx;
           checky = (Ly)/real((R-0.5));
           checky = checky*checky;
           checkz = (Lz)/real((R-0.5));
           checkz = checkz*checkz;
           check2 = checkx+checky+checkz;
           
           if (check2<= 1.d0) then
              NP(i,j,k) = .true.
              drop(i,j,k)=.false.
           endif
        endif   
      enddo
   enddo
enddo

enddo ! number of NP

! select all nodes which are neighbor of NP and choose them as boundary nodes
do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
         if(drop(i,j,k)) then

           if (NP(i+1,j,k).or.NP(i-1,j,k).or.NP(i,j+1,k).or.NP(i,j-1,k).or.NP(i,j,k+1).or.NP(i,j,k-1)) then
              drop(i,j,k) = .false.
              boundary_NP(i,j,k) = .true.
           endif
        endif   
      enddo
   enddo
enddo

! delete boundary nodes located in the NP
do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
         label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2) 
         if (boundary(i,j,k)) then
            indx = lb2index(label)
            deltax = 0 ; deltay = 0 ; deltaz =0
            if (normal0(1,indx)> 0) deltax=1
            if (normal0(1,indx)< 0) deltax=-1
            if (normal0(2,indx)> 0) deltay=1
            if (normal0(2,indx)< 0) deltay=-1
            if (normal0(3,indx)> 0) deltaz=1
            if (normal0(3,indx)< 0) deltaz=-1
            
            if (NP(i+deltax,j,k).or.NP(i,j+deltay,k).or.NP(i,j,k+deltaz)) then
               share(i,j,k) = .true.
               boundary(i,j,k) = .false.
               bnodes = bnodes - 1
            endif
         endif
      enddo
   enddo
enddo

Do ii=1,number_NP !loop for number of NPs
cont_b = 0
cont_np = 0
xs = NP_position(:,ii) 

! select boundary_NP node for each NP 
do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
         if(boundary_NP(i,j,k).or.NP(i,j,k)) then

           Lx = i-xs(1)
           Ly = j-xs(2)
           Lz = k-xs(3)

           checkx = (Lx)/real((R+0.5));
           checkx = checkx*checkx;
           checky = (Ly)/real((R+0.5));
           checky = checky*checky;
           checkz = (Lz)/real((R+0.5));
           checkz = checkz*checkz;
           check2 = checkx+checky+checkz;
           
           if (check2<= 1.d0) then
              if (boundary_NP(i,j,k)) then
                 cont_b = cont_b + 1
                 temp_po(1,cont_b,ii) = i
                 temp_po(2,cont_b,ii) = j
                 temp_po(3,cont_b,ii) = k
              elseif (NP(i,j,k)) then
                 cont_np = cont_np + 1
              endif
           endif
        endif   
      enddo
   enddo
enddo

Print *, ii, cont_np, cont_b
NPnodes(ii) = cont_np 
NPbnodes(ii) = cont_b
enddo

Print *, sum(NPnodes), sum(NPbnodes)

  do k=1,dime(3)
     do j=1,dime(2)
        do i=1,dime(1)

           if (boundary(i,j,k)) write(876,*) "C",i,j,k
           if (boundary_NP(i,j,k)) write(876,*) "O",i,j,k
           if (NP(i,j,k)) write(876,*) "N",i,j,k
           if (share(i,j,k)) write(876,*) "F",i,j,k
           if (bshare(i,j,k)) write(876,*) "K",i,j,k
        enddo
    enddo
 enddo


!write(*,*)sum(shnodes),sum(bshnodes)
End subroutine grid_NP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine find_normal1()
use type
use var 
implicit none
integer :: i, j, k, ii, jj,indx, label
real(dp) :: theta, phi, r, x,y,z
if (allocated(normal1)) deallocate(normal1)
if (allocated(tangent1)) deallocate(tangent1)
allocate(normal1(3,sum(NPbnodes)),tangent1(3,sum(NPbnodes)))
indx = 1
DO ii=1,number_NP
   Do jj=1, NPbnodes(ii)
      i = temp_po(1,jj,ii); j = temp_po(2,jj,ii); k = temp_po(3,jj,ii)
      label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
      lb2index(label) = indx
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
      normal1(1,indx) = sin(theta)*cos(phi)
      normal1(2,indx) = sin(theta)*sin(phi)
      normal1(3,indx) = cos(theta)
      tangent1(1,indx) = cos(theta)*cos(phi)
      tangent1(2,indx) = cos(theta)*sin(phi)
      tangent1(3,indx) = -sin(theta)
      if (phi==0.d0.and.x<0) then
         normal1(:,indx) = -normal1(:,indx)
      endif
      if (phi==Pi/2.d0.and.y<0) then
         normal1(:,indx) = -normal1(:,indx)
      endif
      indx = indx + 1
   enddo
enddo
endsubroutine find_normal1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine anchoring_q()
use type
use var
implicit none
integer :: i, j, k, x, ii, jj, label, indx
real(dp), dimension(3) :: dr, normal, eigenvalues, r
real(dp) :: s, rand(3)
if (allocated(qb1)) deallocate(qb1)
allocate(qb1(6,sum(NPbnodes)))
indx = 1
Do ii=1,number_NP
   Do jj=1, NPbnodes(ii)
      
      if (NPtype=="xx") normal = tangent1(:,indx)
      if (NPtype=="zz") normal = normal1(:,indx)
      
      if (itype==3) then
      call RANDOM_NUMBER( rand )
      rand = rand - 0.5d0
      dr = sqrt( dot_product( rand, rand ) )
      rand = rand/dr
      normal = rand
      endif

      call director_tensor( Sinit, normal, qb1(1:5,indx) )
      call eigen( qb1(:,indx), eigenvalues, r )
   
      s = maxval( eigenvalues )*3.d0/2.d0
      if ( abs(s - Sinit) > 1.d-10 ) then 
      print *, '--> Problem with anchoring 1', i, r
      stop
      end if
      s = sqrt( dot_product( r,r ) )
      if ( abs(s - 1.d0) > 1.d-10 ) then 
      print *, '--> Problem with anchoring 2', i, r
      stop
      end if
      
      qb1(6,indx) = - qb1(1,indx) - qb1(4,indx)
      i = temp_po(1,jj,ii); j = temp_po(2,jj,ii); k = temp_po(3,jj,ii)
      label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
      q(:,label) = qb1(:,indx)
      indx = indx+1
   enddo
enddo
End Subroutine anchoring_q 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine find_boundary_neighbors1()
use var
use type
implicit none
integer :: i,j,k, indx, label
integer :: deltax,deltay,deltaz
if (allocated(boundary_neighbors1)) deallocate(boundary_neighbors1)
allocate(boundary_neighbors1(6,sum(NPbnodes)))
boundary_neighbors1 = 0
indx = 1

do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
         label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2) 
         if (boundary_NP(i,j,k)==.true.) then
            indx = lb2index(label)
            if (normal1(1,indx)>=0) deltax=1
            if (normal1(1,indx)< 0) deltax=-1
            if (normal1(2,indx)>=0) deltay=1
            if (normal1(2,indx)< 0) deltay=-1
            if (normal1(3,indx)>=0) deltaz=1
            if (normal1(3,indx)< 0) deltaz=-1
            boundary_neighbors1(1,indx) = 1 + (i-1+deltax) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
            boundary_neighbors1(2,indx) = 1 + (i-1+2*deltax) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
            boundary_neighbors1(3,indx) = 1 + (i-1) + (j-1+deltay)*dime(1) + (k-1)*dime(1)*dime(2)
            boundary_neighbors1(4,indx) = 1 + (i-1) + (j-1+2*deltay)*dime(1) + (k-1)*dime(1)*dime(2)
            boundary_neighbors1(5,indx) = 1 + (i-1) + (j-1)*dime(1) + (k-1+deltaz)*dime(1)*dime(2)
            boundary_neighbors1(6,indx) = 1 + (i-1) + (j-1)*dime(1) + (k-1+2*deltaz)*dime(1)*dime(2)
         endif
     enddo
   enddo
enddo
end subroutine find_boundary_neighbors1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine correct_neighbors()
use var
use type
implicit none
integer :: i,j,k, indx, label, nfix, nmobile
integer :: deltax,deltay,deltaz
logical,dimension(:,:,:),allocatable :: fix, mobile
real(dp) :: rr
logical :: flag
allocate(fix(dime(1),dime(2),dime(3)))
allocate(mobile(dime(1),dime(2),dime(3)))

fix = .false.
mobile = .false.
nfix = 0
nmobile = 0 

do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)

         if (boundary(i,j,k).or.drop(i,j,k).or.boundary_NP(i,j,k)) then
            mobile(i,j,k) = .true.
            nmobile = nmobile + 1
         else
            fix(i,j,k) = .true.
            nfix = nfix + 1
         endif
      enddo
   enddo
enddo

! select boundary nodes don't have at least one neighbor 
do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
         label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2) 
         flag = .false.
         if (boundary(i,j,k)==.true.) then
            indx = lb2index(label)
            deltax = 0 ; deltay = 0 ; deltaz =0
            if (normal0(1,indx)> 0) deltax=1
            if (normal0(1,indx)< 0) deltax=-1
            if (normal0(2,indx)> 0) deltay=1
            if (normal0(2,indx)< 0) deltay=-1
            if (normal0(3,indx)> 0) deltaz=1
            if (normal0(3,indx)< 0) deltaz=-1
            
            if (fix(i+deltax,j,k).or.fix(i,j+deltay,k).or.fix(i,j,k+deltaz)) then
               bshare(i,j,k) = .true.
            if (fix(i+deltax,j,k)) then 
               flag = .true.
               normal0(1,indx) = 0.d0
            endif

            if (fix(i,j+deltay,k)) then 
               flag = .true.
               normal0(2,indx) = 0.d0
            endif
            
            if (fix(i,j,k+deltaz)) then 
               flag = .true.
               normal0(3,indx) = 0.d0
            endif

            if (flag) then
               rr = sqrt(normal0(1,indx)*normal0(1,indx)+normal0(2,indx)*normal0(2,indx)+normal0(3,indx)*normal0(3,indx))
               normal0(:,indx) = normal0(:,indx) / rr
            endif
            endif
         endif
      enddo
   enddo
enddo

! select boundary nodes which have first neighbor but not second one
! get first accurate derivative
do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
         label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2) 
         if ((boundary(i,j,k)==.true.).and.(bshare(i,j,k)==.false.)) then
            indx = lb2index(label)
            deltax = 0 ; deltay = 0 ; deltaz =0
            if (normal0(1,indx)> 0) deltax=1
            if (normal0(1,indx)< 0) deltax=-1
            if (normal0(2,indx)> 0) deltay=1
            if (normal0(2,indx)< 0) deltay=-1
            if (normal0(3,indx)> 0) deltaz=1
            if (normal0(3,indx)< 0) deltaz=-1
            
            if (fix(i+2*deltax,j,k).or.fix(i,j+2*deltay,k).or.fix(i,j,k+2*deltaz)) then
              bshare(i,j,k) = .true.
            endif

            if (fix(i+2*deltax,j,k)) then 
               boundary_neighbors0(2,indx) = label
               normal0(1,indx) = (1.d0/2.d0)*normal0(1,indx)
            endif

            if (fix(i,j+2*deltay,k)) then 
               boundary_neighbors0(4,indx) = label
               normal0(2,indx) = (1.d0/2.d0)*normal0(2,indx)
            endif
            
            if (fix(i,j,k+2*deltaz)) then 
               boundary_neighbors0(6,indx) = label 
               normal0(3,indx) = (1.d0/2.d0)*normal0(3,indx)
            endif
         endif
     enddo
   enddo
enddo

! select NP boundary nodes don't have at least one neighbor 
do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
         label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2) 
         flag = .false.
         if (boundary_NP(i,j,k)==.true.) then
            indx = lb2index(label)
            deltax = 0 ; deltay = 0 ; deltaz =0
            if (normal1(1,indx)> 0) deltax=1
            if (normal1(1,indx)< 0) deltax=-1
            if (normal1(2,indx)> 0) deltay=1
            if (normal1(2,indx)< 0) deltay=-1
            if (normal1(3,indx)> 0) deltaz=1
            if (normal1(3,indx)< 0) deltaz=-1
            
            if (fix(i+deltax,j,k).or.fix(i,j+deltay,k).or.fix(i,j,k+deltaz)) then
               bshare(i,j,k) = .true.
            if (fix(i+deltax,j,k)) then 
               flag = .true.
               normal1(1,indx) = 0.d0
            endif

            if (fix(i,j+deltay,k)) then 
               flag = .true.
               normal1(2,indx) = 0.d0
            endif
            
            if (fix(i,j,k+deltaz)) then 
               flag = .true.
               normal1(3,indx) = 0.d0
            endif

            if (flag) then
               rr = sqrt(normal1(1,indx)*normal1(1,indx)+normal1(2,indx)*normal1(2,indx)+normal1(3,indx)*normal1(3,indx))
               normal1(:,indx) = normal1(:,indx) / rr
            endif
            endif
         endif
      enddo
   enddo
enddo

! select NP boundary nodes which have first neighbor but not second one
! get first accurate derivative
do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
         label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2) 
         if ((boundary_NP(i,j,k)==.true.).and.(bshare(i,j,k)==.false.)) then
            indx = lb2index(label)
            deltax = 0 ; deltay = 0 ; deltaz =0
            if (normal1(1,indx)> 0) deltax=1
            if (normal1(1,indx)< 0) deltax=-1
            if (normal1(2,indx)> 0) deltay=1
            if (normal1(2,indx)< 0) deltay=-1
            if (normal1(3,indx)> 0) deltaz=1
            if (normal1(3,indx)< 0) deltaz=-1
            
            if (fix(i+2*deltax,j,k).or.fix(i,j+2*deltay,k).or.fix(i,j,k+2*deltaz)) then
               bshare(i,j,k) = .true.
            endif

            if (fix(i+2*deltax,j,k)) then 
               boundary_neighbors1(2,indx) = label
               normal1(1,indx) = (1.d0/2.d0)*normal1(1,indx)
            endif

            if (fix(i,j+2*deltay,k)) then 
               boundary_neighbors1(4,indx) = label
               normal1(2,indx) = (1.d0/2.d0)*normal1(2,indx)
            endif
            
            if (fix(i,j,k+2*deltaz)) then 
               boundary_neighbors1(6,indx) = label 
               normal1(3,indx) = (1.d0/2.d0)*normal1(3,indx)
            endif
         endif
     enddo
   enddo
enddo


 do k=1,dime(3)
     do j=1,dime(2)
        do i=1,dime(1)

           if (boundary(i,j,k)) write(875,*) "C",i,j,k
           if (boundary_NP(i,j,k)) write(875,*) "O",i,j,k
           if (NP(i,j,k)) write(875,*) "N",i,j,k
           if (share(i,j,k)) write(875,*) "F",i,j,k
           if (bshare(i,j,k)) write(875,*) "K",i,j,k
        enddo
    enddo
 enddo

end subroutine correct_neighbors

End Subroutine make_np_droplet

