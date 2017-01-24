subroutine initial()
use type
use var
implicit none

bnodes=0 
boundary_master = .false.
call grid()
call find_neighbors()
call find_normal()
call find_boundary_neighbors()
call q_initial()

call anchoring()

!call qmga(0,0)
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine grid
use type
use var
implicit none
integer(i4b) :: i,j,k,label
integer(i4b) :: cont_b = 0 , cont_d = 0, R
real(dp) :: checkx, checky, checkz, check1, check2, check3
real(dp) :: Lx, Ly, Lz
!********************************************
!loop for surface and inside nodes
center = int(dime/2) + 1
!center = (real(dime)-1.d0)/2.d0  !Ye rep
R = center(1)-2
!R = dime(1)/2.d0-2.d0 !Ye rep
cont_b = 0
cont_d = 0
boundary_master = .false.

do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
         if (channel) then
            if (k==1.or.k==dime(3)) then
                boundary_master(i,j,k) = .true.
                cont_b = cont_b + 1
            endif

         elseif (droplet) then
            Lx = i-center(1)
            Ly = j-center(2)
            Lz = k-center(3)

            checkx = (Lx)/real((R+0.5))
            checkx = checkx*checkx
            checky = (Ly)/real((R+0.5))
            checky = checky*checky
            checkz = (Lz)/real((R+0.5))
            checkz = checkz*checkz
            check1 = checkx+checky+checkz
           
            if (check1<= 1.d0) then
               drop(i,j,k) = .true.
               cont_d = cont_d + 1
            endif

         endif  
      enddo
   enddo
enddo

if (droplet) then
  do k=2,dime(3)-1
     do j=2,dime(2)-1
        do i=2,dime(1)-1

           label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)

           if (drop(i,j,k).and.(drop(i+1,j,k)==.false..or.drop(i-1,j,k)==.false. &
                            .or.drop(i,j+1,k)==.false..or.drop(i,j-1,k)==.false. &
                            .or.drop(i,j,k+1)==.false..or.drop(i,j,k-1)==.false.)) then
              boundary_master(i,j,k) = .true.
              cont_d = cont_d - 1
              cont_b = cont_b + 1
           endif

        enddo
    enddo
 enddo

do k=1,dime(3)
    do j=1,dime(2)
       do i=1,dime(1)

          if (boundary_master(i,j,k)) drop(i,j,k) = .false.

       enddo
   enddo
enddo

endif

bnodes = cont_b
dnodes = cont_d

end subroutine grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine find_neighbors
use type
use var
implicit none
integer(i4b) :: i, j, k, ii, kp 
integer :: x, y, left, right, label_l, label_r, label
integer, dimension(3) :: pos, pos_l, pos_r

neighbors_master = 0 
do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
         label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2) 
         pos(1) = i ; pos(2) = j ; pos(3) = k 
         Do x = 1 , 3
            left = pos(x) - 1
            right = pos(x) + 1
            if (x<=3) then ! Apply periodic boundary
               if (left<1) left = dime(x)
               if (right>dime(x)) right = 1
            endif
            Do y = 1 , 3
               if (y/=x) then
                  pos_l(y) = pos(y)
                  pos_r(y) = pos(y)
               else
                  pos_l(y) = left
                  pos_r(y) = right
               endif
            enddo
            label_l = 1 + (pos_l(1)-1) + (pos_l(2)-1)*dime(1) + (pos_l(3)-1)*dime(1)*dime(2)
            label_r = 1 + (pos_r(1)-1) + (pos_r(2)-1)*dime(1) + (pos_r(3)-1)*dime(1)*dime(2)
            ii = (x-1)*2 + 1
            neighbors_master(ii,label)   = label_l 
            neighbors_master(ii+1,label) = label_r 
         enddo
      enddo
   enddo
enddo
end subroutine find_neighbors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine find_normal()
use var
use type
use parallel
implicit none
integer :: i,j,k, indx, label
real(dp) :: x,y,z
real(dp) :: theta, phi, r
allocate(normal_master(3,bnodes),tangent_master(3,bnodes))
normal_master=0.d0 ; tangent_master=0.d0; indx = 1

if (channel)  then
   do k=1,dime(3)
      do j=1,dime(2)
         do i=1,dime(1)
            if (boundary_master(i,j,k)) then
               if (k==1) normal_master(3,indx) = 1.d0
               if (k==dime(3)) normal_master(3,indx) = -1.d0
               indx = indx + 1
            endif
         enddo
      enddo
   enddo
endif
if (droplet) then
   do k=1,dime(3)
      do j=1,dime(2)
         do i=1,dime(1)
            if (boundary_master(i,j,k)) then
               x = real(i - center(1)); y = real(j - center(2)); z = real(k - center(3))
               r = sqrt(x*x + y*y + z*z)
               theta = acos(z/r)
                  phi = atan2(y,x)
               if (phi==0.d0.and.x<0) then
                  theta = Pi - theta
               endif
               if (phi==Pi/2.d0.and.y<0) then
                  theta = Pi - theta
               endif
               normal_master(1,indx) = sin(theta)*cos(phi)
               normal_master(2,indx) = sin(theta)*sin(phi)
               normal_master(3,indx) = cos(theta)
               tangent_master(1,indx) = cos(theta)*cos(phi)
               tangent_master(2,indx) = cos(theta)*sin(phi)
               tangent_master(3,indx) = -sin(theta)
               if (phi==0.d0.and.x<0) then
                  normal_master(:,indx) = -normal_master(:,indx)
               endif
               if (phi==Pi/2.d0.and.y<0) then
                  normal_master(:,indx) = -normal_master(:,indx)
               endif
               indx = indx + 1
            endif
         enddo
      enddo
   enddo
   normal_master = - normal_master ! the unit vector of droplet inside 
endif

end subroutine find_normal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine find_boundary_neighbors()
use var
use type
implicit none
integer :: i,j,k, label, indx
integer :: deltax,deltay,deltaz
indx = 0

do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
         label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2) 
         if (channel) then
         if (k==1) then
            neighbors_master(5,label) = 1 + (i-1) + (j-1)*dime(1) + (k-1+1)*dime(1)*dime(2)
            neighbors_master(6,label) = 1 + (i-1) + (j-1)*dime(1) + (k-1+2)*dime(1)*dime(2)
         elseif (k==dime(3)) then
            neighbors_master(5,label) = 1 + (i-1) + (j-1)*dime(1) + (k-1-1)*dime(1)*dime(2)
            neighbors_master(6,label) = 1 + (i-1) + (j-1)*dime(1) + (k-1-2)*dime(1)*dime(2)
         endif
         elseif (droplet) then
            if (boundary_master(i,j,k)==.true.) then
               indx = indx + 1
               if (normal_master(1,indx)>=0) deltax=1
               if (normal_master(1,indx)< 0) deltax=-1
               if (normal_master(2,indx)>=0) deltay=1
               if (normal_master(2,indx)< 0) deltay=-1
               if (normal_master(3,indx)>=0) deltaz=1
               if (normal_master(3,indx)< 0) deltaz=-1
               neighbors_master(1,label) = 1 + (i-1+deltax) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
               neighbors_master(2,label) = 1 + (i-1+2*deltax) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
               neighbors_master(3,label) = 1 + (i-1) + (j-1+deltay)*dime(1) + (k-1)*dime(1)*dime(2)
               neighbors_master(4,label) = 1 + (i-1) + (j-1+2*deltay)*dime(1) + (k-1)*dime(1)*dime(2)
               neighbors_master(5,label) = 1 + (i-1) + (j-1)*dime(1) + (k-1+deltaz)*dime(1)*dime(2)
               neighbors_master(6,label) = 1 + (i-1) + (j-1)*dime(1) + (k-1+2*deltaz)*dime(1)*dime(2)

               if ((boundary_master(i+deltax,j,k)==.false..and.drop(i+deltax,j,k)==.false.)) neighbors_master(1,label)=label
               if ((boundary_master(i,j+deltay,k)==.false..and.drop(i,j+deltay,k)==.false.)) neighbors_master(3,label)=label
               if ((boundary_master(i,j,k+deltaz)==.false..and.drop(i,j,k+deltaz)==.false.)) neighbors_master(5,label)=label

               if (boundary_master(i+2*deltax,j,k)==.false..and.drop(i+2*deltax,j,k)==.false.) neighbors_master(2,label)=neighbors_master(1,label)
               if (boundary_master(i,j+2*deltay,k)==.false..and.drop(i,j+2*deltay,k)==.false.) neighbors_master(4,label)=neighbors_master(3,label)
               if (boundary_master(i,j,k+2*deltaz)==.false..and.drop(i,j,k+2*deltaz)==.false.) neighbors_master(6,label)=neighbors_master(5,label)

            endif
         endif
     enddo
   enddo
enddo

end subroutine find_boundary_neighbors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine q_initial
use type
use var
implicit none
integer(i4b) i,j,k, label,N,ADJX,ADJZ
real(dp) max_e, s, dr, rs,qchiralbpI,RR2
real(dp) :: x,y,z, omega,theta,phi,rr,alpha,z2,x2
real(dp), dimension(3):: r = 0.d0, rand = 0.d0, eigenvalues
real(dp), dimension(5) :: q_temp 

real(dp) :: thetax, thetay, thetaz
real(dp)  xxr(3), rot(3,3),xx(3)


!if (itype==6) then
!        if(channel)then
        open(60,file='angle.dat')
        read(60,*)thetax,rs,ADJX,ADJZ
        close(60)
!        endif

!        if(.not.(channel).and..not.(droplet))then
!        open(60,file='angle.dat')
!        read(60,*)thetax,rs,ADJX,ADJZ
!        close(60)
!        endif


!endif


do i=1,dime(1)
   do j=1,dime(2)
     do k=1,dime(3)

         label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
            
            call RANDOM_NUMBER( rand )
            rand = rand - 0.5d0
            dr = sqrt( dot_product( rand, rand ) )
            rand = rand/dr
            if (itype == 2) then
               rand = 0.0
               if (channel) rand(1) = 1.0
               if (droplet) rand(3) = 1.0
            endif

               x = real(i-real(dime(1))/2.d0); y = real(j-real(dime(2))/2.d0); z = real(k-real(dime(3))/2.d0)
!               x = real((i)-real(dime(1)-1)/2.d0); y = real((j)-real(dime(2)-1)/2.d0); z = real((k)-real(dime(3)-1)/2.d0)
                 rr = sqrt(x*x + y*y + z*z)
                 if (rr == 0.d0) then
                 theta=pi/2.d0 !0.d0
                 else
                 theta = acos(z/rr)
                 endif
                 phi = atan2(y,x)

                if(itype==6) then !radial
                 rand = 0.0
                 rand(1)=-sin(phi)
                 rand(2)=cos(phi)

                 elseif(itype==7) then !DSS
                 rand = 0.0
                 omega =qchiral*rr
                rand(1) =  -sin(omega)*sin(phi)+cos(omega)*cos(theta)*cos(phi)
                rand(2) =   sin(omega)*cos(phi)+cos(omega)*cos(theta)*sin(phi)
                rand(3) =  -cos(omega)*sin(theta)

                elseif(itype==8) then !RSS
                 rand = 0.0
                 omega =phi+qchiral*rr
                rand(1) =  -sin(omega)*sin(phi)+cos(omega)*cos(theta)*cos(phi)
                rand(2) =   sin(omega)*cos(phi)+cos(omega)*cos(theta)*sin(phi)
                rand(3) =  -cos(omega)*sin(theta)

                elseif(itype==9) then !Cholesteric
                 rand = 0.0
                rand(1)=cos(z*qchiral)
                rand(2)=sin(z*qchiral)
                rand(3)=0.d0

                elseif(itype==10) then !BS                
                rand = 0.0
                z2=z* sqrt(x*x + y*y)
                RR2 = real(center(1)-2)
                x2=real(RR2)*real(RR2)-z*z                
                theta=atan2(z2,x2) 
                  if (abs(z).ge.RR2)then
                  alpha=pi/2.d0 
                   else
                  alpha=real(RR2)*qchiral*sqrt(x*x + y*y)/sqrt(real(RR2)**2-z**2)
                  endif
                rand(1)=-sin(theta)*cos(alpha)*cos(phi)- sin(alpha)*sin(phi) 
                rand(2)=-sin(theta)*cos(alpha)*sin(phi)+ sin(alpha)*cos(phi)
                rand(3)=cos(theta)*cos(alpha)

                endif


            call director_tensor( Sinit, rand, q_temp )
            ! Checking all the values of the tensor  
            call eigen( q_temp, eigenvalues, r )   ! max_e is the maximun eigenvalue and r is the eigenvector
            s = maxval( eigenvalues )*3.d0/2.d0
            if ( abs(s - Sinit) > 1.d-10 ) then
               print *, '--> Problem with anchoring 1', i, r
               stop
            end if
            dr = sqrt( dot_product( r,r ) )
            if ( abs(dr - 1.d0) > 1.d-10 ) then 
               print *, '--> Problem with anchoring 2', i, r
               stop
            end if

            if (itype==4.or.itype==5.or.itype==6 .or. itype==11) then
!               x = real((i-ADJX)-real(dime(1)-1)/2.d0); y = real((j-ADJX)-real(dime(2)-1)/2.d0); z = real((k-ADJZ)-real(dime(3)-1)/2.d0)
               x = real((i-ADJX)-real(dime(1))/2.d0); y = real((j-ADJX)-real(dime(2))/2.d0); z = real((k-ADJZ)-real(dime(3))/2.d0)

               if (itype==4) then !blue phase I
!                  rs=0.68d0
!               x = real(i-real(dime(1))/2.d0); y = real(j-real(dime(2))/2.d0); z = real(k-real(dime(3))/2.d0)

                 q_temp(1)=0.2d0*(-sin(sqrt(2.d0)*qchiral*rs*y)*cos(sqrt(2.d0)*qchiral*x*rs)+&
                            -sin(sqrt(2.d0)*qchiral*rs*x)*cos(sqrt(2.d0)*qchiral*z*rs)+2.d0*sin(sqrt(2.d0)*qchiral*rs*z)*cos(sqrt(2.d0)*qchiral*y*rs))

                  q_temp(2)=0.2d0*(-sqrt(2.d0)*sin(sqrt(2.d0)*qchiral*rs*x)*sin(sqrt(2.d0)*qchiral*z*rs)+&
                            -sqrt(2.d0)*cos(sqrt(2.d0)*qchiral*rs*y)*cos(sqrt(2.d0)*qchiral*z*rs)+sin(sqrt(2.d0)*qchiral*rs*x)*cos(sqrt(2.d0)*qchiral*y*rs))

                  q_temp(3)=0.2d0*(-sqrt(2.d0)*sin(sqrt(2.d0)*qchiral*rs*z)*sin(sqrt(2.d0)*qchiral*y*rs)+&
                            -sqrt(2.d0)*cos(sqrt(2.d0)*qchiral*rs*x)*cos(sqrt(2.d0)*qchiral*y*rs)+sin(sqrt(2.d0)*qchiral*rs*z)*cos(sqrt(2.d0)*qchiral*x*rs))

                  q_temp(4)=0.2d0*(-sin(sqrt(2.d0)*qchiral*rs*z)*cos(sqrt(2.d0)*qchiral*y*rs)+&
                            -sin(sqrt(2.d0)*qchiral*rs*y)*cos(sqrt(2.d0)*qchiral*x*rs)+2.d0*sin(sqrt(2.d0)*qchiral*rs*x)*cos(sqrt(2.d0)*qchiral*z*rs))

                  q_temp(5)=0.2d0*(-sqrt(2.d0)*sin(sqrt(2.d0)*qchiral*rs*y)*sin(sqrt(2.d0)*qchiral*x*rs)+&
                            -sqrt(2.d0)*cos(sqrt(2.d0)*qchiral*rs*z)*cos(sqrt(2.d0)*qchiral*x*rs)+sin(sqrt(2.d0)*qchiral*rs*y)*cos(sqrt(2.d0)*qchiral*z*rs))

               endif

               if (itype==6 .or. itype==11) then !blue phase I
                    !thetay=thetax

                    rot(1,1) = 1.00
                    rot(1,2) = 0.0
                    rot(1,3) = 0.0
                    rot(2,1) = 0.0
                    rot(2,2) = cos(thetax/180*pi)
                    rot(2,3) = sin(thetax/180*pi)
                    rot(3,1) = 0.0
                    rot(3,2) = -sin(thetax/180*pi)
                    rot(3,3) = cos(thetax/180*pi)

                    xx(1) = x
                    xx(2) = y
                    xx(3) = z

                    xxr=matmul(rot,xx)

                    x = xxr(1)
                    y = xxr(2)
                    z = xxr(3)
   
                  if (itype==6 ) then

                    q_temp(1)=0.2d0*(-sin(sqrt(2.d0)*qchiral*rs*y)*cos(sqrt(2.d0)*qchiral*x*rs)+&
                         -sin(sqrt(2.d0)*qchiral*rs*x)*cos(sqrt(2.d0)*qchiral*z*rs)+2.d0*sin(sqrt(2.d0)*qchiral*rs*z)*cos(sqrt(2.d0)*qchiral*y*rs))

                    q_temp(2)=0.2d0*(-sqrt(2.d0)*sin(sqrt(2.d0)*qchiral*rs*x)*sin(sqrt(2.d0)*qchiral*z*rs)+&
                         -sqrt(2.d0)*cos(sqrt(2.d0)*qchiral*rs*y)*cos(sqrt(2.d0)*qchiral*z*rs)+sin(sqrt(2.d0)*qchiral*rs*x)*cos(sqrt(2.d0)*qchiral*y*rs))

                    q_temp(3)=0.2d0*(-sqrt(2.d0)*sin(sqrt(2.d0)*qchiral*rs*z)*sin(sqrt(2.d0)*qchiral*y*rs)+&
                         -sqrt(2.d0)*cos(sqrt(2.d0)*qchiral*rs*x)*cos(sqrt(2.d0)*qchiral*y*rs)+sin(sqrt(2.d0)*qchiral*rs*z)*cos(sqrt(2.d0)*qchiral*x*rs))
                    
                    q_temp(4)=0.2d0*(-sin(sqrt(2.d0)*qchiral*rs*z)*cos(sqrt(2.d0)*qchiral*y*rs)+&
                         -sin(sqrt(2.d0)*qchiral*rs*y)*cos(sqrt(2.d0)*qchiral*x*rs)+2.d0*sin(sqrt(2.d0)*qchiral*rs*x)*cos(sqrt(2.d0)*qchiral*z*rs))

                    q_temp(5)=0.2d0*(-sqrt(2.d0)*sin(sqrt(2.d0)*qchiral*rs*y)*sin(sqrt(2.d0)*qchiral*x*rs)+&
                            -sqrt(2.d0)*cos(sqrt(2.d0)*qchiral*rs*z)*cos(sqrt(2.d0)*qchiral*x*rs)+sin(sqrt(2.d0)*qchiral*rs*y)*cos(sqrt(2.d0)*qchiral*z*rs))
                  endif

                 if (itype==11 ) then
                 q_temp(1)=0.2d0*(cos(2.d0*qchiral*rs*z)-cos(2.d0*qchiral*y*rs) )
                 q_temp(2)=0.2d0*(sin(2.d0*qchiral*rs*z) )
                 q_temp(3)=0.2d0*(sin(2.d0*qchiral*rs*y) )
                 q_temp(4)=0.2d0*(cos(2.d0*qchiral*rs*x)-cos(2.d0*qchiral*z*rs) )
                 q_temp(5)=0.2d0*(sin(2.d0*qchiral*rs*x) )
                 endif

              endif


               if(itype==5) then !blue phase II
!                 rs=0.86d0
              ! x = real(i - center(1)); y = real(j - center(2)); z = real(k - center(3))
               x = real(i - Int(dime(1)/2)); y = real(j - Int(dime(2)/2)); z = real(k - Int(dime(3)/2))
!               x = real((i-ADJX)-real(dime(1))/2.d0); y = real((j-ADJX)-real(dime(2))/2.d0); z = real((k-ADJZ)-real(dime(3))/2.d0)
                 q_temp(1)=0.2d0*(cos(2.d0*qchiral*rs*z)-cos(2.d0*qchiral*y*rs) )
                 q_temp(2)=0.2d0*(sin(2.d0*qchiral*rs*z) )
                 q_temp(3)=0.2d0*(sin(2.d0*qchiral*rs*y) )
                 q_temp(4)=0.2d0*(cos(2.d0*qchiral*rs*x)-cos(2.d0*qchiral*z*rs) )
                 q_temp(5)=0.2d0*(sin(2.d0*qchiral*rs*x) )
               endif  
            endif

            q_master(:,label) = q_temp
      
      enddo
   enddo
enddo


end subroutine q_initial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine anchoring
use type
use var
implicit none
integer(i4b) i,j,k, label, indx 
real(dp) max_e, s, dr
real(dp), dimension(3):: r, planar = 0.d0, eigenvalues, rand 
real(dp), dimension(5) :: q_temp 

Allocate(qboundary_master(6,bnodes))

indx = 1
do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)

         if (boundary_master(i,j,k) == .true.) then
            
            !type of channel
            planar = 0
            if (channel) then
               if (k==1) then
                  if (Stype/="zz") planar(1) = 1
                  if (Stype=="zz") planar(3) = 1
               elseif (k==dime(3)) then
                  if (Stype=="xx") planar(1) = 1
                  if (Stype=="xy") planar(2) = 1
                  if (Stype=="xz") planar(3) = 1
                  if (Stype=="zz") planar(3) = 1
               endif
            elseif (droplet) then
               if (Stype=="xx") planar = tangent_master(:,indx)
               if (Stype=="zz") planar = normal_master(:,indx)
            endif

            call director_tensor( Sbulk, planar, q_temp)
            ! Checking all the values of the tensor  
            call eigen( q_temp, eigenvalues, r )   ! max_e is the maximun eigenvalue and r is the eigenvector
            s = maxval( eigenvalues )*3.d0/2.d0
            if ( abs(s - Sbulk) > 1.d-10 ) then
               print *, '--> Problem with anchoring 1', i, r
               stop
            end if
            dr = sqrt( dot_product( r,r ) )
            if ( abs(dr - 1.d0) > 1.d-10 ) then 
               print *, '--> Problem with anchoring 2', i, r
               stop
            end if

            qboundary_master(1:5,indx) = q_temp
            qboundary_master(6,indx) = - qboundary_master(1,indx) - qboundary_master(4,indx)
            
            label = 1 + (i-1) + (j-1)*dime(1) + (k - 1)*dime(1)*dime(2)
!            if (itype/=3) q_master(:,label) = q_temp !qboundary(:,indx)
            indx = indx+1
         endif
     
      enddo
   enddo
enddo
  

end subroutine anchoring

end subroutine initial
