subroutine ginzburg_landau(dt_c)
use type
use var
use parallel
implicit none
include 'mpif.h'
integer :: x, label_l, label_r, indx, indx_S, indx_NP, ii
real(dp),dimension(6) :: q_tmp, q_tmp_new, dotd, delta, LdG, dQz, dQx, dQy, dQ, Qtrow,Qprow, Eltot, ChiralE,chir_x,chir_y,chir_z, chir
real(dp),dimension(18) :: grad1
real(dp) :: doubledot, doubledot2, nuQnu
real(dp) :: ElX, Ely, Elz, dt_c

delta(1) = 1.0d0 / 3.0d0 ; delta(4) =  1.0d0/3.0d0 ; delta(6) = 1.0d0/3.d0
delta(2) = 0.d0          ; delta(2) = 0.d0         ; delta(2) = 0.d0
indx_S = 0
indx_NP = 0
ChiralE = 0.d0

do indx=1,length
   q_tmp = dble(q(:,indx))
   if (bulk(indx) == .true.) then
      doubledot = dble(q(1,indx))**2 + 2.d0*dble(q(2,indx))**2 + 2.d0*dble(q(3,indx))**2 + dble(q(4,indx))**2 + 2.d0*dble(q(5,indx))**2 + dble(q(6,indx))**2
      call dotq (q_tmp,dotd)
      LdG = (1.d0 - UldG/3.d0)*q_tmp - ULdG * dotd + Uldg*doubledot * (q_tmp + delta )
     
      do x=1,6 
         label_l = neighbors(1,indx)
         label_r = neighbors(2,indx)
         Elx = ddf(1) *(dble(q(x,label_l)) - 2.d0*q_tmp(x) + dble(q(x,label_r)))
         label_l = neighbors(3,indx)
         label_r = neighbors(4,indx)
         Ely = ddf(2) *(dble(q(x,label_l)) - 2.d0*q_tmp(x) + dble(q(x,label_r)))
         label_l = neighbors(5,indx)
         label_r = neighbors(6,indx)
         Elz = ddf(3) *(dble(q(x,label_l)) - 2.d0*q_tmp(x) + dble(q(x,label_r)))

         Eltot(x) = Elx + Ely + Elz
      enddo

      if (chiral>0) then
         Do x = 1 , 3 
            ii = (x-1)*2 + 1
            label_l = neighbors(ii,indx)
            label_r = neighbors(ii+1,indx)
            grad1(1+(x-1)*6:x*6) =dble(q(:,label_r)-q(:,label_l))*df2(x)
         enddo
         
         ChiralE(1) = 2.0d0*qchiral*(grad1(9)-grad1(14))
         ChiralE(2) = qchiral*(grad1(11)-grad1(16)-grad1(3)+grad1(13))
         ChiralE(3) = qchiral*(grad1(12)-grad1(17)+grad1(2)-grad1(7))
         ChiralE(4) = 2.0d0*qchiral*(grad1(14)-grad1(5))
         ChiralE(5) = qchiral*(grad1(15)-grad1(6)+grad1(4)-grad1(8))
         ChiralE(6) = 2.0d0*qchiral*(grad1(5)-grad1(9))
      endif

      Do x = 1,6
         q_tmp_new (x)=  q_tmp (x) - dt_c * ( LdG(x) - Eltot(x) + 2.0d0*ChiralE(x) ) 
      enddo
    
   elseif (boundary(indx) == .true..and.infinite==.false.) then
      indx_S = indx_S + 1
      dQx=0.d0; dQy=0.d0; dQz=0.d0

      if (droplet) then
         label_l = neighbors(1,indx)
         label_r = neighbors(2,indx)
         dQx = (-3.0d0*q_tmp + 4.d0*dble(q(:,label_l)) - dble(q(:,label_r))) * df2(1)
         label_l = neighbors(3,indx)
         label_r = neighbors(4,indx)
         dQy = (-3.0d0*q_tmp + 4.d0*dble(q(:,label_l)) - dble(q(:,label_r))) * df2(2)
         label_l = neighbors(5,indx)
         label_r = neighbors(6,indx)
         dQz = (-3.0d0*q_tmp + 4.d0*dble(q(:,label_l)) - dble(q(:,label_r))) * df2(3)
         dQ = dQx*abs(normal(1,indx_S)) + dQy*abs(normal(2,indx_S)) + dQz*abs(normal(3,indx_S))
      elseif (channel) then
!         label_l = neighbors(5,indx_S)
!         label_r = neighbors(6,indx_S)
         label_l = neighbors(5,indx)
         label_r = neighbors(6,indx)
         dQz = (-3.0d0*q_tmp + 4.d0*dble(q(:,label_l)) - dble(q(:,label_r))) * df2(3)
         dQ = dQz*abs(normal(3,indx_S))
      endif
      
      
      if (chiral>0) then
          chir_x(1)=  0.d0
          chir_y(1)=  2.d0*q_tmp(3)
          chir_z(1)= -2.d0*q_tmp(2)

          chir_x(2)= -1.d0*q_tmp(3)
          chir_y(2)=  1.d0*q_tmp(5)
          chir_z(2)=  1.d0*(q_tmp(1)-q_tmp(4))

          chir_x(3)=  1.d0*q_tmp(2)
          chir_y(3)=  1.d0*(q_tmp(6)-q_tmp(1))
          chir_z(3)= -1.d0*q_tmp(5)

          chir_x(4)= -2.d0*q_tmp(5)
          chir_y(4)=  0.d0
          chir_z(4)=  2.d0*q_tmp(2)

          chir_x(5)=  1.d0*(q_tmp(4)-q_tmp(6))
          chir_y(5)= -1.d0*q_tmp(2)
          chir_z(5)=  1.d0*q_tmp(3)

          chir_x(6)=  2.d0*q_tmp(5)
          chir_y(6)= -2.d0*q_tmp(3)
          chir_z(6)=  0.d0

          chir = qchiral*(chir_x*(normal(1,indx_S))+chir_y*(normal(2,indx_S))+chir_z*(normal(3,indx_S)))
          dQ = dQ - chir
      endif 
      
      if (nondegenerate) then
         q_tmp_new =  q_tmp - dt_c * ( - dQ + gama(1)*(q_tmp - qboundary(:,indx_S)))
      elseif (degenerate) then
         call surface(q_tmp,normal(:,indx_S),Qtrow,Qprow,nuQnu)
         q_tmp_new  =  q_tmp - dt_c * ( - dQ + 2.d0*gama(1)*((Qtrow-Qprow)-delta*nuQnu))
      endif
   
!   elseif (boundary_NP(i,j,k) == .true..and.NP_infinite==.false.) then
!      indx_NP = indx_NP + 1
!      label_l = local_boundary_neighbors_NP(1,indx_NP)
!      label_r = local_boundary_neighbors_NP(2,indx_NP)
!      dQx = (-3.0d0*q_tmp + 4.d0*dble(q(:,label_l)) - dble(q(:,label_r))) * df2(1)
!      label_l = local_boundary_neighbors_NP(3,indx_NP)
!      label_r = local_boundary_neighbors_NP(4,indx_NP)
!      dQy = (-3.0d0*q_tmp + 4.d0*dble(q(:,label_l)) - dble(q(:,label_r))) * df2(2)
!      label_l = local_boundary_neighbors_NP(5,indx_NP)
!      label_r = local_boundary_neighbors_NP(6,indx_NP)
!      dQz = (-3.0d0*q_tmp + 4.d0*dble(q(:,label_l)) - dble(q(:,label_r))) * df2(3)
!      dQ = dQx*abs(normal_NP(1,indx_NP)) + dQy*abs(normal_NP(2,indx_NP)) + dQz*abs(normal_NP(3,indx_NP))
!
!      if (NP_nondegenerate) then
!         q_tmp_new =  q_tmp - dt_c * ( - dQ + gama(2)*(q_tmp - qboundary_NP(:,indx_NP)))
!      elseif (NP_degenerate) then
!         call surface(q_tmp,normal_NP(:,indx_NP),Qtrow,Qprow,nuQnu)
!         q_tmp_new =  q_tmp - dt_c * ( - dQ + 2.d0*gama(2)*((Qtrow-Qprow)-delta*nuQnu))
!      endif
   
   else
      q_tmp_new =  q(:,indx)
   endif
      q_new(:,indx) = q_tmp_new
enddo

call MPI_BARRIER(MPI_COMM_NEW,ierr)
call MPI_Win_fence(0, win, ierr)
Do ii=1, length
  q(:,ii) = q_new(:,ii) 
enddo
call MPI_Win_fence(0, win, ierr)
call MPI_BARRIER(MPI_COMM_NEW,ierr)

return
contains
subroutine dotq(q,dots)
use type
implicit none
real(dp) , dimension (6):: q
real(dp),dimension(6),intent(out)  :: dots
real (dp),dimension(3,3) :: qo
integer :: i,n

qo(1,1) = q(1)
qo(1,2) = q(2)
qo(1,3) = q(3)

qo(2,1) = q(2)
qo(2,2) = q(4)
qo(2,3) = q(5)

qo(3,1) = q(3)
qo(3,2) = q(5)
qo(3,3) = q(6)

dots = 0.d0
   do i=1,3
      dots(1) =  dots(1) + qo(1,i)* qo(i,1)
      dots(2) =  dots(2) + qo(1,i)* qo(i,2)
      dots(3) =  dots(3) + qo(1,i)* qo(i,3)
      dots(4) =  dots(4) + qo(2,i)* qo(i,2)
      dots(5) =  dots(5) + qo(2,i)* qo(i,3)
      dots(6) =  dots(6) + qo(3,i)* qo(i,3)
   enddo

end subroutine dotq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine surface(Qin,nu,Qt,Qp,nQn)
use var
use type
implicit none
real(dp), dimension(6) ,intent(in)::Qin
real(dp), dimension(3) ,intent(in)::nu
real(dp), dimension(6) ,intent(out)::Qt,Qp
real(dp), intent(out) :: nQn
real(dp), dimension(3,3) :: Qm, P, Qpp, Qtt
integer :: i,j,k,l
Qpp =0.d0;Qtt=0.d0;nQn=0.d0
Qm(1,1)=Qin(1);Qm(1,2)=Qin(2);Qm(1,3)=Qin(3);Qm(2,1)=Qin(2);Qm(2,2)=Qin(4)
Qm(2,3)=Qin(5);Qm(3,1)=Qin(3);Qm(3,2)=Qin(5);Qm(3,3)=Qin(6)

Do i=1,3
   Do j=1,3
      if (i==j) Qtt(i,j) = Qm(i,j) + Sbulk/3.d0
      if (i/=j) Qtt(i,j) = Qm(i,j) 
      
      if (i==j) P(i,j) = 1.d0 - nu(i)*nu(j)
      if (i/=j) P(i,j) = - nu(i)*nu(j)
   enddo
enddo

Do i=1,3
   Do j=1,3
      Do k=1,3
         Do l=1,3
            Qpp(i,j) = Qpp(i,j) + P(i,k) * Qtt(k,l) * P(l,j)
         enddo
      enddo
  enddo
enddo

Do i=1,3
   Do j=1,3
      nQn = nQn + nu(i) * Qtt(i,j) * nu(j)
   enddo
enddo

Qt(1)=Qtt(1,1);Qt(2)=Qtt(1,2);Qt(3)=Qtt(1,3)
Qt(4)=Qtt(2,2);Qt(5)=Qtt(2,3);Qt(6)=Qtt(3,3)


Qp(1)=Qpp(1,1);Qp(2)=Qpp(1,2);Qp(3)=Qpp(1,3)
Qp(4)=Qpp(2,2);Qp(5)=Qpp(2,3);Qp(6)=Qpp(3,3)

endsubroutine surface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
End subroutine ginzburg_landau

