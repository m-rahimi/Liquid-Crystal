subroutine free_energy(energy)
use type
use var
use parallel
implicit none
include 'mpif.h'
integer :: indx, indx_S, indx_NP, ii
real(dp) :: LdG = 0.d0, Elastic(2) = 0.d0, surface = 0.d0, surface_NP = 0.d0
real(dp), intent(out) :: energy(5) 
real(dp), dimension(6) :: q_tem

ii = 0
indx_S = 0
indx_NP = 0
LdG = 0.d0
Elastic = 0.d0
surface = 0.d0
surface_NP = 0.d0

do indx=1,length
   
   q_tem = dble(q(:,indx))
   if (bulk(indx) == .true.) then
      LdG = LdG + LdG_Energy(q_tem)
      Elastic = Elastic + Elastic_Energy(indx,q_tem)
   elseif (boundary(indx) == .true..and.infinite==.false.) then
      indx_S = indx_S + 1
      surface = surface + surface_Energy(indx,indx_S)
   endif

enddo

CALL MPI_BARRIER(MPI_COMM_NEW, ierr)
CALL MPI_ALLREDUCE (LdG, energy(1), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, IERR)
CALL MPI_ALLREDUCE (Elastic(1), energy(2), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, IERR)
CALL MPI_ALLREDUCE (Elastic(2), energy(3), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, IERR)
CALL MPI_ALLREDUCE (surface, energy(4), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, IERR)
CALL MPI_ALLREDUCE (surface_NP, energy(5), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, IERR)

energy(1) =  energy(1) * vd
energy(2) =  energy(2) * vd
energy(3) =  energy(3) * vd
if (nondegenerate) energy(4) =  0.5d0*gama(1)*energy(4)*ad 
if (degenerate) energy(4) =  gama(1)*energy(4)*ad

if (NP_nondegenerate) energy(5) =  0.5d0*gama(2)*energy(5)*ad_NP 
if (NP_degenerate) energy(5) =  gama(2)*energy(5)*ad_NP

energy(6) =  sum(energy(1:5))


return
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Function LdG_Energy(q_tem)
use type
use var
implicit none
real(dp), dimension(6) :: q_tem
real(dp) :: tr1, tr2, tr3
real(dp) :: aux1, a1, a4
real(dp) ::  LdG_Energy

tr1 = sum(fac*q_tem*q_tem)

tr2 = tr1 * tr1

aux1 = 3.d0/2.d0
a1 = -q_tem(1) - q_tem(4)
a4 = q_tem(1) - q_tem(4)

tr3 = 0.75d0*a1**3 - 3.d0*a1*q_tem(2)**2 + aux1*a1*q_tem(3)**2 + aux1*q_tem(3)**2*a4 &
       - 0.5d0*aux1*a1*a4**2 + 6.d0*q_tem(2)*q_tem(3)*q_tem(5) + aux1*a1*q_tem(5)**2 - aux1*a4*q_tem(5)**2
LdG_Energy = 0.5d0*( 1.d0 - ULdG/3.d0 )*tr1 - ULdG/3.d0*tr3 + ULdG/4.d0*tr2
return
end function LdG_Energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Function  Elastic_Energy(label_c,q_t)
use type
use var
use fix
use parallel
implicit none
integer :: ii
integer :: x, y, left, right, label_l, label_r, label_c
real(dp), dimension(18) :: grad, grad1
real(dp), dimension(6) :: q_t
real(dp) :: Elastic_Energy(2), Chiral_Energy

Do x = 1 , 3
   ii = (x-1)*2 + 1
   label_l = neighbors(ii,label_c)
   label_r = neighbors(ii+1,label_c)

   grad1(1+(x-1)*6:x*6) = (q(:,label_r)-q(:,label_l))*df2(x)
   grad(1+(x-1)*6:x*6) = grad1(1+(x-1)*6:x*6)*grad1(1+(x-1)*6:x*6)*fac
enddo

Chiral_Energy = q_t(1) * (-grad1(14)+grad1(9) )+ &
                q_t(2) * ( grad1(13)-grad1(16)-grad1(3)+grad1(11) ) +&
                q_t(3) * (-grad1(7)+grad1(2)-grad1(17)+grad1(12) ) +&
                q_t(4) * (-grad1(5)+grad1(14) ) +&
                q_t(5) * (-grad1(8)+grad1(15)+grad1(4)-grad1(6) )  +&
                q_t(6) * (-grad1(9)+grad1(5) )



Elastic_Energy(1) = 0.5d0*sum(grad)
Elastic_Energy(2) = 2.0d0*qchiral*Chiral_Energy

return
end function Elastic_Energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Function  surface_Energy(label_c,indx)
use type
use var
implicit none
integer,intent(in) ::indx, label_c
real(dp) :: rapini, surface_Energy, galota
real(dp), dimension(6) :: Qin
real(dp), dimension(6) :: q_temp

if (nondegenerate) then
   rapini = 0.d0
   rapini = sum (fac*(q(:,label_c)-qboundary(:,indx))*(q(:,label_c)-qboundary(:,indx)))
   surface_Energy = rapini
elseif (degenerate) then
   galota = 0
   q_temp = q(:,label_c)
   call surface_deg(q_temp,normal(:,indx),Qin)
   galota = sum(fac*Qin*Qin)
  surface_Energy = galota
endif
end function surface_Energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Function  surface_Energy_NP(label_c,indx)
use type
use var
implicit none
integer,intent(in) ::indx, label_c
real(dp) :: rapini, surface_Energy_NP, galota
real(dp), dimension(6) :: Qin
real(dp), dimension(6) :: q_temp

if (NP_nondegenerate) then
   rapini = 0.d0
   rapini = sum (fac*(q(:,label_c)-qboundary_NP(:,indx))*(q(:,label_c)-qboundary_NP(:,indx)))
   surface_Energy_NP = rapini
elseif (NP_degenerate) then
   galota = 0
   q_temp = q(:,label_c)
   call surface_deg(q_temp,normal_NP(:,indx),Qin)
   galota = sum(fac*Qin*Qin)
  surface_Energy_NP = galota
endif
end function surface_Energy_NP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine surface_deg(Qin,nu,Qout)
use var
use type
implicit none
real(dp), dimension(6) ,intent(in)::Qin
real(dp), dimension(3) ,intent(in)::nu
real(dp), dimension(6) ,intent(out)::Qout
real(dp), dimension(3,3) :: Qm, P, Qpp, Qtt
integer :: i,j,k,l
Qpp =0.d0;Qtt=0.d0
Qm(1,1)=Qin(1);Qm(1,2)=Qin(2);Qm(1,3)=Qin(3);Qm(2,1)=Qin(2);Qm(2,2)=Qin(4)
Qm(2,3)=Qin(5);Qm(3,1)=Qin(3);Qm(3,2)=Qin(5);Qm(3,3)=Qin(6)

Do i=1,3
   Do j=1,3
      if (i==j) Qtt(i,j) = Qm(i,j) + Sbulk/3.d0
      if (i/=j) Qtt(i,j) = Qm(i,j) 
      
      if (i==j) P(i,j) = 1 - nu(i)*nu(j)
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

Qout(1)=Qtt(1,1)-Qpp(1,1)
Qout(2)=Qtt(1,2)-Qpp(1,2)
Qout(3)=Qtt(1,3)-Qpp(1,3)
Qout(4)=Qtt(2,2)-Qpp(2,2)
Qout(5)=Qtt(2,3)-Qpp(2,3)
Qout(6)=Qtt(3,3)-Qpp(3,3)
endsubroutine surface_deg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


End subroutine  free_energy


