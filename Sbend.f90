subroutine Sbend(label_c,spbendp)
use type
use var
use fix
implicit none
integer :: ii
integer :: x, y, left, right, label_l, label_r,label_c,label1,label2,label3,label4
real(dp), dimension(18) :: grad,grad1
real(dp) :: Elastic_Energy, Chiral_Energy,q_0,spbendp
real(dp), dimension(6) :: ddQ


!             Do x = 1 ,3 
!             ii = (x-1)*2 + 1


             label_l = neighbors_master(1,label_c)
             label_r = neighbors_master(2,label_c)
             ddQ(1)=dble(q(1, lb2indx(label_r))-2.d0*q(1,lb2indx(label_c))+q(1,lb2indx(label_l)))

             label_l = neighbors_master(3,label_c)
             label_r = neighbors_master(4,label_c)
             ddQ(4)=dble(q(4, lb2indx(label_r))-2.d0*q(4,lb2indx(label_c))+q(4,lb2indx(label_l)))

             label_l = neighbors_master(5,label_c)
             label_r = neighbors_master(6,label_c)
             ddQ(6)=dble(q(6, lb2indx(label_r))-2.d0*q(6,lb2indx(label_c))+q(6,lb2indx(label_l)))

!!!partial Dxy(Qxy)
             label_l = neighbors_master(1,label_c)
             label_r = neighbors_master(2,label_c)

             label1 =  neighbors_master(4,label_r)
             label2 =  neighbors_master(3,label_r)
             label3 =  neighbors_master(4,label_l)
             label4 =  neighbors_master(3,label_l)

             ddQ(2)=dble(q(2,lb2indx(label1))-q(2,lb2indx(label2))-q(2,lb2indx(label3))+q(2,lb2indx(label4)))*(0.5d0*df2(1))

!!!partial Dxz(Qxz)

             label_l = neighbors_master(1,label_c)
             label_r = neighbors_master(2,label_c)

             label1 =  neighbors_master(6,label_r)
             label2 =  neighbors_master(5,label_r)
             label3 =  neighbors_master(6,label_l)
             label4 =  neighbors_master(5,label_l)

             ddQ(3)=dble(q(3,lb2indx(label1))-q(3,lb2indx(label2))-q(3,lb2indx(label3))+q(3,lb2indx(label4)))*(0.5d0*df2(1))

!!!partial Dyz(Qyz)
             label_l = neighbors_master(3,label_c)
             label_r = neighbors_master(4,label_c)

             label1 =  neighbors_master(6,label_r)
             label2 =  neighbors_master(5,label_r)
             label3 =  neighbors_master(6,label_l)
             label4 =  neighbors_master(5,label_l)

             ddQ(5)=dble(q(5,lb2indx(label1))-q(5,lb2indx(label2))-q(5,lb2indx(label3))+q(5,lb2indx(label4)))*(0.5d0*df2(1))

!             grad1(1+(x-1)*6:x*6)
!             =dble(q(:,label_r)-q(:,label_l))*df2(x)
!             grad(1+(x-1)*6:x*6) =
!             grad1(1+(x-1)*6:x*6)*grad1(1+(x-1)*6:x*6)*fac
!             enddo
spbendp=ddQ(1)+ddQ(4)+ddQ(6)+2.d0*(ddQ(2)+ddQ(3)+ddQ(5))



end subroutine Sbend

