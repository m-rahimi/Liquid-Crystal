!----------------------------------------------------------
!----------------------------------------------------------
subroutine qmga(cpu,step)
  use type
  use var
  implicit none
  integer (i4b) ,intent(in) :: step,cpu
  integer (i4b) :: i,j,k,label, coun = 0, N
  real(dp),dimension(3) :: e,e_ref,rotvec,director
  real(dp) :: qw,qx,qy,qz,pdot,norm_a,norm_b,angle,jj,kk,mm,norm_vec
  real(dp), dimension(3) :: scalar
  real(dp), dimension(5) ::   q_tmp
   integer (i4b) ::  iout=50, zero = 0
!*********************************************************************
  Character(len=10) :: namefile1, namefile2
  Character(len=20), Parameter :: sfile_1= 'traj'
  
  write(namefile1,'(I10)')int(cpu)
  write(namefile2,'(I10)')int(step/traj)
  open(iout,FILE=trim(adjustr(sfile_1))//trim(adjustl(namefile1))//"."//trim(adjustl(namefile2)),STATUS='replace')

!**************************************************************************************       
  
  write(iout,100)'ITEM: TIMESTEP'
  write(iout,102)step
  write(iout,100)'ITEM: NUMBER OF ATOMS'
  write(iout,102) nodes + bnodes + number_NP + NPbnodes 
  write(iout,100)'ITEM: BOX BOUND ff ff ff'
  write(iout,101) real(zero), real(dime(1))
  write(iout,101) real(Zero), real(dime(2))
  write(iout,101) real(zero), real(dime(3))
  write(iout,100)'ITEM: ATOMS id type x y z c_qdat[1] c_qdat[2] c_qdat[3] c_qdat[4] '

  coun = 0
  e_ref(1)=0.d0
  e_ref(2)=0.d0
  e_ref(3)=1.d0


  do k=1,dime(3)
     do j=1,dime(2)
        do i=1,dime(1)
           coun = coun + 1
           label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
           q_tmp = q_master(1:5,label)
           call eigen( q_tmp, scalar,director(:))

           e=director(:) 
           norm_a = sqrt(e(1)**2 + e(2)**2 + e(3)**2)
! Normalize the vectors

           e(1)=e(1)/norm_a
           e(2)=e(2)/norm_a
           e(3)=e(3)/norm_a

           if (e_ref(1) == e(1).and.e_ref(2)== e(2).and.(e(3)==e_ref(3).or.e(3)==-e_ref(3)))then
              rotvec(1) = 1.d0
              rotvec(2) = 0.d0
              rotvec(3) = 0.d0
           else 
              rotvec(1) =   e_ref(2)*e(3)- e_ref(3)*e(2)
              rotvec(2) = -(e_ref(1)*e(3)- e_ref(3)*e(1))
              rotvec(3) =   e_ref(1)*e(2)- e_ref(2)*e(1)
           endif
      
           norm_vec = sqrt(rotvec(1)**2 + rotvec(2)**2 + rotvec(3)**2)
      
           rotvec(1) = rotvec(1)/norm_vec
           rotvec(2) = rotvec(2)/norm_vec
           rotvec(3) = rotvec(3)/norm_vec
      
           angle = acos(e(1)*e_ref(1) + e(2)*e_ref(2) + e(3)*e_ref(3))
      
           qw = cos(angle/2.d0)
           qx = rotvec(1) * sin(angle/2.d0)
           qy = rotvec(2) * sin(angle/2.d0)
           qz = rotvec(3) * sin(angle/2.d0)
     
           if(drop(i,j,k))then
!              write(iout,458)coun,1,real(i),real(j),real(k),qw,qx,qy,qz,' '
           elseif (boundary_master(i,j,k)) then
              write(iout,458)coun,2,real(i),real(j),real(k),qw,qx,qy,qz,' '
           endif
      end do
    end do
  end do
  if (number_NP > 0) then
  Do i=1, number_NP
    coun = coun + 1
    write(iout,458)coun,3,real(NP_position(1,i)),real(NP_position(2,i)),real(NP_position(3,i)),qw,qx,qy,qz,' '
  enddo
  endif
  close(iout)
458 format(I10,1x,I3,1x,7(f10.3),A)        
  
100 format(A)
101 format(2f7.2)
102 format (I8)
end subroutine qmga

!----------------------------------------------------------
!----------------------------------------------------------

subroutine paraview_eigen_new(cpu,step)
  use type
  use var
  implicit none
  integer (i4b) ,intent(in) :: cpu,step
  
  integer (i4b) :: i,iout=50,j,k,label
  real(dp), dimension(5) ::   q_tmp
  real(dp) :: r2,lbox,cx,scalar(3),x,y,z,eta,splaybend
  real(dp),dimension(3) :: dx,director,xs
!*********************************************************************
  Character(len=10) :: namefile1, namefile2
  Character(len=20), Parameter :: sfile_3= 'surface',end='.vtk'  
  
  write(namefile1,'(I10)')int(cpu)
  write(namefile2,'(I10)')int(step)
  open(iout,FILE=trim(adjustr(sfile_3))//trim(adjustl(namefile1))//"."//trim(adjustl(namefile2))//trim(adjustl(end)),STATUS='replace')
  
  xs(1) = 0.5d0*xmax
  xs(2) = 0.5d0*ymax
  xs(3) = 0.5d0*zmax

!**************************************************************************************

  write(iout,100)'# vtk DataFile Version 3.0'
  write(iout,100)'vtk output'
  write(iout,100)'ASCII'
  write(iout,100)'DATASET STRUCTURED_GRID'
  write(iout,110)'DIMENSIONS',dime(1),dime(2),dime(3)
  write(iout,120)'POINTS',dime(1)*dime(2)*dime(3) ,'  float'

  do k=1,dime(3)
     do j=1,dime(2)
        do i=1,dime(1)
           write(iout,151)  real(i),real(j),real(k),'  '
        enddo
     enddo
  enddo


  write(iout,*)
  write(iout,120)'POINT_DATA', dime(1)*dime(2)*dime(3) 
  write(iout,100)'SCALARS scalars float'
  write(iout,100)'LOOKUP_TABLE default'

  
  do k=1,dime(3)
     do j=1,dime(2)
        do i=1,dime(1)
              label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
              q_tmp = q(1:5,label)
              call eigen ( q_tmp(:), scalar, director(:))
           if (global_NP(i,j,k)) then 
              write(iout,152 )  1.000,'  '
           else
              write(iout,152 )  maxval(scalar)*1.5d0,'  '
           endif
        enddo
     enddo
  enddo

  write(iout,*)
  write(iout,100)'SCALARS biaxial float'
  write(iout,100)'LOOKUP_TABLE default'

  
  do k=1,dime(3)
     do j=1,dime(2)
        do i=1,dime(1)
              label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
              q_tmp = q(1:5,label)
              call eigen(q_tmp(:), scalar,director(:))
              write(iout,152 ) scalar(2)*3.d0 + maxval(scalar)*1.5d0,'  ' 
        enddo
     enddo
  enddo


  write(iout,*)
  write(iout,100)'VECTORS directors float'

  do k=1,dime(3)
     do j=1,dime(2)
        do i=1,dime(1)
              label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
              q_tmp = q(1:5,label)
              call eigen ( q_tmp(:), scalar, director(:))
              if (director(3)>0) then 
                write(iout,153 )  director,'  '
              else
                 write(iout,153 )  -director,'  '
              endif
        enddo
     enddo
  enddo

 write(iout,*)
  write(iout,100)'SCALARS SplayBend float'
  write(iout,100)'LOOKUP_TABLE default'

  do k=1,dime(3)
     do j=1,dime(2)
        do i=1,dime(1)
              label = 1 + (i-1) + (j-1)*dime(1) + (k-1)*dime(1)*dime(2)
              q_tmp = q(1:5,label)
  !            if(bulk(i,j,k)) then
              call Sbend (label,splaybend)
              write(iout,152 ) splaybend,'  '
   !        else
    !          write(iout,152 )  -1.0,'  '
   !        endif
        enddo
     enddo
  enddo

  
  close (iout)
!*****************************************************************************************
!*****************************************************************************************
100 format('',A)
120 format('',A,I12,A)
110 format('',A, I6, 1X, I6, 1X, I6)
151 format (3f12.6 ,A)
152 format (f12.6 ,A)
153 format (3f12.6 ,A)
end subroutine paraview_eigen_new
