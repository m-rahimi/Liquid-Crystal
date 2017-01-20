subroutine convert
use var
use type
use parallel
implicit none
integer (i4b) :: i,j,k, N, label, indx

Mnodes = nodes + bnodes
length = int(Mnodes/size) + 1

allocate(lb2indx(tnodes))
allocate(q_reduced(6,length*size),neighbors_reduced(6,length*size))
allocate(bulk_reduced(length*size), boundary_reduced(length*size))

q_reduced = 0.0; lb2indx = 0
bulk_reduced = .false. ; boundary_reduced = .false.


label = 0
indx  = 0
do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
         label = label + 1
         if (boundary_master(i,j,k).or.bulk_master(i,j,k)) then
            indx = indx + 1
            q_reduced(1:5,indx) = q_master(:,label)
            q_reduced(6,indx) = -q_reduced(1,indx)-q_reduced(4,indx)
            lb2indx(label) = indx
            if (boundary_master(i,j,k)) boundary_reduced(indx) = .true.
            if (bulk_master(i,j,k)) bulk_reduced(indx) = .true.

         endif
      enddo
   enddo
enddo


label = 0
indx  = 0
do k=1,dime(3)
   do j=1,dime(2)
      do i=1,dime(1)
         label = label + 1
         if (boundary_master(i,j,k).or.bulk_master(i,j,k)) then
            indx = indx + 1
            Do N = 1, 6
               neighbors_reduced(N,indx) = lb2indx(neighbors_master(N,label))
            enddo
         endif
      enddo
   enddo
enddo

!deallocate(q_master,neighbors_master,boundary_master,bulk_master)
End subroutine convert
