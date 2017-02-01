subroutine scater
use var
use type
use parallel
implicit none
include 'mpif.h'
integer :: i, j, Wsize, indx
integer :: receive, send(size), displs(size+1)

Call MPI_BCAST  (ad,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_NEW,ierr)
Call MPI_BCAST  (vd,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_NEW,ierr)
Call MPI_BCAST  (df,3,MPI_DOUBLE_PRECISION,master,MPI_COMM_NEW,ierr)
Call MPI_BCAST  (df2,3,MPI_DOUBLE_PRECISION,master,MPI_COMM_NEW,ierr)
Call MPI_BCAST  (ddf,3,MPI_DOUBLE_PRECISION,master,MPI_COMM_NEW,ierr)
Call MPI_BCAST  (frame,1,MPI_INT,master,MPI_COMM_NEW,ierr)
Call MPI_BCAST  (length,1,MPI_INT,master,MPI_COMM_NEW,ierr)

!Wsize = dp*6*length
Wsize = i4b*6*length !use this for sp
CALL MPI_Win_allocate_shared(Wsize, 1, MPI_INFO_NULL, MPI_COMM_NEW, Pq, win,ierr)

call MPI_Win_fence(0, win, ierr)
Do i=1,length
   q(:,i) = 0
Enddo
call MPI_Win_fence(0, win, ierr)

Call MPI_scatter(q_reduced,6*length,MPI_REAL,q,6*length,MPI_REAL,master,MPI_COMM_NEW,ierr)

allocate(neighbors(6,length))
Call MPI_scatter(neighbors_reduced,6*length,MPI_INT,neighbors,6*length,MPI_INT,master,MPI_COMM_NEW,ierr)
Do i=1,length
   neighbors(:,i) = neighbors(:,i) - length*ID
Enddo

allocate(bulk(length))
bulk = .false.
Call MPI_scatter(bulk_reduced,length,MPI_BYTE,bulk,length,MPI_BYTE,master,MPI_COMM_NEW,ierr)

allocate(boundary(length))
boundary = .false.
Call MPI_scatter(boundary_reduced,length,MPI_BYTE,boundary,length,MPI_BYTE,master,MPI_COMM_NEW,ierr)


indx = 0
Do i=1,length
   if (boundary(i)) indx = indx + 1
enddo
receive = indx
CALL MPI_GATHER(receive, 1, MPI_INTEGER, send, 1, MPI_INTEGER, master, MPI_COMM_NEW, IERR) 
Allocate(qboundary(6,receive))
If (ID == master) then
   Displs = 0
   Do i = 1, size + 1
      Do j = 1, i-1
         Displs(i) = Displs(i) + send(j)
      enddo
   enddo
endif
CALL MPI_SCATTERV (qboundary_master, 6*send, 6*Displs, MPI_DOUBLE_PRECISION, qboundary, 6*receive, MPI_DOUBLE_PRECISION, master, MPI_COMM_NEW, IERR)

Allocate(normal(3,receive))
CALL MPI_SCATTERV (normal_master, 3*send, 3*Displs, MPI_DOUBLE_PRECISION, normal, 3*receive, MPI_DOUBLE_PRECISION, master, MPI_COMM_NEW, IERR)

If (ID==master) then
   Print *, "Q tensors are distributed to all processors"
   Print *, "Number of nodes per processor =", length
   deallocate(q_reduced,neighbors_reduced,qboundary_master)
   deallocate(bulk_reduced,boundary_reduced)
endif
End subroutine scater
