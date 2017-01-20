subroutine allocate_master
use var
use type
use fix
implicit none

allocate(q_master(5,tnodes),neighbors_master(6,tnodes))
allocate(boundary_master(dime(1),dime(2),dime(3)),bulk_master(dime(1),dime(2),dime(3)))

If (droplet) then
   allocate(drop(dime(1),dime(2),dime(3)))
endif
end subroutine allocate_master 
