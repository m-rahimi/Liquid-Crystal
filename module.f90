Module type
	integer, parameter :: dp  = kind(1.d00)
	integer, parameter :: sp  = kind(1.0)
	integer, parameter :: i4b = selected_int_kind (8)
    integer, parameter :: i2b = selected_int_kind (4)
    integer, parameter :: i1b = selected_int_kind (2)
end module type

Module Fix
	use type
	integer, parameter :: NI = 5 ! 5 independent parametrs in Q tensor
	real(dp), parameter :: pi = 3.141592653589793238462643383279502884
   real(dp), parameter :: kbT = 4.11437276d-21   ! joules/kelvin
   real(dp), parameter :: nematic_l = 7.15d0 ! nm 
   real(dp), parameter :: Elastic_L1 = 6.0d-12  ! J/m
   real(dp),dimension(6) ::    fac=(/ 1.d0, 2.d0, 2.d0, 1.d0, 2.d0, 1.d0 /)
End module Fix	
 
Module Var
	use type
	use fix
	
	 namelist /indata/ system,stype,anchor,box,dime,NPtype,NPanchor,radius,number_NP, lhat, gama, ULdG, chiral, dt,    &
                      accuracy, iseed ,sample, traj, Itype, restart, output, pitch 
    
!    logical :: restart
    Character(10) :: system, stype, nptype, anchor, NPanchor 
     
    integer (i4b), dimension(3) :: dime ! number of nodes in each direction
    integer (i4b), dimension(3) :: center ! center of droplet
!    real (dp), dimension(3) :: center ! center of droplet
    integer (i1b) :: number_NP ! number of NP
    integer (i4b) :: mc_total, sample, traj, iseed(3), sycle, frame
    integer (i4b) :: Tnodes, Bnodes, nodes, Dnodes, Mnodes, TNPbnodes ! number of nodes 
    integer (i4b) :: NPnodes, NPbnodes ! number of nodes
    integer (i1b) :: method ! method for simulation 1 MC; 2 GB ; 3 MC-GB
    integer (i1b) :: Itype ! initial configuration 1 random; 2 uniform
    integer (i1b) :: restart ! restart simulation 0 no; 1 yes
    integer (i1b) :: NPmove  ! Max radius of NP move
    integer (i1b) :: output ! type of output
    integer (i4b), dimension(:,:), allocatable :: neighbors_master
    integer (i4b), dimension(:,:), allocatable :: neighbors_reduced
    integer (i4b), dimension(:,:), allocatable :: neighbors
    integer (i4b), dimension(:), allocatable :: lb2indx
        
    real, dimension(:,:), allocatable :: NP_position

    real(dp) :: accuracy
    integer :: radius ! radius of NP
    real (dp) :: ULdG
    real (dp) :: qchiral, pitch ! chiral pich
    real (dp) :: beta ! kT for MC NP move
    real (dp) :: chiral ! winding number
    real (dp), dimension (4) :: lhat !elastic constant
    real (dp), dimension (2) :: gama !anchoring
    real (dp), dimension (3) :: dt 
    
    real(dp) :: Sbulk, Sinit
    real(dp) :: vd, vd_g, Volume
    real(dp) :: Escale, Escale_s
    real(dp) :: ad, ad_NP ! difrential of surface
    real(dp) :: xmax, ymax, zmax
    real(dp), dimension(3) :: df, df2, ddf, box
    
    real(sp), dimension(:,:), allocatable :: q_master, q_reduced, q_new ! Q tensor and the transformation of Q
    real(dp), dimension(:,:), allocatable :: qboundary_master, qboundary, qboundary_NP ! Q tensor for boundary
    real(dp), dimension(:), allocatable :: gEnergy ! ghost Energy of each node
    real(dp), dimension(:,:), allocatable :: trj_E
    real(dp), dimension(:,:), allocatable :: normal, normal_master, normal_NP, tangent_master, tangent_NP
    
    !SHARED MEMORY
    real(sp) :: q(6,*)
    Pointer (Pq, q)

    integer(i4b), dimension(:,:), allocatable :: local_neighbors
    
    logical(i1b), dimension(:,:,:), allocatable :: boundary_master !true for boundary nodes
    logical(i1b), dimension(:,:,:), allocatable :: Pbox ! true for active nodes for each Processors
    logical(i1b), dimension(:,:,:), allocatable :: NP !true for nodes inside NP
    logical(i1b), dimension(:,:,:), allocatable :: drop !true for nodes inside droplet
    logical(i1b), dimension(:,:,:), allocatable :: bulk_master !true for nodes inside droplet
    logical(i1b), dimension(:,:,:), allocatable :: global_NP, global_boundary_NP
    logical(i1b), dimension(:), allocatable :: boundary_reduced, boundary !true for boundary nodes
    logical(i1b), dimension(:), allocatable :: bulk_reduced, bulk !true for boundary nodes

    logical :: channel,droplet,PBC3
    logical :: infinite,nondegenerate,degenerate
    logical :: NP_infinite, NP_nondegenerate, NP_degenerate
	
end module var
module parallel
use type
integer :: ierr,size,id, master = 0
integer :: length
integer :: MPI_COMM_NEW, Win
end module parallel
