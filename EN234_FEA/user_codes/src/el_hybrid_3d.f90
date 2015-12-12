!     Subroutines for basic 3D linear elastic elements 



!==========================SUBROUTINE el_linelast_3dbasic ==============================
subroutine el_hybrid_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D

    use Element_Utilities, only : m =>  hybrid_shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D

    use Element_Utilities, only : dmdxi =>  hybrid_shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D

    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D

    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : dxdxi_red => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx=>vol_avg_shape_function_derivatives_3D

    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          

    ! Local Variables
    integer      :: n_points,kint,P_POINTS,n_points_red

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6),D_DEV(6,6),G(6,6)                           ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array),BM(6,8)            ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  B_red(6,length_dof_array)
   ! real (prec)  :: M(1),DMDXI(1,3)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
     real (prec)  ::  dxidx_red(3,3), determinant_red
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12                              ! Material properties
    real (prec)  :: K_UU(length_dof_array,length_dof_array),K_PP(8,8)
    real (prec)  :: K_UP(length_dof_array,8),k_PU(8,length_dof_array)
    real (prec)  ::K_PP_INV(8,8),K_PP_DET
    real (prec)  ::p(8),pint
     real (prec)  ::s(4,4),s_inv(4,4),sdet,s2(8,8),s2_inv(8,8),s2det
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    if (n_nodes == 4) n_points_red = 1
    if (n_nodes == 10) n_points_red = 4
    if (n_nodes == 8) n_points_red = 1
    if (n_nodes == 20) n_points_red = 8


    if (n_nodes == 4) p_points = 1
    if (n_nodes == 10) p_points = 4
    if (n_nodes == 8) p_points = 1
    if (n_nodes == 20) p_points = 8

    call initialize_integration_points(n_points, n_nodes, xi, w)



    element_residual = 0.d0
    element_stiffness = 0.d0
	p=0.d0
    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu)
    d11 = (1.D0-xnu)*E/( (1.d0+xnu)*(1.d0-2.D0*xnu) )
    d12 = xnu*E/( (1.d0+xnu)*(1.d0-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44

      D_DEV=0.D0
    D_DEV(1,1) = 2.d0
    D_DEV(2,2) = 2.d0
    D_DEV(3,3) = 2.d0
    D_DEV(4,4) = 1.d0
    D_DEV(5,5) = 1.d0
    D_DEV(6,6) = 1.d0
    D_DEV=D_DEV*E/2.d0/(1.D0+XNU)
    G=0.D0
    G(1:3,1:3)=1.D0
    D_DEV=D_DEV-G*E/3.D0/(1.D0+XNU)
  
    !     --  Loop over integration points
  K_UP=0.D0
  K_PU=0.D0
  K_PP=0.D0
  k_uu=0.d0

    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)

        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        K_UU(1:3*N_NODES,1:3*N_NODES)= K_uu(1:3*N_NODES,1:3*N_NODES)+&
         matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D_DEV,B(1:6,1:3*n_nodes)))*w(kint)*determinant


        bm=0.d0
        if (p_points>1) then

        call calculate_shapefunctions(xi(1:3,kint),p_points,m,dmdxi)


        BM(1,1:p_points)=m(1:p_points)
        BM(2,1:p_points)=m(1:p_points)
        BM(3,1:p_points)=m(1:p_points)

        else
        m(1)=1.d0
        BM(1:3,1:p_points)=1.d0

         end if
        ! call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
!write(*,*) m

        K_UP(1:3*N_NODES,1:P_POINTS)= K_UP(1:3*N_NODES,1:P_POINTS)+&
        MATMUL(TRANSPOSE(B(1:6,1:3*N_NODES)),BM(1:6,1:P_POINTS))*w(kint)*determinant
        K_PU(1:P_POINTS,1:3*N_NODES)= K_pu(1:P_POINTS,1:3*N_NODES)+&
        MATMUL(TRANSPOSE(BM(1:6,1:P_POINTS)),B(1:6,1:3*N_NODES))*w(kint)*determinant
        K_PP(1:P_POINTS,1:P_POINTS)= K_pp(1:P_POINTS,1:P_POINTS)+&
        spread(M(1:p_points),dim=2,ncopies=P_POINTS)*spread(M(1:p_points),dim=1,ncopies=P_POINTS)&
        *-3.D0*(1.D0-2.D0*XNU)/E*w(kint)*determinant




    end do
  !  call invert_small(K_PP(1:P_POINTS,1:P_POINTS),K_PP_INV(1:P_POINTS,1:P_POINTS),K_PP_DET)
if (p_points==1) then
k_pp_inv(1,1)=1.d0/k_pp(1,1)
else if (p_points==4)then
s(1:4,1:4)=k_pp(1:4,1:4)
call invert_small(s,s_inv,sdet)
k_pp_inv(1:4,1:4)=s_inv(1:4,1:4)
else if (p_points==8)then
s2(1:8,1:8)=k_pp(1:8,1:8)
call invert_small(s2,s2_inv,s2det)
k_pp_inv(1:8,1:8)=s2_inv(1:8,1:8)
end if
 p(1:p_points)=-1.d0*matmul(matmul(k_pp_inv(1:p_points,1:p_points),k_pu(1:p_points,1:3*n_nodes)),&
 dof_total(1:3*n_nodes)+dof_increment(1:3*n_nodes))




    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
         B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)



        if (p_points>1) then

        call calculate_shapefunctions(xi(1:3,kint),p_points,m,dmdxi)

        else
        m(1)=1.d0

         end if

        pint=dot_product(p(1:p_points),m(1:p_points))
        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)

        stress = matmul(D_dev,strain+dstrain)
        stress(1)=stress(1)+pint
        stress(2)=stress(2)+pint
        stress(3)=stress(3)+pint
        !write(*,*)stress(1),strain(1)+dstrain(1),pint,lmn,kint
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

        end do


 element_stiffness(1:3*n_nodes,1:3*n_nodes)=K_UU(1:3*N_NODES,1:3*N_NODES)-&
  MATMUL(MATMUL(K_UP(1:3*N_NODES,1:P_POINTS),K_PP_INV(1:P_POINTS,1:P_POINTS)),K_PU(1:P_POINTS,1:3*N_NODES))

    return
end subroutine el_hybrid_3dbasic


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_hybrid_3dbasic_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
    n_state_variables, initial_state_variables, &                                                        ! Input variables
    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
               
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          
    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    !
    !     Subroutine to compute element force vector for a linear elastodynamic problem
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
	
    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu) 
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44
  
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
      
        stress = matmul(D,strain+dstrain)
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

    end do
  
    return
end subroutine el_hybrid_3dbasic_dynamic


!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_hybrid_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : m =>  hybrid_shape_functions_3D
    use Element_Utilities, only : dmdxi =>  hybrid_shape_function_derivatives_3D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx=>vol_avg_shape_function_derivatives_3D
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
  
    ! Local Variables
    logical      :: strcmp
  


    real (prec)  ::  sdev(6)                           ! Deviatoric stress



    real (prec)  ::  smises                          ! Pressure and Mises stress


     integer      :: n_points,kint,P_POINTS,k

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6),D_DEV(6,6),G(6,6)                           ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array),BM(6,8)            ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  B_red(6,length_dof_array)
   ! real (prec)  :: M(1),DMDXI(1,3)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
     real (prec)  ::  dxidx_red(3,3), determinant_red
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12                              ! Material properties
    real (prec)  :: K_UU(length_dof_array,length_dof_array),K_PP(8,8)
    real (prec)  :: K_UP(length_dof_array,8),k_PU(8,length_dof_array)
    real (prec)  ::K_PP_INV(8,8),K_PP_DET
    real (prec)  ::p(8),pint
    real (prec)  ::s(4,4),s_inv(4,4),sdet,s2(8,8),s2_inv(8,8),s2det
    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio





    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27




    if (n_nodes == 4) p_points = 1
    if (n_nodes == 10) p_points = 4
    if (n_nodes == 8) p_points = 1
    if (n_nodes == 20) p_points = 8

    call initialize_integration_points(n_points, n_nodes, xi, w)



    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu)
    d11 = (1.D0-xnu)*E/( (1.d0+xnu)*(1.d0-2.D0*xnu) )
    d12 = xnu*E/( (1.d0+xnu)*(1.d0-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44
    D_DEV=0.D0
    D_DEV(1,1) = 2.d0
    D_DEV(2,2) = 2.d0
    D_DEV(3,3) = 2.d0
    D_DEV(4,4) = 1.d0
    D_DEV(5,5) = 1.d0
    D_DEV(6,6) = 1.d0
    D_DEV=D_DEV*E/2.d0/(1.D0+XNU)
    G=0.D0
    G(1:3,1:3)=1.D0
    D_DEV=D_DEV-G*E/3.D0/(1.D0+XNU)
  
    !     --  Loop over integration points
K_UP=0.D0
K_PU=0.D0
K_PP=0.D0
k_uu=0.d0

    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)

        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        K_UU(1:3*N_NODES,1:3*N_NODES)= K_uu(1:3*N_NODES,1:3*N_NODES)+&
         matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D_DEV,B(1:6,1:3*n_nodes)))*w(kint)*determinant


        bm=0.d0
        if (p_points>1) then

        call calculate_shapefunctions(xi(1:3,kint),p_points,m,dmdxi)


        BM(1,1:p_points)=m(1:p_points)
        BM(2,1:p_points)=m(1:p_points)
        BM(3,1:p_points)=m(1:p_points)

        else
        m(1)=1.d0
        BM(1:3,1:p_points)=1.d0

         end if


        K_UP(1:3*N_NODES,1:P_POINTS)= K_UP(1:3*N_NODES,1:P_POINTS)+&
        MATMUL(TRANSPOSE(B(1:6,1:3*N_NODES)),BM(1:6,1:P_POINTS))*w(kint)*determinant
        K_PU(1:P_POINTS,1:3*N_NODES)= K_pu(1:P_POINTS,1:3*N_NODES)+&
        MATMUL(TRANSPOSE(BM(1:6,1:P_POINTS)),B(1:6,1:3*N_NODES))*w(kint)*determinant
        K_PP(1:P_POINTS,1:P_POINTS)= K_pp(1:P_POINTS,1:P_POINTS)+&
        spread(M(1:p_points),dim=2,ncopies=P_POINTS)*spread(M(1:p_points),dim=1,ncopies=P_POINTS)&
        *-3.D0*(1.D0-2.D0*XNU)/E*w(kint)*determinant

!write(*,*)n


    end do
  if (p_points==1) then
k_pp_inv(1,1)=1.d0/k_pp(1,1)
else if (p_points==4)then
s(1:4,1:4)=k_pp(1:4,1:4)
call invert_small(s,s_inv,sdet)
k_pp_inv(1:4,1:4)=s_inv(1:4,1:4)
else if (p_points==8)then
s2(1:8,1:8)=k_pp(1:8,1:8)
call invert_small(s2,s2_inv,s2det)
k_pp_inv(1:8,1:8)=s2_inv(1:8,1:8)
end if

 p(1:p_points)=-1.d0*matmul(matmul(k_pp_inv(1:p_points,1:p_points),k_pu(1:p_points,1:3*n_nodes)),&
 dof_total(1:3*n_nodes)+dof_increment(1:3*n_nodes))



 nodal_fieldvariables = 0.d0

        do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
         B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)



        if (p_points>1) then

        call calculate_shapefunctions(xi(1:3,kint),p_points,m,dmdxi)

        else
        m(1)=1.d0

         end if
         call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
!write(*,*) n
        pint=dot_product(p(1:p_points),m(1:p_points))
        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)


        stress = matmul(D_dev,strain+dstrain)
        sdev=stress

        stress(1)=stress(1)+pint
        stress(2)=stress(2)+pint
        stress(3)=stress(3)+pint



        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            end if
        end do
    end do


  
    return
end subroutine fieldvars_hybrid_3dbasic

