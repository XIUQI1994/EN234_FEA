!     Subroutines for basic 3D linear elastic elements 



!==========================SUBROUTINE el_linelast_3dbasic ==============================
subroutine el_hypoelastic_3d(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
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
    integer      :: n_points,kint

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6),D1(6,6),D2(6,6),D3(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dNdx_average(6,length_dof_array),Bmodifer(6,length_dof_array)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: s0,e0,n0,k             ! Material properties
    real (prec)  :: vol
    real (prec)  :: evol,emises,se,Es,Et
    real (prec)  :: edev(6)


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

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
	
    D = 0.d0
    D2=0.D0
    D2(1,1)=2
    D2(2,2)=2
    D2(3,3)=2
    D2(4,4)=1
    D2(5,5)=1
    D2(6,6)=1
    D3=0.D0
    D3(1:3,1:3)=1
    s0= element_properties(1)
    e0 = element_properties(2)
    n0= element_properties(3)
    k= element_properties(4)
!    d44 = 0.5D0*E/(1.d0+xnu)
!    d11 = (1.D0-xnu)*E/( (1.d0+xnu)*(1.d0-2.D0*xnu) )
!    d12 = xnu*E/( (1.d0+xnu)*(1.d0-2.D0*xnu) )
!    D(1:3,1:3) = d12
!    D(1,1) = d11
!    D(2,2) = d11
!    D(3,3) = d11
!    D(4,4) = d44
!    D(5,5) = d44
!    D(6,6) = d44
  
    !     --  Loop over integration points
    dNbardx=0.d0
    vol=0.d0
    do kint=1,n_points
    call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        dNbardx=dNbardx+dNdx*w(kint)*determinant
        vol=vol+w(kint)*determinant
    end do
        dNbardx=dNbardx/vol/3.d0
        dNdx_average=0.d0
        dNdx_average(1,1:3*n_nodes-2:3)=dNbardx(1:n_nodes,1)
        dNdx_average(1,2:3*n_nodes-1:3)=dNbardx(1:n_nodes,2)
        dNdx_average(1,3:3*n_nodes:3)=dNbardx(1:n_nodes,3)
        dNdx_average(2,1:3*n_nodes-2:3)=dNbardx(1:n_nodes,1)
        dNdx_average(2,2:3*n_nodes-1:3)=dNbardx(1:n_nodes,2)
        dNdx_average(2,3:3*n_nodes:3)=dNbardx(1:n_nodes,3)
        dNdx_average(3,1:3*n_nodes-2:3)=dNbardx(1:n_nodes,1)
        dNdx_average(3,2:3*n_nodes-1:3)=dNbardx(1:n_nodes,2)
        dNdx_average(3,3:3*n_nodes:3)=dNbardx(1:n_nodes,3)

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

        if ( element_identifier ==1002) then
        Bmodifer=0.d0
        Bmodifer(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)/3.d0
        Bmodifer(1,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)/3.d0
        Bmodifer(1,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)/3.d0
        Bmodifer(2,1:3*n_nodes)=Bmodifer(1,1:3*n_nodes)
        Bmodifer(3,1:3*n_nodes)=Bmodifer(1,1:3*n_nodes)

        B=B-Bmodifer+dNdx_average
        end if



        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
      !


        evol=0.d0
        edev=0.d0
        emises=0.d0
        se=0.d0

        strain=strain+dstrain
        evol = sum(strain(1:3))/3.d0
        edev = strain
        edev(1:3) = edev(1:3)-evol
        edev(4:6)=edev(4:6)/2.d0
        emises = ( (dot_product(edev(1:3),edev(1:3)) + 2.d0*dot_product(edev(4:6),edev(4:6)) )/1.5d0)**0.5d0
        !emises = dsqrt( (dot_product(edev(1:3),edev(1:3)) + 2.d0*dot_product(edev(4:6),edev(4:6)) )/1.5d0)
        if (emises>e0)then
        se=((emises/e0)**(1.d0/n0))*s0
        Et=s0/e0/n0
        Et=Et*((emises/e0)**(1.d0/n0-1))
        else
        se=(((1.d0+n0*n0)/(n0-1.d0)/(n0-1.d0)-(n0/(n0-1.d0)-emises/e0)**2.d0)**0.5d0-1.d0/(n0-1.d0))*s0
       !Et=(((1.d0+n0*n0)/(n0-1.d0)/(n0-1.d0)-(n0/(n0-1.d0)-emises/e0)**2.d0)**-0.5d0)*(n0/(n0-1.d0)-emises/e0)*s0/e0
        Et = s0*(n0/(n0-1.d0)-emises/e0)/e0
        Et = Et/dsqrt( (1.d0+n0*n0)/(n0-1.d0)/(n0-1.d0) - (n0/((n0-1.d0)-emises/e0)**2.d0 ))
        end if

        if(emises>0.d0) then
        Es=se/emises


!
        stress=se/emises*edev*2.d0/3.d0
        stress(1:3)=stress(1:3)+evol*3.d0*k
        else
           stress = 0.d0
           Es = Et
        endif

  !caculate D
 write(IOW,*) stress
        if (emises==0.d0)then
        D=Et/3.d0*D2+(k-2.D0*Et/9.D0)*D3
        else

        D1 = spread(edev,dim=2,ncopies=6)*spread(edev,dim=1,ncopies=6)
        D=4.D0/9.D0/emises/emises*(Et-Es)*D1+Es/3.d0*D2+(k-2.D0*Es/9.D0)*D3
        end if


   !     call hypoelastic_material(strain+dstrain,properties,n_properties,stress,D)

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D,B(1:6,1:3*n_nodes)))*w(kint)*determinant

    end do
 ! write(IOW,*) element_residual
    return
end subroutine el_hypoelastic_3d


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_hypoelastic_3d_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
end subroutine el_hypoelastic_3d_dynamic


!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_hypoelastic_3d(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
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
  
    integer      :: n_points,kint,k0

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6),D1(6,6),D2(6,6),D3(6,6)    ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dNdx_average(6,length_dof_array),Bmodifer(6,length_dof_array)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: s0,e0,n0,k             ! Material properties
    real (prec)  :: vol
    real (prec)  :: evol,emises,se,Es,Et,p,smises
    real (prec)  :: edev(6)
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0

     D = 0.d0
    D2=0.D0
    D2(1,1)=2
    D2(2,2)=2
    D2(3,3)=2
    D2(4,4)=1
    D2(5,5)=1
    D2(6,6)=1
    D3=0.D0
    D3(1:3,1:3)=1
    s0= element_properties(1)
    e0 = element_properties(2)
    n0= element_properties(3)
    k= element_properties(4)
	
!    D = 0.d0
!    E = element_properties(1)
!    xnu = element_properties(2)
!    d44 = 0.5D0*E/(1+xnu)
!    d11 = (1.D0-xnu)*E/( (1.d0+xnu)*(1.d0-2.D0*xnu) )
!    d12 = xnu*E/( (1.d0+xnu)*(1.d0-2.D0*xnu) )
!    D(1:3,1:3) = d12
!    D(1,1) = d11
!    D(2,2) = d11
!    D(3,3) = d11
!    D(4,4) = d44
!    D(5,5) = d44
!    D(6,6) = d44
  
    !     --  Loop over integration points
    dNbardx=0.d0
    vol=0.d0
    do kint=1,n_points
    call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        dNbardx=dNbardx+dNdx*w(kint)*determinant
        vol=vol+w(kint)*determinant
    end do
        dNbardx=dNbardx/vol/3.d0
        dNdx_average=0.d0
        dNdx_average(1,1:3*n_nodes-2:3)=dNbardx(1:n_nodes,1)
        dNdx_average(1,2:3*n_nodes-1:3)=dNbardx(1:n_nodes,2)
        dNdx_average(1,3:3*n_nodes:3)=dNbardx(1:n_nodes,3)
        dNdx_average(2,1:3*n_nodes-2:3)=dNbardx(1:n_nodes,1)
        dNdx_average(2,2:3*n_nodes-1:3)=dNbardx(1:n_nodes,2)
        dNdx_average(2,3:3*n_nodes:3)=dNbardx(1:n_nodes,3)
        dNdx_average(3,1:3*n_nodes-2:3)=dNbardx(1:n_nodes,1)
        dNdx_average(3,2:3*n_nodes-1:3)=dNbardx(1:n_nodes,2)
        dNdx_average(3,3:3*n_nodes:3)=dNbardx(1:n_nodes,3)

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

        if ( element_identifier ==1002) then
        Bmodifer=0.d0
        Bmodifer(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)/3.d0
        Bmodifer(1,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)/3.d0
        Bmodifer(1,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)/3.d0
        Bmodifer(2,1:3*n_nodes)=Bmodifer(1,1:3*n_nodes)
        Bmodifer(3,1:3*n_nodes)=Bmodifer(1,1:3*n_nodes)

        B=B-Bmodifer+dNdx_average
        end if

          strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
      !


        evol=0.d0
        edev=0.d0
        emises=0.d0
        se=0.d0

        strain=strain+dstrain
        evol = sum(strain(1:3))/3.d0
        edev = strain
        edev(1:3) = edev(1:3)-evol
        edev(4:6)=edev(4:6)/2.d0
        emises = dsqrt( dot_product(edev(1:3),edev(1:3)) + 2.d0*dot_product(edev(4:6),edev(4:6)) )*dsqrt(1.5d0)
        if (emises>e0)then
        se=((emises/e0)**(1.d0/n0))*s0
        Et=s0/e0/n0
        Et=Et*((emises/e0)**(1.d0/n0-1.d0))
        else
        se=(dsqrt((1.d0+n0*n0)/(n0-1.d0)/(n0-1.d0)-(n0/(n0-1.d0)-emises/e0)**2.d0)-1.d0/(n0-1.d0))*s0
        Et=(((1.d0+n0*n0)/(n0-1.d0)/(n0-1.d0)-(n0/(n0-1.d0)-emises/e0)**2.d0)**-0.5d0)*(n0/(n0-1.d0)-emises/e0)*s0/e0
        end if
        if(emises>0.d0) then
        Es=se/emises


!
        stress=se/emises*edev*2.d0/3.d0
        stress(1:3)=stress(1:3)+evol*3.d0*k
        else
           stress = 0.d0
           Es = Et
        endif


        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k0 = 1,n_field_variables
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
            endif
        end do
 
    end do
  
    return
end subroutine fieldvars_hypoelastic_3d

!subroutine hypoelastic_material(strain,properties,n_properties,stress,D)
!
!   use Types
!   use ParamIO
!
!   implicit none
!
!
!
!   real (prec), intent ( in )  :: strain(6)
!   real (prec), intent (in)  ::  properties(n_properties)
!
!
!   end subroutine hypoelastic_material
