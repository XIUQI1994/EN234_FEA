!     Subroutines for basic 3D explicit dynamic elements



!==========================SUBROUTINE el_linelast_3dbasic ==============================
subroutine el_fracture_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dNdx_average(6,length_dof_array),Bmodifer(6,length_dof_array)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    real (prec)  :: vol
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
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1.d0+xnu)
    d11 = (1.D0-xnu)*E/( (1.d0+xnu)*(1.d0-2.D0*xnu) )
    d12 = xnu*E/( (1.d0+xnu)*(1.d0-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44
  
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

        if ( element_identifier ==1001) then
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
      
        stress = matmul(D,strain+dstrain)
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D,B(1:6,1:3*n_nodes)))*w(kint)*determinant

    end do
  
    return
end subroutine el_fracture_3dbasic


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_fracture_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
    integer      ::  Ndelete
    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6),stress0(6),stress1(6),dstress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]

    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node

    real (prec)  :: I(3,3)
    real (prec)  :: dof_2d(3,length_coord_array/3),dof_inc2d(3,length_coord_array/3)
    real (prec)  :: delta_F(3,3),mid_F(3,3),mid_F_inver(3,3),J,delta_L(3,3),delta_L_vol
    real (prec)  :: dNdy(length_dof_array/3,3),dNbardy(length_dof_array/3,3)
    real (prec)  :: vol,eta,delta_eta
    real (prec)  ::  dNdy_average(6,length_dof_array),Bmodifer(6,length_dof_array)
    real (prec)  :: dstrain_2d(3,3),dspin(3,3),drotation(3,3),r(3,3),r_inv(3,3),rdeter
    real (prec)  ::updated_state_variables_currentpoint(8),initial_state_variables_currentpoint(8)
    real (prec)  :: E, xnu, v,D44, D11, D12,D(6,6),vf,ematrix

 real (prec)  :: stress0_2d(3,3),dev_stress0(3,3)
  real (prec)  :: dev_dstrain(3,3)
 real (prec)  :: predictor_devstress(3,3),predictor_devstress_1d(6),predictor_volstress


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

  I=0.d0
  I(1,1)=1.d0
  I(2,2)=1.d0
  I(3,3)=1.d0

  vol=0.d0
  eta=0.d0
  dNbardy=0.d0
  delta_eta=0.d0


    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    v=xnu
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
        dof_2d=reshape(dof_total,(/3,length_coord_array/3/))
        dof_inc2d = reshape(dof_increment,(/3,length_coord_array/3/))
        delta_F(1:3,1:3)=matmul(dof_inc2d(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))

        mid_F=matmul(0.5d0*dof_inc2d(1:3,1:n_nodes)+dof_2d(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))

        mid_F(1,1)=mid_F(1,1)+1.d0
        mid_F(2,2)=mid_F(2,2)+1.d0
        mid_F(3,3)=mid_F(3,3)+1.d0

        call invert_small(mid_F,mid_F_inver,J)
        dNdy(1:n_nodes,1:3)=matmul(dNdx(1:n_nodes,1:3),mid_F_inver(1:3,1:3))

        delta_L=matmul(delta_F,mid_F_inver)

        vol=vol+w(kint)*determinant
        eta=eta+J*w(kint)*determinant
        delta_eta=delta_eta+J*(delta_L(1,1)+delta_L(2,2)+delta_L(3,3))*w(kint)*determinant
        dNbardy=dNbardy+dNdy*J*w(kint)*determinant
       end do

       eta=eta/vol
       dNbardy=dNbardy/eta/vol/3.d0
       delta_eta=delta_eta/eta/vol
        dNdy_average=0.d0
        dNdy_average(1,1:3*n_nodes-2:3)=dNbardy(1:n_nodes,1)
        dNdy_average(1,2:3*n_nodes-1:3)=dNbardy(1:n_nodes,2)
        dNdy_average(1,3:3*n_nodes:3)=dNbardy(1:n_nodes,3)
        dNdy_average(2,1:3*n_nodes-2:3)=dNbardy(1:n_nodes,1)
        dNdy_average(2,2:3*n_nodes-1:3)=dNbardy(1:n_nodes,2)
        dNdy_average(2,3:3*n_nodes:3)=dNbardy(1:n_nodes,3)
        dNdy_average(3,1:3*n_nodes-2:3)=dNbardy(1:n_nodes,1)
        dNdy_average(3,2:3*n_nodes-1:3)=dNbardy(1:n_nodes,2)
        dNdy_average(3,3:3*n_nodes:3)=dNbardy(1:n_nodes,3)
        ndelete=0
        do kint = 1, n_points

        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        dof_2d=reshape(dof_total,(/3,length_coord_array/3/))
        dof_inc2d = reshape(dof_increment,(/3,length_coord_array/3/))
        delta_F(1:3,1:3)=matmul(dof_inc2d(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))

        mid_F=matmul(0.5d0*dof_inc2d(1:3,1:n_nodes)+dof_2d(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))

        mid_F(1,1)=mid_F(1,1)+1.d0
        mid_F(2,2)=mid_F(2,2)+1.d0
        mid_F(3,3)=mid_F(3,3)+1.d0

        call invert_small(mid_F,mid_F_inver,J)
        dNdy(1:n_nodes,1:3)=matmul(dNdx(1:n_nodes,1:3),mid_F_inver(1:3,1:3))
!dNdy(1:n_nodes,1:3)=dNdx(1:n_nodes,1:3)
        delta_L=matmul(delta_F,mid_F_inver)
        delta_L_vol=delta_L(1,1)+delta_L(2,2)+delta_L(3,3)
        delta_L(1,1)=delta_L(1,1)+(delta_eta-delta_L_vol)/3.d0
        delta_L(2,2)=delta_L(2,2)+(delta_eta-delta_L_vol)/3.d0
        delta_L(3,3)=delta_L(3,3)+(delta_eta-delta_L_vol)/3.d0

        dstrain_2d=(delta_L+transpose(delta_L))/2.d0
        dspin=(delta_L-transpose(delta_L))/2.d0

        r=I-dspin/2.d0
        call invert_small(r,r_inv,rdeter)
        drotation=matmul(r,I+dspin/2.d0)


        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)


        Bmodifer=0.d0
        Bmodifer(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)/3.d0
        Bmodifer(1,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)/3.d0
        Bmodifer(1,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)/3.d0
        Bmodifer(2,1:3*n_nodes)=Bmodifer(1,1:3*n_nodes)
        Bmodifer(3,1:3*n_nodes)=Bmodifer(1,1:3*n_nodes)
       B=B-Bmodifer+dNdy_average

  initial_state_variables_currentpoint(1:8)=initial_state_variables(((kint-1)*8+1):kint*8)
!ematrix = initial_state_variables_currentpoint(7)
 !Vf = initial_state_variables_currentpoint(8)
!stress(1:6)=initial_state_variables_currentpoint(1:6)
!dstrain(1)=dstrain_2d(1,1)
!dstrain(2)=dstrain_2d(2,2)
!dstrain(3)=dstrain_2d(3,3)
!dstrain(4)=dstrain_2d(1,2)
!dstrain(5)=dstrain_2d(1,3)
!dstrain(6)=dstrain_2d(2,3)
!dstress = matmul(D,dstrain)
!stress=stress+dstress
!updated_state_variables_currentpoint(1:6)=stress(1:6)
!updated_state_variables_currentpoint(7)=ematrix
!updated_state_variables_currentpoint(8)=vf


      call gurson(element_properties,n_properties,n_state_variables,initial_state_variables_currentpoint, &
              updated_state_variables_currentpoint,dstrain_2d,dRotation,stress,ndelete)

       updated_state_variables(((kint-1)*8+1):kint*8)=updated_state_variables_currentpoint(1:8)

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

    end do
 if (Ndelete ==  n_points) then
    element_deleted = .true.
    else
    element_deleted = .false.
    endif
        return
end subroutine el_fracture_dynamic


!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_fracture_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
  
    integer      :: n_points,kint,k

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node

    real (prec)  ::updated_state_variables_currentpoint(8),initial_state_variables_currentpoint(8)

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



nodal_fieldvariables=0.d0
  
    !     --  Loop over integration points
    do kint = 1, n_points
call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)

updated_state_variables_currentpoint(1:8)=updated_state_variables(((kint-1)*8+1):kint*8)
stress(1:6)=updated_state_variables_currentpoint(1:6)


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
            else if (strcmp(field_variable_names(k),'vf',2) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + &
                updated_state_variables_currentpoint(8)*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
 
    end do
  
    return
end subroutine fieldvars_fracture_3dbasic



subroutine gurson(element_properties,n_properties,n_state_variables,initial_state_variables, &
updated_state_variables,dstrain,dRot,stress1,ndelete)
 use Types
 use ParamIO
 use Globals, only : TIME, DTIME
 use Element_Utilities, only :rotatesymvec
 use Element_Utilities, only : invert_small
 implicit none
 integer, intent( in ) :: n_properties
 integer, intent( in ) :: n_state_variables
 integer, intent( inout ) :: ndelete
 real (prec), intent( in ) :: element_properties(n_properties)
 real (prec), intent( in ) :: initial_state_variables(n_state_variables)
 real (prec), intent( in ) :: dstrain(3,3)
 real (prec), intent( in ) :: dRot(3,3)
 real (prec), intent( out ) :: stress1(6)
 real (prec), intent( out ) :: updated_state_variables(n_state_variables)

 real (prec)  :: stress0(6),stress0_2d(3,3),dev_stress0(3,3)
 real (prec)  :: ematrix,dematrix,ematrix1
 real (prec)  :: vf,f,vf1
 real (prec)  :: dev_dstrain(3,3)
 real (prec)  :: predictor_devstress(3,3),predictor_devstress_1d(6),predictor_volstress
 real (prec)  :: fFbar
 real (prec)  :: E,v,Y,eps,m,q1,q2,q3,fn,epsn,sn,fc,fF
 real (prec)  :: phi,phi2,F1,F2,se0,p0,se,p,t_inc,tol,err1
 real (prec)  :: w(2,1),dw(2,1),b(2,1),K(2,2),K_inv(2,2),determinant,ee_inc,epsv_inc
 real (prec)  :: dphi_dse,dphi_dp,ddphi_dse_dse,ddphi_dp_dp
 real (prec)  :: dF1_ddphi_dse,dF1_ddphi_dp,dF2_ddphi_dse,dF2_ddphi_dp
 real (prec)  :: dF1_dphi,dF2_dphi,ddphi_dse_dphi,ddphi_dp_dphi
 real (prec)  ::dse_dee_inc,dp_depsv_inc,dF1_dee_inc,dF1_depsv_inc,dF2_dee_inc,dF2_depsv_inc
integer      ::  maxit,nit


E=element_properties(1)
v=element_properties(2)
Y=element_properties(3)
eps=element_properties(4)
m=element_properties(5)
q1=element_properties(6)
q2=element_properties(7)
q3=element_properties(8)
fn=element_properties(9)
epsn=element_properties(10)
sn=element_properties(11)
fc=element_properties(12)
fF=element_properties(13)

 stress0 = initial_state_variables(1:6) ! Stress at start of increment
 ematrix = initial_state_variables(7)
 Vf = initial_state_variables(8)

 stress0_2d=0.d0
 stress0_2d(1,1)=stress0(1)
 stress0_2d(2,2)=stress0(2)
 stress0_2d(3,3)=stress0(3)
 stress0_2d(1,2)=stress0(4)
 stress0_2d(1,3)=stress0(5)
 stress0_2d(2,3)=stress0(6)
 stress0_2d(2,1)=stress0(4)
 stress0_2d(3,1)=stress0(5)
 stress0_2d(3,2)=stress0(6)


Dev_stress0=stress0_2d
 Dev_stress0(1,1)=Dev_stress0(1,1)-(stress0(1)+stress0(2)+stress0(3))/3.d0
 Dev_stress0(2,2)=Dev_stress0(2,2)-(stress0(1)+stress0(2)+stress0(3))/3.d0
 Dev_stress0(3,3)=Dev_stress0(3,3)-(stress0(1)+stress0(2)+stress0(3))/3.d0


 Dev_dstrain=dstrain
 Dev_dstrain(1,1)=Dev_dstrain(1,1)-(dstrain(1,1)+dstrain(2,2)+dstrain(3,3))/3.d0
 Dev_dstrain(2,2)=Dev_dstrain(2,2)-(dstrain(1,1)+dstrain(2,2)+dstrain(3,3))/3.d0
 Dev_dstrain(3,3)=Dev_dstrain(3,3)-(dstrain(1,1)+dstrain(2,2)+dstrain(3,3))/3.d0
!write(IOW,*) drot

 predictor_devstress=E/(1.d0+v)*Dev_dstrain+matmul(dRot,matmul(Dev_stress0,transpose(dRot)))
 predictor_volstress=(stress0(1)+stress0(2)+stress0(3))/3.d0+E/3.D0/(1-2.d0*v)*&
 (dstrain(1,1)+dstrain(2,2)+dstrain(3,3))

 predictor_devstress_1d(1)=predictor_devstress(1,1)
 predictor_devstress_1d(2)=predictor_devstress(2,2)
 predictor_devstress_1d(3)=predictor_devstress(3,3)
 predictor_devstress_1d(4)=predictor_devstress(1,2)
 predictor_devstress_1d(5)=predictor_devstress(1,3)
 predictor_devstress_1d(6)=predictor_devstress(2,3)

 fFbar=(q1+dsqrt(q1*q1-q3))/q3

 if (Vf<fc) then
 f=Vf
 else if (Vf<fF) then
 f=fc+(fFbar-fc)/(fF-fc)*(Vf-fc)
 else
 f=fFbar
 ndelete=ndelete+1
 end if


se0=dsqrt( dot_product(predictor_devstress_1d(1:3),predictor_devstress_1d(1:3)) +&
  2.d0*dot_product(predictor_devstress_1d(4:6),predictor_devstress_1d(4:6)) )*dsqrt(1.5d0)
p0=predictor_volstress

phi2=se0*se0/Y/Y+2.D0*q1*f*cosh(1.5d0*q2*p0/Y)-(1.d0+q3*f*f)

!write(IOW,*) Phi2

if (phi2<0.00000001d0) then
stress1=predictor_devstress_1d
stress1(1)=stress1(1)+predictor_volstress
stress1(2)=stress1(2)+predictor_volstress
stress1(3)=stress1(3)+predictor_volstress
dematrix=0.d0
else
!write(*,*) F

phi=dsqrt(phi2)

t_inc=Dtime
tol=0.00001d0
maxit=500
w=0.d0
nit=0
ee_inc=w(1,1)
epsv_inc=w(2,1)
err1=1.D0
do while ((err1>tol) .AND.  (nit<maxit))          ! Newton Raphson loop
nit = nit + 1
if (nit==maxit) then
    write(6,*) se0
    write(6,*) p0
    stop
endif

se=se0-1.5d0*E/(1.d0+v)*w(1,1)
p=p0-1.d0/3.d0*E/(1.d0-2*v)*w(2,1)

phi2=se*se/Y/Y+2.D0*q1*f*cosh(1.5d0*q2*p/Y)-(1.d0+q3*f*f)

if (phi2 .lt. 1.d-10) then
     ee_inc = ee_inc/10.d0
     epsv_inc = epsv_inc/10.d0
     w = w/10.d0
     nit = nit+1
     write(IOW,*) ' funny iteration'
     if (nit>maxit) then
         write(6,*) ' Newton iterations failed'
         write(6,*) 'Have  a nice day'
         stop
     endif
     cycle
 end if
phi=dsqrt(phi2)

dphi_dse=1.d0/phi*se/Y/Y
dphi_dp=0.5d0/phi*2.d0*q1*f*sinh(1.5d0*q2*p/Y)*1.5d0*q2/Y

F1=dsqrt(dphi_dse*dphi_dse+2.d0/9.d0*dphi_dp*dphi_dp)*(ee_inc/t_inc/eps)-dphi_dse*(phi**m)
F2=dsqrt(dphi_dse*dphi_dse+2.d0/9.d0*dphi_dp*dphi_dp)*(epsv_inc/t_inc/eps)-dphi_dp*(phi**m)

dF1_ddphi_dse=1.d0/dsqrt(dphi_dse*dphi_dse+2.d0/9.d0*dphi_dp*dphi_dp)*&
(ee_inc/t_inc/eps)*2.d0*dphi_dse-(phi**m)
dF1_ddphi_dp=1.d0/dsqrt(dphi_dse*dphi_dse+2.d0/9.d0*dphi_dp*dphi_dp)*&
(ee_inc/t_inc/eps)*4.d0/9.d0*dphi_dp
dF1_dphi=-1.d0*m*dphi_dse*(phi**(m-1.d0))

dF2_ddphi_dp=1.d0/dsqrt(dphi_dse*dphi_dse+2.d0/9.d0*dphi_dp*dphi_dp)*&
(epsv_inc/t_inc/eps)*4.d0/9.d0*dphi_dp-(phi**m)
dF2_ddphi_dse=1.d0/dsqrt(dphi_dse*dphi_dse+2.d0/9.d0*dphi_dp*dphi_dp)*&
(epsv_inc/t_inc/eps)*2.d0*dphi_dp
dF2_dphi=-1.d0*m*dphi_dp*(phi**(m-1.d0))

ddphi_dse_dphi=-1.d0/phi/phi*se/Y/Y
ddphi_dp_dphi=-0.5d0/phi/phi*2.d0*q1*f*sinh(1.5d0*q2*p/Y)*1.5d0*q2/Y

dse_dee_inc=-1.5d0*E/(1.d0+v)
dp_depsv_inc=-1.d0/3.d0*E/(1.d0-2.d0*v)

ddphi_dse_dse=1.d0/phi/Y/Y
ddphi_dp_dp=0.5d0/phi*2.d0*q1*f*cosh(1.5d0*q2*p/Y)*1.5d0*q2/Y*1.5d0*q2/Y

dF1_dee_inc=dsqrt(dphi_dse*dphi_dse+2.d0/9.d0*dphi_dp*dphi_dp)*(1.d0/t_inc/eps)
dF1_dee_inc=dF1_dee_inc+dF1_dphi*dphi_dse*dse_dee_inc
dF1_dee_inc=dF1_dee_inc+dF1_ddphi_dse*ddphi_dse_dphi*dphi_dse*dse_dee_inc
dF1_dee_inc=dF1_dee_inc+dF1_ddphi_dp*ddphi_dp_dphi*dphi_dse*dse_dee_inc
dF1_dee_inc=dF1_dee_inc+dF1_ddphi_dse*ddphi_dse_dse*dse_dee_inc

dF1_depsv_inc=dF1_dphi*dphi_dp*dp_depsv_inc
dF1_depsv_inc=dF1_depsv_inc+dF1_ddphi_dse*ddphi_dse_dphi*dphi_dp*dp_depsv_inc
dF1_depsv_inc=dF1_depsv_inc+dF1_ddphi_dp*ddphi_dp_dphi*dphi_dp*dp_depsv_inc
dF1_depsv_inc=dF1_depsv_inc+dF1_ddphi_dp*ddphi_dp_dp*dp_depsv_inc

dF2_depsv_inc=dsqrt(dphi_dse*dphi_dse+2.d0/9.d0*dphi_dp*dphi_dp)*(1.d0/t_inc/eps)
dF2_depsv_inc=dF2_depsv_inc+dF2_dphi*dphi_dp*dp_depsv_inc
dF2_depsv_inc=dF2_depsv_inc+dF2_ddphi_dp*ddphi_dp_dphi*dphi_dp*dp_depsv_inc
dF2_depsv_inc=dF2_depsv_inc+dF2_ddphi_dse*ddphi_dse_dphi*dphi_dp*dp_depsv_inc
dF2_depsv_inc=dF2_depsv_inc+dF2_ddphi_dp*ddphi_dp_dp*dp_depsv_inc

dF2_dee_inc=dF2_dphi*dphi_dse*dse_dee_inc
dF2_dee_inc=dF2_dee_inc+dF2_ddphi_dp*ddphi_dp_dphi*dphi_dse*dse_dee_inc
dF2_dee_inc=dF2_dee_inc+dF2_ddphi_dse*ddphi_dse_dphi*dphi_dse*dse_dee_inc
dF2_dee_inc=dF2_dee_inc+dF2_ddphi_dse*ddphi_dse_dse*dse_dee_inc

K(1,1)=dF1_dee_inc
K(1,2)=dF1_depsv_inc
K(2,1)=dF2_dee_inc
K(2,2)=dF2_depsv_inc
b(1,1) = -1.d0*F1
b(2,1)=-1.d0*F2
call invert_small(K,K_inv,determinant)
dw=matmul(K_inv,b)
w = w + dw
err1 = dsqrt((dw(1,1)*dw(1,1)+dw(2,1)*dw(2,1))/(w(1,1)*w(1,1)+w(2,1)*w(2,1)))
!if (w(1,1)>0.d0) then
ee_inc=w(1,1)
!else
!ee_inc=ee_inc/10.d0
!w(1,1) = ee_inc
!nit = nit+1
!     write(*,*) ' funny iteration1'
!
!endif
!if (w(2,1)<0.d0) then
!write(*,*) ' funny iteration3' ,epsv_inc, PHI2
!epsv_inc = epsv_inc/10.d0
!w(2,1)=epsv_inc
!!nit = nit+1
!
!
!
!else
epsv_inc=w(2,1)
!endif
if (nit==maxit) then
    write(6,*) se0
    write(6,*) p0
    stop
endif

end do

!write(IOW,'(//A25,I5,A26)') '    Inner NR Step number ',nit,'         For Inner NR     '
     write(IOW,'(A36,I5)')       '    nit                             ',nit
     write(IOW,'(A36,D12.5)')    '    err                             ',err1
     write(IOW,'(A36,D12.5)')    '    ee_inc                          ',ee_inc
     write(IOW,'(A36,D12.5)')    '    depsv_inc                       ',epsv_inc
     write(IOW,'(A36,D12.5)')    '    se                              ',se
     write(IOW,'(A36,D12.5)')    '    p                               ',p
!    write(IOW,'(A36,D12.5)')    '    PHI2                            ',phi2
!if (phi2<0.00000001d0) then
!stress1=predictor_devstress_1d
!stress1(1)=stress1(1)+predictor_volstress
!stress1(2)=stress1(2)+predictor_volstress
!stress1(3)=stress1(3)+predictor_volstress
!dematrix=0.d0
!else
stress1=predictor_devstress_1d-ee_inc*E/(1.d0+v)*1.5d0*predictor_devstress_1d/se0
stress1(1)=stress1(1)+(p0-E/3.d0/(1.d0-2.d0*v)*epsv_inc)
stress1(2)=stress1(2)+(p0-E/3.d0/(1.d0-2.d0*v)*epsv_inc)
stress1(3)=stress1(3)+(p0-E/3.d0/(1.d0-2.d0*v)*epsv_inc)
dematrix=eps*dtime/(1.D0-vf)*(phi**m)/dsqrt(dphi_dse*dphi_dse+2.d0/9.d0*dphi_dp*dphi_dp)*&
(dphi_dse*se+1.d0/3.d0*dphi_dp*p)

end if


write(IOW,*) ' ematrix ',ematrix

10 vf1=1.d0+(vf-1.d0)*exp(-1.d0*epsv_inc)+fN*dematrix/sn/dsqrt(2.d0*pi)*exp(-1.d0/2.d0*&
   (((ematrix-epsn)/sn)**2.D0))
ematrix1=ematrix+dematrix
updated_state_variables(1:6)=stress1(1:6)
updated_state_variables(7)=ematrix1
updated_state_variables(8)=vf1
write(IOW,*) 'vf1',vf1
!write(IOW,'(A36,D12.5)')    '    ematrix1                       ',ematrix1
!write(IOW,'(A36,D12.5)')    '    dematrix                       ',dematrix
write(*,*) nit

return

 end subroutine gurson

