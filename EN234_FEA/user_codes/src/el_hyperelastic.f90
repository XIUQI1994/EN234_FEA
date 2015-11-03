!     Subroutines for basic 3D linear elastic elements 



!==========================SUBROUTINE el_linelast_3dbasic ==============================
subroutine el_hyperelastic_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
    integer      :: n_points,kint,q

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  kirchhoffstress(6),cauchystress(6), kirchhoffstress_2d(3,3),cauchystress_2d(3,3)                        ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array),B1(9,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dNdx_average(6,length_dof_array),Bmodifer(6,length_dof_array)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  dofdis(3,length_dof_array/3)
    real (prec)  :: E, xnu, D44, D11, D12 ,u1,K1             ! Material properties
    real (prec)  :: vol
    real (prec)  ::  Bdefor(3,3),Bdefor_1d(6),Bdefor_INVER(3,3),Bdefor_INVER_1d(6),bdeterminat
    real (prec)  :: F(3,3),F_inver(3,3),J
    real (prec)  :: I0(6,6),I(6),G(6,9)
    real (prec)  :: dNdy(length_dof_array/3,3)
    real (prec)  :: IXI(6,6),BXBINV(6,6),IXBINV(6,6)
    real (prec)  ::  s(3,length_dof_array/3),pvec(length_dof_array),svec(length_dof_array)
    real (prec)  ::smat(length_dof_array,length_dof_array)
    real (prec)  ::pmat(length_dof_array,length_dof_array),sigma(length_dof_array,length_dof_array)

    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    I0=0.D0
    I0(1,1)=1.D0
    I0(2,2)=1.D0
    I0(3,3)=1.D0
    I0(4,4)=1.D0/2.D0
    I0(5,5)=1.D0/2.D0
    I0(6,6)=1.D0/2.D0
    I=0.D0
    I(1:3)=1.D0

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
	
    D = 0.d0
    U1 = element_properties(1)
    K1 = element_properties(2)
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
!    dNbardx=0.d0
!    vol=0.d0
!    do kint=1,n_points
!    call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
!        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
!        call invert_small(dxdxi,dxidx,determinant)
!        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
!        dNbardx=dNbardx+dNdx*w(kint)*determinant
!        vol=vol+w(kint)*determinant
!    end do
!        dNbardx=dNbardx/vol/3.d0
!        dNdx_average=0.d0
!        dNdx_average(1,1:3*n_nodes-2:3)=dNbardx(1:n_nodes,1)
!        dNdx_average(1,2:3*n_nodes-1:3)=dNbardx(1:n_nodes,2)
!        dNdx_average(1,3:3*n_nodes:3)=dNbardx(1:n_nodes,3)
!        dNdx_average(2,1:3*n_nodes-2:3)=dNbardx(1:n_nodes,1)
!        dNdx_average(2,2:3*n_nodes-1:3)=dNbardx(1:n_nodes,2)
!        dNdx_average(2,3:3*n_nodes:3)=dNbardx(1:n_nodes,3)
!        dNdx_average(3,1:3*n_nodes-2:3)=dNbardx(1:n_nodes,1)
!        dNdx_average(3,2:3*n_nodes-1:3)=dNbardx(1:n_nodes,2)
!        dNdx_average(3,3:3*n_nodes:3)=dNbardx(1:n_nodes,3)

    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
!        B = 0.d0
!        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
!        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
!        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
!        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
!        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
!        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
!        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
!        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
!        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)
!


!        if ( element_identifier ==1002) then
!        Bmodifer=0.d0
!        Bmodifer(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)/3.d0
!        Bmodifer(1,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)/3.d0
!        Bmodifer(1,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)/3.d0
!        Bmodifer(2,1:3*n_nodes)=Bmodifer(1,1:3*n_nodes)
!        Bmodifer(3,1:3*n_nodes)=Bmodifer(1,1:3*n_nodes)
!
!        B=B-Bmodifer+dNdx_average
!        end if

        !calulation of B dydx


        dofdis = reshape(dof_total+dof_increment,(/3,length_coord_array/3/))
        F(1:3,1:3)=matmul(dofdis(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
        F(1,1)=F(1,1)+1
        F(2,2)=F(2,2)+1
        F(3,3)=F(3,3)+1
      !  write(IOW,*)F
        Bdefor=matmul(F(1:3,1:3),transpose(F(1:3,1:3)))
        Bdefor_1d(1)=Bdefor(1,1)
        Bdefor_1d(2)=Bdefor(2,2)
        Bdefor_1d(3)=Bdefor(3,3)
        Bdefor_1d(4)=Bdefor(1,2)
        Bdefor_1d(5)=Bdefor(1,3)
        Bdefor_1d(6)=Bdefor(2,3)
        call invert_small(Bdefor,Bdefor_inver,Bdeterminat)
        Bdefor_inver_1d(1)=Bdefor_inver(1,1)
        Bdefor_inver_1d(2)=Bdefor_inver(2,2)
        Bdefor_inver_1d(3)=Bdefor_inver(3,3)
        Bdefor_inver_1d(4)=Bdefor_inver(1,2)
        Bdefor_inver_1d(5)=Bdefor_inver(1,3)
        Bdefor_inver_1d(6)=Bdefor_inver(2,3)

        call invert_small(F,F_inver,J)
        dNdy(1:n_nodes,1:3)=matmul(dNdx(1:n_nodes,1:3),F_inver(1:3,1:3))
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

        B1 = 0.d0
        B1(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B1(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B1(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B1(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B1(5,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B1(6,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B1(7,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B1(8,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B1(9,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)
!find cauchy stress
        cauchystress=Bdefor_1d
        cauchystress(1:3)= cauchystress(1:3)-sum(Bdefor_1d(1:3))/3.d0
        cauchystress= cauchystress*u1/J**(5.D0/3.D0)
        cauchystress(1:3)= cauchystress(1:3)+K1*(J-1)
        cauchystress_2d(1,1)=cauchystress(1)
        cauchystress_2d(2,2)=cauchystress(2)
        cauchystress_2d(3,3)=cauchystress(3)
        cauchystress_2d(1,2)=cauchystress(4)
        cauchystress_2d(1,3)=cauchystress(5)
        cauchystress_2d(2,3)=cauchystress(6)
        cauchystress_2d(2,1)=cauchystress(4)
        cauchystress_2d(3,1)=cauchystress(5)
        cauchystress_2d(3,2)=cauchystress(6)

!transform to kirchstress
 !      kirchhoffstress_2d(1:3,1:3)=J*matmul(F_inver(1:3,1:3),cauchystress_2d(1:3,1:3))
       kirchhoffstress_2d(1:3,1:3)=J*cauchystress_2d(1:3,1:3)
       kirchhoffstress=0.d0
       kirchhoffstress(1)=kirchhoffstress_2d(1,1)
       kirchhoffstress(2)=kirchhoffstress_2d(2,2)
       kirchhoffstress(3)=kirchhoffstress_2d(3,3)
       kirchhoffstress(4)=kirchhoffstress_2d(1,2)
       kirchhoffstress(5)=kirchhoffstress_2d(1,3)
       kirchhoffstress(6)=kirchhoffstress_2d(2,3)

!define G
       G=0.D0
       G(1,1)=Bdefor(1,1)
       G(1,4)=Bdefor(1,2)
       G(1,6)=Bdefor(1,3)
       G(2,2)=Bdefor(2,2)
       G(2,5)=Bdefor(1,2)
       G(2,8)=Bdefor(2,3)
       G(3,3)=Bdefor(3,3)
       G(3,7)=Bdefor(1,3)
       G(3,9)=Bdefor(1,3)
       G(4,1)=Bdefor(1,2)
       G(4,2)=Bdefor(1,2)
       G(4,4)=Bdefor(2,2)
       G(4,5)=Bdefor(1,1)
       G(4,6)=Bdefor(2,3)
       G(4,8)=Bdefor(1,3)
       G(5,1)=Bdefor(1,3)
       G(5,3)=Bdefor(1,3)
       G(5,4)=Bdefor(2,3)
       G(5,6)=Bdefor(3,3)
       G(5,7)=Bdefor(1,1)
       G(5,9)=Bdefor(1,2)
       G(6,2)=Bdefor(2,3)
       G(6,3)=Bdefor(2,3)
       G(6,5)=Bdefor(1,3)
       G(6,7)=Bdefor(1,2)
       G(6,8)=Bdefor(3,3)
       G(6,9)=Bdefor(2,2)
       G=G*2.D0



 !parts of D

IXI = spread(I,dim=2,ncopies=6)*spread(I,dim=1,ncopies=6)
BXBINV=spread(Bdefor_1d,dim=2,ncopies=6)*spread(Bdefor_inver_1d,dim=1,ncopies=6)
IXBINV=spread(I,dim=2,ncopies=6)*spread(Bdefor_inver_1d,dim=1,ncopies=6)
D=I0+sum(Bdefor_1d(1:3))*IXBINV/9.d0-IXI/3.d0-BXBINV/3.d0
D=D*U1/(J**(2.D0/3.D0))
D=D+K1*J*(J-0.5D0)*IXBINV
!write(IOW,*)matmul(G(1:6,1:9),B1(1:9,1:3*n_nodes))

S = reshape(matmul(transpose(B),kirchhoffstress),(/3,length_dof_array/3/))
do q = 1,n_nodes
Pvec = reshape(spread(transpose(dNdx(q:q,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
Pmat(3*q-2:3*q,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
Svec = reshape(spread(S(1:3,q:q),dim=2,ncopies=n_nodes),(/3*n_nodes/))
Smat(3*q-2:3*q,1:3*n_nodes) = spread(Svec,dim=1,ncopies=3)
end do
Sigma = Pmat*transpose(Smat)

!define stiffness


!        strain = matmul(B,dof_total)
!        dstrain = matmul(B,dof_increment)
!
!        stress = matmul(D,strain+dstrain)
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),kirchhoffstress)*w(kint)*determinant

        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + (-Sigma+matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D(1:6,1:6),matmul(G(1:6,1:9),B1(1:9,1:3*n_nodes)))))&
            *w(kint)*determinant

    end do
  
    return
end subroutine el_hyperelastic_3dbasic


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_hyperelastic_3dbasic_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
end subroutine el_hyperelastic_3dbasic_dynamic


!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_hyperelastic_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
    use Element_Utilities, only :principalvals33
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
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dNdx_average(6,length_dof_array),Bmodifer(6,length_dof_array)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    real (prec)  :: p, smises                          ! Pressure and Mises stress
    real (prec)  :: vol
    real (prec)  ::  dofdis(3,length_dof_array/3)
    real (prec)  :: u1,K1             ! Material properties
    real (prec)  ::  Bdefor(3,3),Bdefor_1d(6),Bdefor_INVER(3,3),Bdefor_INVER_1d(6),bdeterminat
    real (prec)  :: F(3,3),F_inver(3,3),J
    real (prec)  :: I0(6,6),I(6),G(6,9)
    real (prec)  :: dNdy(length_dof_array/3,3)
    real (prec)  ::  kirchhoffstress(6),cauchystress(6), kirchhoffstress_2d(3,3),cauchystress_2d(3,3)
    real (prec)  :: principalstress(3)
    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    U1 = element_properties(1)
    K1 = element_properties(2)

    nodal_fieldvariables = 0.d0
	
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
!
    !     --  Loop over integration points
!    dNbardx=0.d0
!    vol=0.d0
       do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
!        B = 0.d0
!        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
!        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
!        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
!        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
!        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
!        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
!        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
!        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
!        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)
!


!        if ( element_identifier ==1002) then
!        Bmodifer=0.d0
!        Bmodifer(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)/3.d0
!        Bmodifer(1,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)/3.d0
!        Bmodifer(1,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)/3.d0
!        Bmodifer(2,1:3*n_nodes)=Bmodifer(1,1:3*n_nodes)
!        Bmodifer(3,1:3*n_nodes)=Bmodifer(1,1:3*n_nodes)
!
!        B=B-Bmodifer+dNdx_average
!        end if

        !calulation of B dydx


        dofdis = reshape(dof_total+dof_increment,(/3,length_coord_array/3/))
        F(1:3,1:3)=matmul(dofdis(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
        F(1,1)=F(1,1)+1
        F(2,2)=F(2,2)+1
        F(3,3)=F(3,3)+1
      !  write(IOW,*)F
        Bdefor=matmul(F(1:3,1:3),transpose(F(1:3,1:3)))
        Bdefor_1d(1)=Bdefor(1,1)
        Bdefor_1d(2)=Bdefor(2,2)
        Bdefor_1d(3)=Bdefor(3,3)
        Bdefor_1d(4)=Bdefor(1,2)
        Bdefor_1d(5)=Bdefor(1,3)
        Bdefor_1d(6)=Bdefor(2,3)
        call invert_small(Bdefor,Bdefor_inver,Bdeterminat)
        Bdefor_inver_1d(1)=Bdefor_inver(1,1)
        Bdefor_inver_1d(2)=Bdefor_inver(2,2)
        Bdefor_inver_1d(3)=Bdefor_inver(3,3)
        Bdefor_inver_1d(4)=Bdefor_inver(1,2)
        Bdefor_inver_1d(5)=Bdefor_inver(1,3)
        Bdefor_inver_1d(6)=Bdefor_inver(2,3)

        call invert_small(F,F_inver,J)
        dNdy(1:n_nodes,1:3)=matmul(dNdx(1:n_nodes,1:3),F_inver(1:3,1:3))
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

!        B1 = 0.d0
!        B1(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
!        B1(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
!        B1(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
!        B1(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
!        B1(5,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
!        B1(6,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
!        B1(7,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
!        B1(8,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
!        B1(9,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)
!!find cauchy stress
        cauchystress=Bdefor_1d
        cauchystress(1:3)= cauchystress(1:3)-sum(Bdefor_1d(1:3))/3.d0
        cauchystress= cauchystress*u1/J**(5.D0/3.D0)
        cauchystress(1:3)= cauchystress(1:3)+K1*(J-1)
        cauchystress_2d(1,1)=cauchystress(1)
        cauchystress_2d(2,2)=cauchystress(2)
        cauchystress_2d(3,3)=cauchystress(3)
        cauchystress_2d(1,2)=cauchystress(4)
        cauchystress_2d(1,3)=cauchystress(5)
        cauchystress_2d(2,3)=cauchystress(6)
        cauchystress_2d(2,1)=cauchystress(4)
        cauchystress_2d(3,1)=cauchystress(5)
        cauchystress_2d(3,2)=cauchystress(6)

!transform to kirchstress
 !      kirchhoffstress_2d(1:3,1:3)=J*matmul(F_inver(1:3,1:3),cauchystress_2d(1:3,1:3))
       kirchhoffstress_2d(1:3,1:3)=J*cauchystress_2d(1:3,1:3)
       kirchhoffstress=0.d0
       kirchhoffstress(1)=kirchhoffstress_2d(1,1)
       kirchhoffstress(2)=kirchhoffstress_2d(2,2)
       kirchhoffstress(3)=kirchhoffstress_2d(3,3)
       kirchhoffstress(4)=kirchhoffstress_2d(1,2)
       kirchhoffstress(5)=kirchhoffstress_2d(1,3)
       kirchhoffstress(6)=kirchhoffstress_2d(2,3)

        stress=cauchystress
        p = sum(stress(1:3))/3.d0
        sdev = stress
        principalstress=principalvals33(stress)
        principalstress=principalstress
        sdev(1:3) = sdev(1:3)-p
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
            ELSE    if (strcmp(field_variable_names(k),'SP1',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + &
                                        principalstress(1)*N(1:n_nodes)&
             *determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SP2',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + &
                                            principalstress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SP3',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                                           principalstress(3)*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
 
    end do
  
    return
end subroutine fieldvars_hyperelastic_3dbasic

