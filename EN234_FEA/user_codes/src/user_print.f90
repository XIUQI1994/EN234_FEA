subroutine user_print(n_steps)
  use Types
  use ParamIO
  use Globals, only : TIME, DTIME
!  use Mesh
  use Printparameters, only : n_user_print_files                  ! No. files specified by the user
  use Printparameters, only : n_user_print_parameters             ! No. user supplied print parameters
  use Printparameters, only : user_print_units                    ! Unit numbers
  use Printparameters, only : user_print_parameters               ! List of user supplied parameters
  use User_Subroutine_Storage, only : length_state_variable_array ! Max no. state variables on any element
  implicit none
  
  integer, intent(in) :: n_steps                                 ! Current step number
  
  integer ::  lmn,start_element,end_element
  integer ::  status
  integer ::  n_state_vars_per_intpt                                         ! No. state variables per integration point
  real (prec) ::  J_integral_value,J_integral_value_area,Radius
  real (prec) ::   vol_averaged_strain(6),vol_averaged_stress(6)                                  ! Volume averaged strain in an element
  real (prec), allocatable ::   vol_averaged_state_variables(:)              ! Volume averaged state variables in an element

write(user_print_units(1),'(A36,I5)')'n_steps',n_steps
  write(user_print_units(1),*)'time',time+dtime

!
!  Use this file to process or print time histories of the solution, or to print a non-standard mesh.
!
!  As an example, this subroutine computes the volume averaged infinitesimal strain and the volume average
!  element state variables (if they exist) in an element.   The element is specified by the user.
!  The first six state variables (which are usually the stresses) are printed along with the strains.
!
!
!
!start_element = int(user_print_parameters(1))     ! The element number
!  end_element= int(user_print_parameters(2))
!  Radius=user_print_parameters(3)
!  do lmn =  start_element, end_element
!call compute_J_integral(lmn,Radius,J_integral_value)
!  J_integral_value_area=J_integral_value_area+J_integral_value
!end do
!  write(user_print_units(1),*) J_integral_value_area
do lmn=1,2
call hybrid_print(lmn,n_steps)
end do
!   allocate(vol_averaged_state_variables(length_state_variable_array), stat=status)
!
!   if (status/=0) then
!      write(IOW,*) ' Error in subroutine user_print'
!      write(IOW,*) ' Unable to allocate memory for state variables '
!      stop
!   endif
!
!   lmn = int(user_print_parameters(1))     ! The element number
!
!
!   call compute_element_volume_average_3D(lmn,vol_averaged_strain,vol_averaged_state_variables,length_state_variable_array, &
!                                                       n_state_vars_per_intpt)
!   call compute_element_stress_average_3D(lmn,vol_averaged_stress,vol_averaged_state_variables,length_state_variable_array, &
!                                                       n_state_vars_per_intpt)
!write(user_print_units(1),*) TIME+DTIME,vol_averaged_strain(1)
!write(user_print_units(2),*) TIME+DTIME,vol_averaged_stress(1)


!write(user_print_units(1),*) vol_averaged_strain(1)
!write(user_print_units(2),*) vol_averaged_stress(1)

!    if (TIME<1.d-12) then
!      if (n_state_vars_per_intpt<6) then
!        write(user_print_units(1),'(A)') 'VARIABLES = TIME,e11,e22,e33,e12,e13,e23'
!      else
!         write(user_print_units(1),'(A)') 'VARIABLES = TIME,e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23'
!      endif
!    endif
!
!   if (n_state_vars_per_intpt<6) then
!      write(user_print_units(1),'(7(1x,D12.5))') TIME+DTIME,vol_averaged_strain(1:6)
!   else
!      vol_averaged_state_variables(1:3) = vol_averaged_state_variables(1:3) + vol_averaged_state_variables(7)
!      write(user_print_units(1),'(13(1x,D12.5))') TIME+DTIME,vol_averaged_strain(1:6),vol_averaged_state_variables(1:6)
!   endif
!

!      allocate storage for field variables (see printing.f90) (i.e. P vector - called field_variables)
!      call your assemble_field_projection (constructs the P vector)
!      call your assemble_lumped_mass_matrix
!       P(n) = P(n)/lumped_mass_matrix(n) for all nodes n  (or solve M*P=P)

!       Print tecplot header
!       Print nodes (coords, DOF, field vars for each node one at a time)
!       Print element connectivity




end subroutine user_print

subroutine compute_element_volume_average_3D(lmn,vol_averaged_strain,vol_averaged_state_vars,length_output_array, &
                                                                                                       n_state_vars_per_intpt)
    use Types
    use ParamIO
    use Mesh, only : extract_element_data
    use Mesh, only : extract_node_data
    use User_Subroutine_Storage
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx=>vol_avg_shape_function_derivatives_3D
    implicit none

    integer, intent ( in )      :: lmn                                          ! Element number
    integer, intent ( in )      :: length_output_array

    real (prec), intent( out )  ::  vol_averaged_strain(6)
    real (prec), intent( out )  ::  vol_averaged_state_vars(length_output_array)

    integer, intent( out )      :: n_state_vars_per_intpt

    ! Local variables

    integer    :: node_identifier                              ! Flag identifying node type
    integer    :: element_identifier                           ! Flag identifying element type (specified in .in file)
    integer    :: n_nodes                                      ! # nodes on the element
    integer    :: n_properties                                 ! # properties for the element
    integer    :: n_state_variables                            ! # state variables for the element


    integer, allocatable    :: node_list(:)                                ! List of nodes on the element (connectivity)

    real( prec ), allocatable   :: element_properties(:)                  ! Element or material properties, stored in order listed in input file
    real( prec ), allocatable   :: initial_state_variables(:)             ! Element state variables at start of step.
    real( prec ), allocatable   :: updated_state_variables(:)             ! State variables at end of time step

    real( prec ), allocatable   :: x(:,:)                                  ! Nodal coords x(i,a) is ith coord of ath node
    real( prec ), allocatable   :: dof_increment(:)                        ! DOF increment, using usual element storage convention
    real( prec ), allocatable   :: dof_total(:)                            ! accumulated DOF, using usual element storage convention

    integer      :: n_points,kint,i
    integer      :: n_coords, n_dof
    integer      :: iof
    integer      :: status

    real (prec)  ::  el_vol
    real (prec), allocatable  ::  B(:,:)               ! strain = B*(dof_total+dof_increment)
    real (prec), allocatable  ::  dNdx_average(:,:),Bmodifer(:,:)
    real (prec)  ::  strain(6),stress(6)               ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  dstrain(6)                        ! Strain increment vector
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    !
    !  Allocate memory to store element data.
    !  The variables specifying the size of the arrays are stored in the module user_subroutine_storage
    !  They are initialized when the input file is read, and specify the size of the arrays required to store data
    !  for any element in the mesh.  Some elements may require less storage.

    allocate(node_list(length_node_array), stat=status)
    allocate(element_properties(length_property_array), stat=status)
    allocate(initial_state_variables(length_state_variable_array), stat=status)
    allocate(updated_state_variables(length_state_variable_array), stat=status)
    allocate(x(3,length_coord_array/3), stat=status)
    allocate(dof_increment(length_dof_array), stat=status)
    allocate(dof_total(length_dof_array), stat=status)
    allocate(B(6,length_dof_array), stat=status)
    allocate(dNdx_average(6,length_dof_array), stat=status)
    allocate(Bmodifer(6,length_dof_array), stat=status)


    if (status/=0) then
       write(IOW,*) ' Error in subroutine compute_volume_average_3D'
       write(IOW,*) ' Unable to allocate memory for element variables '
       stop
    endif
    !
    ! Extract element and node data from global storage (see module Mesh.f90 for the source code for these subroutines)

    call extract_element_data(lmn,element_identifier,n_nodes,node_list,n_properties,element_properties, &
                                            n_state_variables,initial_state_variables,updated_state_variables)

    do i = 1, n_nodes
        iof = 3*(i-1)+1     ! Points to first DOF for the node in the dof_increment and dof_total arrays
        call extract_node_data(node_list(i),node_identifier,n_coords,x(1:3,i),n_dof, &
                                                 dof_increment(iof:iof+2),dof_total(iof:iof+2))
    end do

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    vol_averaged_strain = 0.d0
   ! vol_averaged_stress = 0.d0
    vol_averaged_state_vars = 0.d0
    el_vol = 0.d0
    n_state_vars_per_intpt = n_state_variables/n_points

    if (n_state_vars_per_intpt>size(vol_averaged_state_vars)) then
       write(IOW,*) ' Error detected in subroutine compute_element_volume_average_3d '
       write(IOW,*) ' The element contains ',n_state_vars_per_intpt
       write(IOW,*) ' but the array storing averaged state variables has length ',size(vol_averaged_state_vars)
       stop
    endif
    !     --  Loop over integration points
    dNbardx=0.d0
    do kint=1,n_points
    call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        dNbardx=dNbardx+dNdx*w(kint)*determinant
        el_vol=el_vol+w(kint)*determinant
    end do
        dNbardx=dNbardx/el_vol/3.d0
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

        vol_averaged_strain(1:6) = vol_averaged_strain(1:6) + (strain(1:6)+dstrain(1:6))*w(kint)*determinant

        if (n_state_vars_per_intpt>0) then
           vol_averaged_state_vars(1:n_state_vars_per_intpt) = vol_averaged_state_vars(1:n_state_vars_per_intpt) &
                              + updated_state_variables(iof:iof+n_state_vars_per_intpt-1)*w(kint)*determinant
        endif



    end do

    vol_averaged_strain = vol_averaged_strain/el_vol
    vol_averaged_state_vars = vol_averaged_state_vars/el_vol

    deallocate(node_list)
    deallocate(element_properties)
    deallocate(initial_state_variables)
    deallocate(updated_state_variables)
    deallocate(x)
    deallocate(dof_increment)
    deallocate(dof_total)
    deallocate(B)
    deallocate(dNdx_average)
    deallocate(Bmodifer)

    return




end subroutine compute_element_volume_average_3D




subroutine hybrid_print(lmn,n_steps)
     use Globals, only : TIME, DTIME
    use Types
    use ParamIO
    use Mesh, only : extract_element_data
    use Mesh, only : extract_node_data
    use User_Subroutine_Storage
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx=>vol_avg_shape_function_derivatives_3D
    use Printparameters, only : n_user_print_files                  ! No. files specified by the user
    use Printparameters, only : n_user_print_parameters             ! No. user supplied print parameters
    use Printparameters, only : user_print_units                    ! Unit numbers
    use Printparameters, only : user_print_parameters               ! List of user supplied parameters
    use User_Subroutine_Storage, only : length_state_variable_array ! Max no. state variables on any element

    implicit none

    integer, intent ( in )      :: lmn                                          ! Element number
    integer, intent(in) :: n_steps
!    integer, intent ( in )      :: length_output_array
!
!    real (prec), intent( out )  ::  vol_averaged_strain(6)
!    real (prec), intent( out )  ::  vol_averaged_state_vars(length_output_array)
!
!    integer, intent( out )      :: n_state_vars_per_intpt

    ! Local variables

    integer    :: node_identifier                              ! Flag identifying node type
    integer    :: element_identifier                           ! Flag identifying element type (specified in .in file)
    integer    :: n_nodes                                      ! # nodes on the element
    integer    :: n_properties                                 ! # properties for the element
    integer    :: n_state_variables                            ! # state variables for the element


    integer, allocatable    :: node_list(:)                                ! List of nodes on the element (connectivity)

    real( prec ), allocatable   :: element_properties(:)                  ! Element or material properties, stored in order listed in input file
    real( prec ), allocatable   :: initial_state_variables(:)             ! Element state variables at start of step.
    real( prec ), allocatable   :: updated_state_variables(:)             ! State variables at end of time step

    real( prec ), allocatable   :: x(:,:)                                  ! Nodal coords x(i,a) is ith coord of ath node
    real( prec ), allocatable   :: dof_increment(:)                        ! DOF increment, using usual element storage convention
    real( prec ), allocatable   :: dof_total(:)                            ! accumulated DOF, using usual element storage convention

    integer      :: n_points,kint,i,n_nodes_real
    integer      :: n_coords, n_dof,l
    integer      :: iof
    integer      :: status

    real (prec)  ::  el_vol,xi_g(3)
    real (prec), allocatable  ::  B(:,:)               ! strain = B*(dof_total+dof_increment)
    real (prec), allocatable  ::  dNdx_average(:,:),Bmodifer(:,:)
    real (prec)  ::  strain(6),stress(6)               ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  dstrain(6)                        ! Strain increment vector
    real (prec)  ::  dxidx(3,3), determinant,dof(3)           ! Jacobian inverse and determinant
    !
    !  Allocate memory to store element data.
    !  The variables specifying the size of the arrays are stored in the module user_subroutine_storage
    !  They are initialized when the input file is read, and specify the size of the arrays required to store data
    !  for any element in the mesh.  Some elements may require less storage.

    allocate(node_list(length_node_array), stat=status)
    allocate(element_properties(length_property_array), stat=status)
    allocate(initial_state_variables(length_state_variable_array), stat=status)
    allocate(updated_state_variables(length_state_variable_array), stat=status)
    allocate(x(3,length_coord_array/3), stat=status)
    allocate(dof_increment(length_dof_array), stat=status)
    allocate(dof_total(length_dof_array), stat=status)
    allocate(B(6,length_dof_array), stat=status)
    allocate(dNdx_average(6,length_dof_array), stat=status)
    allocate(Bmodifer(6,length_dof_array), stat=status)


    if (status/=0) then
       write(IOW,*) ' Error in subroutine compute_volume_average_3D'
       write(IOW,*) ' Unable to allocate memory for element variables '
       stop
    endif
    !
    ! Extract element and node data from global storage (see module Mesh.f90 for the source code for these subroutines)

    call extract_element_data(lmn,element_identifier,n_nodes,node_list,n_properties,element_properties, &
                                            n_state_variables,initial_state_variables,updated_state_variables)

    do i = 1, n_nodes
        iof = 3*(i-1)+1     ! Points to first DOF for the node in the dof_increment and dof_total arrays
        call extract_node_data(node_list(i),node_identifier,n_coords,x(1:3,i),n_dof, &
                                                 dof_increment(iof:iof+2),dof_total(iof:iof+2))
    end do

    if (n_nodes == 5) n_nodes_real = 4
    if (n_nodes == 14)n_nodes_real = 10
    if (n_nodes == 9) n_nodes_real = 8
    if (n_nodes == 24) n_nodes_real = 20

   if (n_nodes_real == 4) n_points = 1
    if (n_nodes_real == 10) n_points = 4
    if (n_nodes_real == 8) n_points = 8
    if (n_nodes_real == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes_real, xi, w)


    !     --  Loop over integration points
 write(user_print_units(1),*)'lmn,node,x1,x2,x3,u1,u2,u3'
do l=1,n_nodes_real
dof=dof_total((l-1)*3+1:l*3)+dof_increment((l-1)*3+1:l*3)
write(user_print_units(1),'(I5,I5,D12.5,D12.5,D12.5,D12.5,D12.5,D12.5)')lmn,l,x(1:3,l)+dof,dof
end do
  write(user_print_units(1),*)'lmn,kint,x1,x2,x3,s11,s22,s33,s12,s13,s23'


    do kint=1,n_points
    call calculate_shapefunctions(xi(1:3,kint),n_nodes_real,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes_real),dNdxi(1:n_nodes_real,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes_real,1:3) = matmul(dNdxi(1:n_nodes_real,1:3),dxidx)
       xi_g=matmul(x,N(1:n_nodes_real))
       stress=updated_state_variables((kint-1)*6+1:kint*6)


     write(user_print_units(1),'(I5,I5,D12.5,D12.5,D12.5,D12.5,D12.5,D12.5,D12.5,D12.5,D12.5)')lmn,kint,xi_g,stress
    end do

    deallocate(node_list)
    deallocate(element_properties)
    deallocate(initial_state_variables)
    deallocate(updated_state_variables)
    deallocate(x)
    deallocate(dof_increment)
    deallocate(dof_total)
    deallocate(B)
    deallocate(dNdx_average)
    deallocate(Bmodifer)

    return




end subroutine  hybrid_print







subroutine compute_element_stress_average_3D(lmn,vol_averaged_stress,vol_averaged_state_vars,length_output_array, &
                                                                                                       n_state_vars_per_intpt)
    use Types
    use ParamIO
    use Mesh, only : extract_element_data
    use Mesh, only : extract_node_data
    use User_Subroutine_Storage
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx=>vol_avg_shape_function_derivatives_3D
    implicit none

    integer, intent ( in )      :: lmn                                          ! Element number
    integer, intent ( in )      :: length_output_array

    real (prec), intent( out )  ::  vol_averaged_stress(6)
    real (prec), intent( out )  ::  vol_averaged_state_vars(length_output_array)

    integer, intent( out )      :: n_state_vars_per_intpt

    ! Local variables

    integer    :: node_identifier                              ! Flag identifying node type
    integer    :: element_identifier                           ! Flag identifying element type (specified in .in file)
    integer    :: n_nodes                                      ! # nodes on the element
    integer    :: n_properties                                 ! # properties for the element
    integer    :: n_state_variables                            ! # state variables for the element


    integer, allocatable    :: node_list(:)                                ! List of nodes on the element (connectivity)

    real( prec ), allocatable   :: element_properties(:)                  ! Element or material properties, stored in order listed in input file
    real( prec ), allocatable   :: initial_state_variables(:)             ! Element state variables at start of step.
    real( prec ), allocatable   :: updated_state_variables(:)             ! State variables at end of time step

    real( prec ), allocatable   :: x(:,:)                                  ! Nodal coords x(i,a) is ith coord of ath node
    real( prec ), allocatable   :: dof_increment(:)                        ! DOF increment, using usual element storage convention
    real( prec ), allocatable   :: dof_total(:)                            ! accumulated DOF, using usual element storage convention

    integer      :: n_points,kint,i
    integer      :: n_coords, n_dof
    integer      :: iof
    integer      :: status

    real (prec)  ::  el_vol
    real (prec), allocatable  ::  B(:,:)               ! strain = B*(dof_total+dof_increment)
    real (prec), allocatable  ::  dNdx_average(:,:),Bmodifer(:,:)
    real (prec)  ::  strain(6),stress(6)               ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  dstrain(6)                        ! Strain increment vector
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  :: evol,emises,se,Es,Et,p,smises
    real (prec)  :: edev(6)
    real (prec)  :: s0,e0,n0,k
    real (prec)  ::  D(6,6),D1(6,6),D2(6,6),D3(6,6)
    !
    !  Allocate memory to store element data.
    !  The variables specifying the size of the arrays are stored in the module user_subroutine_storage
    !  They are initialized when the input file is read, and specify the size of the arrays required to store data
    !  for any element in the mesh.  Some elements may require less storage.

    allocate(node_list(length_node_array), stat=status)
    allocate(element_properties(length_property_array), stat=status)
    allocate(initial_state_variables(length_state_variable_array), stat=status)
    allocate(updated_state_variables(length_state_variable_array), stat=status)
    allocate(x(3,length_coord_array/3), stat=status)
    allocate(dof_increment(length_dof_array), stat=status)
    allocate(dof_total(length_dof_array), stat=status)
    allocate(B(6,length_dof_array), stat=status)
    allocate(dNdx_average(6,length_dof_array), stat=status)
    allocate(Bmodifer(6,length_dof_array), stat=status)


    if (status/=0) then
       write(IOW,*) ' Error in subroutine compute_volume_average_3D'
       write(IOW,*) ' Unable to allocate memory for element variables '
       stop
    endif
    !
    ! Extract element and node data from global storage (see module Mesh.f90 for the source code for these subroutines)

    call extract_element_data(lmn,element_identifier,n_nodes,node_list,n_properties,element_properties, &
                                            n_state_variables,initial_state_variables,updated_state_variables)

    do i = 1, n_nodes
        iof = 3*(i-1)+1     ! Points to first DOF for the node in the dof_increment and dof_total arrays
        call extract_node_data(node_list(i),node_identifier,n_coords,x(1:3,i),n_dof, &
                                                 dof_increment(iof:iof+2),dof_total(iof:iof+2))
    end do

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

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

    call initialize_integration_points(n_points, n_nodes, xi, w)

    vol_averaged_stress = 0.d0
   ! vol_averaged_stress = 0.d0
    vol_averaged_state_vars = 0.d0
    el_vol = 0.d0
    n_state_vars_per_intpt = n_state_variables/n_points

    if (n_state_vars_per_intpt>size(vol_averaged_state_vars)) then
       write(IOW,*) ' Error detected in subroutine compute_element_volume_average_3d '
       write(IOW,*) ' The element contains ',n_state_vars_per_intpt
       write(IOW,*) ' but the array storing averaged state variables has length ',size(vol_averaged_state_vars)
       stop
    endif
    !     --  Loop over integration points
    dNbardx=0.d0
    do kint=1,n_points
    call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        dNbardx=dNbardx+dNdx*w(kint)*determinant
        el_vol=el_vol+w(kint)*determinant
    end do
        dNbardx=dNbardx/el_vol/3.d0
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
        evol=0.d0
        edev=0.d0
        emises=0.d0
        se=0.d0

        strain=strain+dstrain
        evol = sum(strain(1:3))/3.d0
        edev = strain
        edev(1:3) = edev(1:3)-evol
        edev(4:6)=edev(4:6)/2.d0
        emises = dsqrt( dot_product(edev(1:3),edev(1:3)) + 2.d0*dot_product(edev(4:6),edev(4:6)) )*dsqrt(1/1.5d0)
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

        vol_averaged_stress(1:6) = vol_averaged_stress(1:6) + stress(1:6)*w(kint)*determinant

        if (n_state_vars_per_intpt>0) then
           vol_averaged_state_vars(1:n_state_vars_per_intpt) = vol_averaged_state_vars(1:n_state_vars_per_intpt) &
                              + updated_state_variables(iof:iof+n_state_vars_per_intpt-1)*w(kint)*determinant
        endif



    end do

    vol_averaged_stress = vol_averaged_stress/el_vol
    vol_averaged_state_vars = vol_averaged_state_vars/el_vol

    deallocate(node_list)
    deallocate(element_properties)
    deallocate(initial_state_variables)
    deallocate(updated_state_variables)
    deallocate(x)
    deallocate(dof_increment)
    deallocate(dof_total)
    deallocate(B)
    deallocate(dNdx_average)
    deallocate(Bmodifer)

    return


end subroutine compute_element_stress_average_3D







subroutine compute_J_integral(lmn,Radius,J_integral_value)
    use Types
    use ParamIO
    use Mesh, only : extract_element_data
    use Mesh, only : extract_node_data
    use Mesh, only : zone,zone_list
    use User_Subroutine_Storage
    use Element_Utilities, only : N => shape_functions_2D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_2D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_2D
    use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
    use Element_Utilities, only : dxdxi => jacobian_2D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

!    integer, intent ( in )      :: area_integral_option                                    ! Element number
    integer, intent ( in )      :: lmn
    real (prec), intent( in )   :: Radius
!    real (prec), intent( in )   :: Height

    real (prec), intent( out )  ::  J_integral_value


    ! Local variables

    integer    :: node_identifier                              ! Flag identifying node type
    integer    :: element_identifier                           ! Flag identifying element type (specified in .in file)
    integer    :: n_nodes                                      ! # nodes on the element
    integer    :: n_properties                                 ! # properties for the element
    integer    :: n_state_variables                            ! # state variables for the element
    integer    :: n_coords                                     ! No. coords for a node
    integer    :: n_dof                                        ! No. DOFs for a node

    integer      :: status
    integer      :: iof
!    integer      :: lmn               ! Element number
    integer      :: i                 ! Loop counter

    integer      :: n_points,kint,k
    real (prec)  ::  strain(3), dstrain(3)             ! Strain vector contains [e11, e22, 2e12]
    real (prec)  ::  stress(3)                         ! Stress vector contains [s11, s22, s12]
    real (prec)  ::  D(3,3)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
    real (prec)  :: E, xnu, D33, D11, D12
    real (prec)  :: xi1,xi2,du1dx2
!   The arrays below have to be given dimensions large enough to store the data. It doesnt matter if they are too large.

    integer, allocatable    :: node_list(:)                                ! List of nodes on the element (connectivity)

    real( prec ), allocatable   :: element_properties(:)                  ! Element or material properties, stored in order listed in input file
    real( prec ), allocatable   :: initial_state_variables(:)             ! Element state variables at start of step.
    real( prec ), allocatable   :: updated_state_variables(:)             ! State variables at end of time step

    real( prec ), allocatable   :: x(:,:)                                  ! Nodal coords x(i,a) is ith coord of ath node
    real( prec ), allocatable   :: dof_increment(:)                        ! DOF increment, using usual element storage convention
    real( prec ), allocatable   :: dof_total(:)                            ! accumulated DOF, using usual element storage convention
    real( prec ), allocatable   :: nodal_q(:)                              ! Nodal values of q

    real (prec), allocatable  ::  B(:,:)                                   ! strain = B*(dof_total+dof_increment)
    real (prec), allocatable  ::  dndx2(:)
    real (prec), allocatable  ::  dofu1(:)
    real (prec), allocatable  ::  xx1(:)
    real (prec), allocatable  ::  xx2(:)
    !
    !  The variables specifying the sizes of the arrays listed below are determined while reading the input file
    !  They specify the array dimensions required to store the relevant variables for _any_ element or node in the mesh
    !  The actual data will vary depending on the element or node selected
    !
    allocate(node_list(length_node_array), stat=status)
    allocate(element_properties(length_property_array), stat=status)
    allocate(initial_state_variables(length_state_variable_array), stat=status)
    allocate(updated_state_variables(length_state_variable_array), stat=status)
    allocate(x(2,length_coord_array/2), stat=status)
    allocate(dof_increment(length_dof_array), stat=status)
    allocate(dof_total(length_dof_array), stat=status)
    allocate(dofu1(length_dof_array/2), stat=status)
    allocate(nodal_q(length_node_array), stat=status)
    allocate(B(3,length_dof_array), stat=status)
    allocate(dndx2(length_dof_array/2), stat=status)
     allocate(xx2(length_dof_array/2), stat=status)
      allocate(xx1(length_dof_array/2), stat=status)

  !  Write your code to calculate the J integral here

!    start_element = zone_list(2)%start_element

  !  The two subroutines below extract data for elements and nodes (see module Mesh.f90 for the source code for these subroutines)

    call extract_element_data(lmn,element_identifier,n_nodes,node_list,n_properties,element_properties, &
                                            n_state_variables,initial_state_variables,updated_state_variables)

    do i = 1, n_nodes
        iof = 2*(i-1)+1     ! Points to first DOF for the node in the dof_increment and dof_total arrays
        call extract_node_data(node_list(i),node_identifier,n_coords,x(1:2,i),n_dof, &
                                                 dof_increment(iof:iof+2),dof_total(iof:iof+2))
    end do


    if (n_nodes == 3) n_points = 1
    if (n_nodes == 4) n_points = 4
    if (n_nodes == 6) n_points = 4
    if (n_nodes == 8) n_points = 9
    if (n_nodes == 9) n_points = 9

    call initialize_integration_points(n_points, n_nodes, xi, w)


 D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d33 = 0.5D0*E/(1+xnu)
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1,2) = d12
    D(2,1) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d33
      J_integral_value=0

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))

        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
        xi1=0.d0
        xx1(1:n_nodes)=x(1,1:n_nodes)
        xx2(1:n_nodes)=x(2,1:n_nodes)
        xi1=dot_product(xx1,N)

        xi2=0.d0
        Xi2=dot_product(xx2,N)


        B = 0.d0
        B(1,1:2*n_nodes-1:2) = dNdx(1:n_nodes,1)
        B(2,2:2*n_nodes:2) = dNdx(1:n_nodes,2)
        B(3,1:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
        B(3,2:2*n_nodes:2) = dNdx(1:n_nodes,1)
        dndx2=0.d0
        dofu1(1:n_nodes)=dof_total(1:2*n_nodes-1:2)+dof_increment(1:2*n_nodes-1:2)

         dndx2(1:n_nodes) = dNdx(1:n_nodes,2)
         du1dx2=dot_product(dofu1,dndx2)





        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
         strain=strain+dstrain

        stress = matmul(D,strain)
!        J_integral_value= J_integral_value-(1/0.0006/dsqrt(xi1*xi1+xi2*xi2))


        J_integral_value= J_integral_value-((1./Radius/dsqrt(xi1*xi1+xi2*xi2))*&
        (stress(3)*du1dx2*xi2+stress(1)*du1dx2*xi1+stress(3)*strain(2)*xi1+stress(2)*&
        strain(2)*xi2-1/2.d0*dot_product(stress,strain)*xi2)*w(kint)*determinant)

!   write(IOW,*) lmn,J_integral_value

    end do



    deallocate(node_list)
    deallocate(element_properties)
    deallocate(initial_state_variables)
    deallocate(updated_state_variables)
    deallocate(x)
    deallocate(dof_increment)
    deallocate(dof_total)
    deallocate(nodal_q)
    deallocate(B)


    return




end subroutine compute_J_integral
!the J value I get is  1.2286213310963479E-003
