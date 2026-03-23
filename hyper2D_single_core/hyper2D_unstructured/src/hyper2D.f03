! ######################################################################
! ######################################################################
! ###########     _                           ____  ____      ##########   
! ###########    | |__  _   _ _ __   ___ _ __(___ \|  _ \     ##########
! ###########    | '_ \| | | | '_ \ / _ \ '__| __) | | | |    ##########
! ###########    | | | | |_| | |_) |  __/ |   / __/| |_| |    ##########
! ###########    |_| |_|\__, | .__/ \___|_|  /_____|____/     ##########
! ###########           |___/|_|                              ##########
! ###########                                                 ########## 
! ######################################################################
! ######################################################################

program hyper2D

   use global_module   ! Simulation parameters: domain size, number of cells etc
   use grid            ! Loads and stores the grid
   use pde             ! Definition of the system of equations
   use integration     ! Functions for integrating in time and numerical fluxes
   use tools           ! Output subroutines etc etc
   use initialization

   implicit none

   ! The solution is a 3D matrix. 
   ! The first index "eqID" represents the equation 
   ! (1: density, 2: x-momentum, 3: y-momentum, 4: total energy)
   ! these are repeated for each species in the multifluid formulation
   ! the second index represents the cell index.

   real(kind=8), dimension(:,:), allocatable :: U, U_new

   integer      :: t_ID ! Variables for time integration
   real(kind=8) :: dt, t_now, CFL_now


   CALL READINPUT()

   write(*,*) "Initializing solution..."
   allocate(U(N_SPECIES_FLUID*Neq,NCELLS), U_new(N_SPECIES_FLUID*Neq,NCELLS))
   call initialize_solution(U) ! See the pde.f03 module
   
   write(*,*) "Writing solution at time step", 0, "..."
   !call export_sol_vtk(0, U)
   CALL GRID_SAVE(0, 0.d0, U)

   ! $$$$$$$$$$$ Integrate in time $$$$$$$$$$$$$

   ! Nt = ceiling(t_end/dt) ! t_end and dt are defined in global_module.f90
   ! do t_ID = 1, Nt

   t_now = 0.0d0   ! Init
   dt    = 1.0d-10  ! Very first time step
   t_ID  = 0.0     ! Init
   
   do while( t_now .le. t_end)

      ! ------ Prepare variables ------
      t_ID  = t_ID + 1
      t_now = t_now + dt

      invdt_adv = 0.0 ! Global variable, init to zero
      invdt_cond = 0.0
      invdt_diff = 0.0
      invdt_coll = 0.0

      ! ------ Integrate by dt ------
      call forward_Euler_step(U, U_new, dt)

      ! ----- Write solution to VTK, every ... time steps -----
      if ( (t_ID .GE. DUMP_GRID_START) .AND. (mod(t_ID, DUMP_GRID_EVERY) .EQ. 0) ) then 
         write(*,*) "Writing solution at time step", t_ID, "..."
         !call export_sol_vtk(t_ID, U)
         CALL GRID_SAVE(t_ID, t_now, U)
      end if

      ! ------ Estimate current Courant number and update time step -----
      CFL_now = MAX(invdt_adv,invdt_cond, invdt_diff, invdt_coll)*dt
      IF ( mod(t_ID, STATS_EVERY) .EQ. 0 ) THEN
         write(*,'(A,EN15.5,A,F10.5,A,ES14.7,A,A,F10.5,A,F10.5,A,F10.5)') 'Time', t_now, &
         ' [s]. Current CFL: ', CFL_now, '. dt = ', dt, '[s]', &
         ' CFL Advection: ', invdt_adv*dt, ' CFL Conduction: ', invdt_cond*dt, &
         ' CFL Diffusion: ', invdt_diff*dt, ' CFL Collisions: ', invdt_coll*dt
      END IF
      dt = MIN(dt*CFL_target/CFL_now, dtmax)

   end do

end program
