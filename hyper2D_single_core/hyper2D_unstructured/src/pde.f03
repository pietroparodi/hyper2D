module pde

   use global_module
   use grid

   implicit none

   contains

   ! ============================================================

   subroutine initialize_solution(U)
   
      implicit none

      real(kind=8), dimension(:,:), intent(inout) :: U

      real(kind=8) :: P0, rho, FRAC
      INTEGER :: I, J, SP_ID, eleID

      U = 0.d0

      DO I = 1, MIXTURES(INITIAL_MIX_ID)%N_COMPONENTS
         SP_ID = MIXTURES(INITIAL_MIX_ID)%COMPONENTS(I)%ID
         J = (SP_ID-1)*Neq

         FRAC = MIXTURES(INITIAL_MIX_ID)%COMPONENTS(I)%MOLFRAC

         rho = INITIAL_NRHO*FRAC*SPECIES(SP_ID)%MOLECULAR_MASS
         P0   = INITIAL_NRHO*FRAC*kB*INITIAL_TEMP ! [Pa] gas pressure

         ! Initialize all cells
         DO eleID = 1, NCELLS
            ! Skip cells that are not fluid
            IF (U2D_GRID%CELL_PG(eleID) .NE. -1) THEN
               IF (GRID_BC(U2D_GRID%CELL_PG(eleID))%VOLUME_BC == SOLID) CYCLE
            END IF
            U(J+1,eleID) = rho     ! Density
            U(J+2,eleID) = rho*INITIAL_UX ! Momentum along x
            U(J+3,eleID) = rho*INITIAL_UY ! Momentum along y
            U(J+4,eleID) = rho*(INITIAL_UX**2 + INITIAL_UY**2)/2.0 + P0/(SPECIES(SP_ID)%GAMMA-1.0) ! total energy
         END DO
      END DO


   end subroutine


   ! ============================================================

   subroutine compute_noslip_state(U, Twall, U_wall)

      implicit none

      real(kind=8), dimension(:), intent(in)  :: U
      real(kind=8), dimension(:), intent(out) :: U_wall
      real(kind=8), intent(in)  :: Twall

      real(kind=8), dimension(Neq) :: prim
      INTEGER :: I, FIRST, LAST
      real(kind=8) :: ux, uy, T

      DO I = 1, N_SPECIES_FLUID
         FIRST = (I-1)*Neq+1
         LAST = I*Neq

         call compute_primitive_from_conserved(U(FIRST:LAST), prim, I)

         ux = prim(2)
         uy = prim(3)
         T = prim(4)
         ! Compose new state
         ! Densities are unchanged
         prim(2) = -ux
         prim(3) = -uy
         prim(4) = Twall !2*Twall - T

         call compute_conserved_from_primitive(prim, U_wall(FIRST:LAST), I)

      END DO

   end subroutine 

   ! ============================================================

   subroutine compute_sym_state(U, nx, ny, U_sym)

      implicit none

      real(kind=8), dimension(:), intent(in)  :: U
      real(kind=8), dimension(:), intent(out) :: U_sym
      real(kind=8), intent(in) :: nx, ny

      real(kind=8), dimension(Neq) :: prim

      real(kind=8) :: ux, uy, u_norm

      INTEGER :: I, FIRST, LAST

      DO I = 1, N_SPECIES_FLUID
         FIRST = (I-1)*Neq+1
         LAST = I*Neq

         call compute_primitive_from_conserved(U(FIRST:LAST), prim, I)

         ux = prim(2)
         uy = prim(3)

         ! Mirror the normal velocity component
         ! (this is equivalent to removing the normal component two times)
         u_norm = ux*nx + uy*ny
         ux = ux - 2.0*u_norm*nx
         uy = uy - 2.0*u_norm*ny

         ! Compose new state
         prim(2) = ux
         prim(3) = uy

         call compute_conserved_from_primitive(prim, U_sym(FIRST:LAST), I)

      END DO

   end subroutine 


   ! ============================================================

   subroutine compute_moving_state(U, nx, ny, Twall, ux_wall, uy_wall, U_sym)

      implicit none

      real(kind=8), dimension(:), intent(in)  :: U
      real(kind=8), dimension(:), intent(out) :: U_sym
      real(kind=8), intent(in) :: nx, ny, Twall, ux_wall, uy_wall

      real(kind=8), dimension(Neq) :: prim

      INTEGER :: I, FIRST, LAST

      REAL(KIND=8) :: u_norm, ux_wall_par, uy_wall_par, ux, uy, T

      u_norm = ux_wall*nx + uy_wall*ny
      ux_wall_par = ux_wall - nx*u_norm
      uy_wall_par = uy_wall - ny*u_norm

      DO I = 1, N_SPECIES_FLUID
         FIRST = (I-1)*Neq+1
         LAST = I*Neq

         call compute_primitive_from_conserved(U(FIRST:LAST), prim, I)
         
         ux = prim(2)
         uy = prim(3)
         T = prim(4)

         ! Compose new state
         prim(2) = 2.*ux_wall_par - ux
         prim(3) = 2.*uy_wall_par - uy
         prim(4) = 2*Twall - T

         call compute_conserved_from_primitive(prim, U_sym(FIRST:LAST), I)

      END DO

   end subroutine 

   ! ============================================================

   subroutine compute_primitive_from_conserved(U, prim, SP_ID)

      ! Computes vector of primitive variables "prim" from the conserved variables "U"

      implicit none

      real(kind=8), dimension(Neq), intent(in)  :: U
      real(kind=8), dimension(Neq), intent(out) :: prim
      INTEGER,                      intent(in)  :: SP_ID

      ! Working variables
      real(kind=8) :: rho, ux, uy, P

      ! Extract primitive variables
      rho = U(1)
      ux  = U(2)/(rho + 1.0d-25) ! Use a small tolerance, since rho may be zero
      uy  = U(3)/(rho + 1.0d-25) ! Use a small tolerance, since rho may be zero
      P   = (SPECIES(SP_ID)%GAMMA - 1.0)*( U(4) - rho*(ux**2.0 + uy**2.0)/2.0 )

      ! Compose array of primitive variables
      prim(1) = rho
      prim(2) = ux
      prim(3) = uy 
      prim(4) = P/(kB*rho/SPECIES(SP_ID)%MOLECULAR_MASS) ! T = P/(n*kB) = P/(rho*Ri), with Ri=kB/M the gas constant


   end subroutine

   ! ============================================================

   subroutine compute_conserved_from_primitive(prim, U, SP_ID)

      ! Computes vector of conserved variables "U" from the primitive variables "prim"

      implicit none

      real(kind=8), dimension(Neq), intent(in)  :: prim
      real(kind=8), dimension(Neq), intent(out) :: U
      INTEGER,                      intent(in)  :: SP_ID

      ! Working variables
      real(kind=8) :: rho, ux, uy, T, P

      ! Extract primitive variables
      rho = prim(1)
      ux  = prim(2)
      uy  = prim(3)
      T   = prim(4)

      P = rho/SPECIES(SP_ID)%MOLECULAR_MASS*kB*T

      ! Compose array of conserved variables
      U(1) = rho    ! Density
      U(2) = rho*ux ! Momentum along x
      U(3) = rho*uy ! Momentum along y
      U(4) = rho*(ux**2 + uy**2)/2.0 + P/(SPECIES(SP_ID)%GAMMA-1.0) ! total energy

   end subroutine
   
   ! ============================================================

   subroutine compute_flux_ws(U, F_dot_n, nx, ny, ws_max, ws_min, SP_ID)

      ! Computes the convective flux along x,
      ! and also the maximum and minimum wave speeds (required by some numerical flux schemes)

      implicit none

      real(kind=8), dimension(Neq), intent(in)  :: U
      real(kind=8),                 intent(in)  :: nx, ny
      real(kind=8), dimension(Neq), intent(out) :: F_dot_n
      real(kind=8),                 intent(out) :: ws_max, ws_min
      INTEGER,                      intent(in)  :: SP_ID

      real(kind=8), dimension(Neq) :: prim, Fx, Fy
      real(kind=8) :: rho, ux, uy, T, P, rhoE, u_dot_n
      
      ! Compute primitive variables
      call compute_primitive_from_conserved(U, prim, SP_ID)

      rho = prim(1)
      ux  = prim(2)
      uy  = prim(3)
      T   = prim(4)

      P    = rho*kB/SPECIES(SP_ID)%MOLECULAR_MASS*T ! Compute pressure, P = n*kB*T
      rhoE = rho*(ux*ux + uy*uy)/2.0 + P/(SPECIES(SP_ID)%GAMMA-1)

      ! Assemble fluxes Fx
      Fx(1) = rho*ux
      Fx(2) = rho*ux*ux + P
      Fx(3) = rho*ux*uy
      Fx(4) = rhoE*ux + P*ux
   
      ! Assemble fluxes Fy
      Fy(1) = rho*uy
      Fy(2) = rho*ux*uy
      Fy(3) = rho*uy*uy + P
      Fy(4) = rhoE*uy + P*uy

      ! Rotate in the direction of the normal
      F_dot_n = Fx*nx + Fy*ny

      ! Maximum and minimum wave speeds, normal to the interface (eigenvalues of the Euler system)
      u_dot_n  = ux*nx + uy*ny
      ws_max   = u_dot_n + sqrt(SPECIES(SP_ID)%GAMMA*P/rho)
      ws_min   = u_dot_n - sqrt(SPECIES(SP_ID)%GAMMA*P/rho)

   end subroutine



   subroutine compute_cell_centered_gradients_green_gauss(U, gradU)

      ! This subroutine computes cell-centered gradients using the Green-Gauss
      ! method. This method is very simple, but is inaccurate for skewed cells
      ! or neighboring cells of different size.

      implicit none
   
      real(kind=8), dimension(:,:), intent(in)  :: U
      real(kind=8), dimension(:,:,:), intent(inout)  :: gradU

      INTEGER :: I, J, neigh, FACE_PG, SP_ID, FIRST, LAST, NEIGHBORPG
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Uface
      REAL(KIND=8) :: nx, ny, Aface, Vcell
      REAL(KIND=8), DIMENSION(:), allocatable :: U_adj

      LOGICAL :: ONAXIS, FLUIDBOUNDARY


      ! Compute gradient in each cell from nodal values
      gradU = 0.d0
      ALLOCATE(Uface(N_SPECIES_FLUID*Neq))
      ALLOCATE(U_adj(N_SPECIES_FLUID*Neq))

      DO I = 1, NCELLS
         ! Skip cells that are not fluid
         IF (U2D_GRID%CELL_PG(I) .NE. -1) THEN
            IF (GRID_BC(U2D_GRID%CELL_PG(I))%VOLUME_BC == SOLID) CYCLE
         END IF

         ONAXIS = .FALSE.
         DO J = 1, 3 ! The cell face
            ! Extract data
            Aface = U2D_GRID%CELL_FACES_AREA(J,I)
            nx = U2D_GRID%EDGE_NORMAL(1,J,I)
            ny = U2D_GRID%EDGE_NORMAL(2,J,I)
            Vcell = U2D_GRID%CELL_VOLUMES(I)

            FLUIDBOUNDARY = .FALSE.
            neigh = U2D_GRID%CELL_NEIGHBORS(J,I)
            IF (neigh == -1) THEN
               FLUIDBOUNDARY = .TRUE.
            ELSE
               NEIGHBORPG = U2D_GRID%CELL_PG(neigh)
               IF (NEIGHBORPG .NE. -1) THEN
                  IF (GRID_BC(NEIGHBORPG)%VOLUME_BC == SOLID) FLUIDBOUNDARY = .TRUE.
               END IF
            END IF

            IF (.NOT. FLUIDBOUNDARY) THEN
               Uface = 0.5*(U(:,I) + U(:,neigh))
            else
               FACE_PG = U2D_GRID%CELL_EDGES_PG(J,I)
               if (GRID_BC(FACE_PG)%PARTICLE_BC == STATE) then ! ++++++++ STATE BOUNDARY +++++++++++++++++++
                  DO SP_ID = 1, N_SPECIES_FLUID
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq
                     CALL compute_primitive_from_conserved(GRID_BC(FACE_PG)%U_BOUND(FIRST:LAST), Uface(FIRST:LAST), SP_ID)
                  END DO
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == NOSLIP) then ! ++++++++ WALL NO-SLIP BOUNDARY ++++++++++++++++++++
                  Uface = U(:,I)
                  DO SP_ID = 1, N_SPECIES_FLUID
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq
                     Uface(FIRST+1) = 0.d0
                     Uface(FIRST+2) = 0.d0
                     Uface(FIRST+3) = GRID_BC(FACE_PG)%TEMP
                  END DO
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == MOVING) then ! ++++++++ MOVING BOUNDARY ++++++++++++++++++++
                  Uface = U(:,I)
                  DO SP_ID = 1, N_SPECIES_FLUID
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq
                     Uface(FIRST+1) = GRID_BC(FACE_PG)%UX
                     Uface(FIRST+2) = GRID_BC(FACE_PG)%UY
                     Uface(FIRST+3) = GRID_BC(FACE_PG)%TEMP
                  END DO
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == SYMMETRY) then ! ++++++++ SYM BOUNDARY ++++++++++++++++++++
                  ONAXIS = .TRUE.
                  Uface = U(:,I)
                  DO SP_ID = 1, N_SPECIES_FLUID
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq
                     Uface(FIRST+1) = U(FIRST+1,I) - (U(FIRST+1,I)*nx + U(FIRST+2,I)*ny)*nx
                     Uface(FIRST+2) = U(FIRST+2,I) - (U(FIRST+1,I)*nx + U(FIRST+2,I)*ny)*ny
                  END DO
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == KINETIC) then ! ++++++++ KINETIC BOUNDARY ++++++++++++++++++++
                  Uface = U(:,I)
                  DO SP_ID = 1, N_SPECIES_FLUID
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq
                     Uface(FIRST+1) = GRID_BC(FACE_PG)%UX
                     Uface(FIRST+2) = GRID_BC(FACE_PG)%UY
                     Uface(FIRST+3) = GRID_BC(FACE_PG)%TEMP
                  END DO
               else
                  print*, "ERROR! UNKNOWN BOUNDARY TYPE ", neigh, " for element ", I, &
                  " Check the mesh or the pre-processing."
                  print*, "ABORTING!"
                  STOP
               end if

            end if

            gradU(1,:,I) = gradU(1,:,I) + Uface*nx*Aface/Vcell
            gradU(2,:,I) = gradU(2,:,I) + Uface*ny*Aface/Vcell

         END DO

         !IF (ONAXIS) gradU(2,:,I) = 0.d0

      END DO

      DEALLOCATE(Uface)
      DEALLOCATE(U_adj)

   end subroutine


   subroutine compute_cell_centered_gradients_weighted_least_squares(U, gradU)

      ! This subroutine computes cell-centered gradients of the primitive variables using the 
      ! weighted-least-squares method of White [https://doi.org/10.2514/6.2019-0127]
      ! Coefficients for the least-squares solution are precomputed and stored in the
      ! matrix LSTSQ_COEFFS after the mesh is read.

      implicit none
   
      real(kind=8), dimension(:,:), intent(in)  :: U
      real(kind=8), dimension(:,:,:), intent(inout)  :: gradU

      INTEGER :: I, J, neigh, FACE_PG, SP_ID, FIRST, LAST, NEIGHBORPG
      REAL(KIND=8) :: nx, ny, Aface, Vcell
      REAL(KIND=8), DIMENSION(:), allocatable :: U_adj
      REAL(KIND=8) :: DX, DY, WJ

      LOGICAL :: ONAXIS, FLUIDBOUNDARY


      INTEGER :: VERT
      REAL(KIND=8) :: XVERT, YVERT, DIST, XC, YC
      REAL(KIND=8), DIMENSION(3) :: C1
      REAL(KIND=8) :: THETA, Twall

      ! Compute gradient in each cell from nodal values
      gradU = 0.d0
      ALLOCATE(U_adj(N_SPECIES_FLUID*Neq))

      DO I = 1, NCELLS
         ! Skip cells that are not fluid
         IF (U2D_GRID%CELL_PG(I) .NE. -1) THEN
            IF (GRID_BC(U2D_GRID%CELL_PG(I))%VOLUME_BC == SOLID) CYCLE
         END IF
         
         ONAXIS = .FALSE.
         DO J = 1, 3 ! The cell face
            ! Extract data
            Aface = U2D_GRID%CELL_FACES_AREA(J,I)
            nx = U2D_GRID%EDGE_NORMAL(1,J,I)
            ny = U2D_GRID%EDGE_NORMAL(2,J,I)
            Vcell = U2D_GRID%CELL_VOLUMES(I)

            FLUIDBOUNDARY = .FALSE.
            neigh = U2D_GRID%CELL_NEIGHBORS(J,I)
            IF (neigh == -1) THEN
               FLUIDBOUNDARY = .TRUE.
            ELSE
               NEIGHBORPG = U2D_GRID%CELL_PG(neigh)
               IF (NEIGHBORPG .NE. -1) THEN
                  IF (GRID_BC(NEIGHBORPG)%VOLUME_BC == SOLID) FLUIDBOUNDARY = .TRUE.
               END IF
            END IF

            IF (.NOT. FLUIDBOUNDARY) THEN
               U_adj = U(:,neigh)
            else
               FACE_PG = U2D_GRID%CELL_EDGES_PG(J,I)
               if (GRID_BC(FACE_PG)%PARTICLE_BC == STATE) then ! ++++++++ STATE BOUNDARY +++++++++++++++++++
                  DO SP_ID = 1, N_SPECIES_FLUID
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq
                     CALL compute_primitive_from_conserved(GRID_BC(FACE_PG)%U_BOUND(FIRST:LAST), U_adj(FIRST:LAST), SP_ID)
                  END DO
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == NOSLIP) then ! ++++++++ WALL NO-SLIP BOUNDARY ++++++++++++++++++++
                  U_adj = U(:,I)
                  DO SP_ID = 1, N_SPECIES_FLUID
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq
                     U_adj(FIRST+1) = -U(FIRST+1,I)
                     U_adj(FIRST+2) = -U(FIRST+2,I)
                     U_adj(FIRST+3) = 2.*GRID_BC(FACE_PG)%TEMP - U(FIRST+3,I)
                  END DO
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == MOVING) then ! ++++++++ MOVING BOUNDARY ++++++++++++++++++++
                  U_adj = U(:,I)
                  DO SP_ID = 1, N_SPECIES_FLUID
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq
                     U_adj(FIRST+1) = 2.*GRID_BC(FACE_PG)%UX - U(FIRST+1,I)
                     U_adj(FIRST+2) = 2.*GRID_BC(FACE_PG)%UY - U(FIRST+2,I)
                     U_adj(FIRST+3) = 2.*GRID_BC(FACE_PG)%TEMP - U(FIRST+3,I)
                  END DO
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == SYMMETRY) then ! ++++++++ SYM BOUNDARY ++++++++++++++++++++
                  ONAXIS = .TRUE.
                  U_adj = U(:,I)
                  DO SP_ID = 1, N_SPECIES_FLUID
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq
                     U_adj(FIRST+1) = U(FIRST+1,I) - 2.0*(U(FIRST+1,I)*nx + U(FIRST+2,I)*ny)*nx
                     U_adj(FIRST+2) = U(FIRST+2,I) - 2.0*(U(FIRST+1,I)*nx + U(FIRST+2,I)*ny)*ny
                  END DO
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == KINETIC) then ! ++++++++ KINETIC BOUNDARY ++++++++++++++++++++
                  U_adj = U(:,I)
                  
                  C1 = U2D_GRID%CELL_CENTROIDS(:,I)
                  THETA = ATAN2(C1(2), C1(1))
                  Twall = GRID_BC(FACE_PG)%TEMP + 200.d0*(SIN(THETA))**2
                  DO SP_ID = 1, N_SPECIES_FLUID
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq
                     U_adj(FIRST+1) = 2.*GRID_BC(FACE_PG)%UX - U(FIRST+1,I)
                     U_adj(FIRST+2) = 2.*GRID_BC(FACE_PG)%UY - U(FIRST+2,I)
                     U_adj(FIRST+3) = 2.*Twall - U(FIRST+3,I)
                  END DO
               else
                  print*, "ERROR! UNKNOWN BOUNDARY TYPE ", neigh, " for element ", I, &
                  " Check the mesh or the pre-processing."
                  print*, "ABORTING!"
                  STOP
               end if
            end if


            IF (FLUIDBOUNDARY) THEN
               XC = U2D_GRID%CELL_CENTROIDS(1, I)
               YC = U2D_GRID%CELL_CENTROIDS(2, I)

               VERT = U2D_GRID%CELL_NODES(J,I)
               XVERT = U2D_GRID%NODE_COORDS(1, VERT)
               YVERT = U2D_GRID%NODE_COORDS(2, VERT)

               DIST = -2.*((XC-XVERT)*nx + (YC-YVERT)*ny)
               WJ = 1./DIST
            ELSE
               DX = U2D_GRID%CELL_CENTROIDS(1, neigh) - U2D_GRID%CELL_CENTROIDS(1, I)
               DY = U2D_GRID%CELL_CENTROIDS(2, neigh) - U2D_GRID%CELL_CENTROIDS(2, I)
               WJ = 1./SQRT(DX*DX + DY*DY)
            END IF

            gradU(1,:,I) = gradU(1,:,I) + WJ*(U_adj - U(:,I))*U2D_GRID%LSTSQ_COEFFS(1,J,I)
            gradU(2,:,I) = gradU(2,:,I) + WJ*(U_adj - U(:,I))*U2D_GRID%LSTSQ_COEFFS(2,J,I)

         END DO

         !IF (ONAXIS) gradU(2,:,I) = 0.d0

      END DO

      DEALLOCATE(U_adj)

   end subroutine


end module
