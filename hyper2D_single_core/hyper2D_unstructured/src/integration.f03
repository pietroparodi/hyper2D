module integration

   use pde
   use grid
   use global_module
   use tools

   implicit none
   
   contains


   subroutine forward_Euler_step(U, U_new, dt)

      ! This function performs one step of the Forward Euler explicit time integrator 

      implicit none
   
      real(kind=8), dimension(:,:), intent(inout)  :: U, U_new

      real(kind=8), intent(in) :: dt
      integer :: eleID, eqID, intID
      
      real(kind=8), dimension(:,:), ALLOCATABLE :: gradUprim_neigh
      real(kind=8), dimension(:), allocatable :: F_dot_n_hyper, F_dot_n_diff, F_dot_n_wall, S
      real(kind=8) :: nx, ny, Aface, Acell, Vcell, dLR
      integer      :: neigh, FACE_PG, I, FIRST, LAST, NEIGHBORPG
      real(kind=8), dimension(3) :: C1, C2

      real(kind=8), dimension(:,:), ALLOCATABLE :: Uprim
      real(kind=8), dimension(:), ALLOCATABLE :: Uneigh
      real(kind=8), dimension(:,:,:), ALLOCATABLE :: gradUprim
      REAL(KIND=8) :: THETA, Twall

      LOGICAL :: FLUIDBOUNDARY

      ALLOCATE(F_dot_n_hyper(N_SPECIES_FLUID*Neq))
      ALLOCATE(F_dot_n_diff(N_SPECIES_FLUID*Neq))

      ALLOCATE(F_dot_n_wall(N_SPECIES_FLUID*Neq))
      ALLOCATE(S(N_SPECIES_FLUID*Neq))
      ALLOCATE(Uprim(N_SPECIES_FLUID*Neq,NCELLS))
      ALLOCATE(Uneigh(N_SPECIES_FLUID*Neq))
      ALLOCATE(gradUprim(2,N_SPECIES_FLUID*Neq,NCELLS))
      ALLOCATE(gradUprim_neigh(2,N_SPECIES_FLUID*Neq))

      DO eleID = 1, NCELLS
         ! Skip cells that are not fluid
         IF (U2D_GRID%CELL_PG(eleID) .NE. -1) THEN
            IF (GRID_BC(U2D_GRID%CELL_PG(eleID))%VOLUME_BC == SOLID) CYCLE
         END IF

         DO I = 1, N_SPECIES_FLUID
            FIRST = (I-1)*Neq+1
            LAST = I*Neq
            call compute_primitive_from_conserved(U(FIRST:LAST,eleID), Uprim(FIRST:LAST,eleID), I)
         END DO
      END DO

      CALL compute_cell_centered_gradients_weighted_least_squares(Uprim, gradUprim)
      !CALL compute_cell_centered_gradients_green_gauss(Uprim, gradUprim)

      do eleID = 1, NCELLS
         ! Skip cells that are not fluid
         IF (U2D_GRID%CELL_PG(eleID) .NE. -1) THEN
            IF (GRID_BC(U2D_GRID%CELL_PG(eleID))%VOLUME_BC == SOLID) CYCLE
         END IF

         Acell = U2D_GRID%CELL_AREAS(eleID)
         Vcell = U2D_GRID%CELL_VOLUMES(eleID)
         U_new(:,eleID) = U(:,eleID) ! Init

         do intID = 1, 3 ! Only triangular element supported

            F_dot_n_hyper = 0.d0
            F_dot_n_diff = 0.d0
            F_dot_n_wall = 0.d0

            ! Extract data
            Aface = U2D_GRID%CELL_FACES_AREA(intID,eleID)
            nx = U2D_GRID%EDGE_NORMAL(1,intID,eleID)
            ny = U2D_GRID%EDGE_NORMAL(2,intID,eleID)

            FLUIDBOUNDARY = .FALSE.
            neigh = U2D_GRID%CELL_NEIGHBORS(intID,eleID)
            IF (neigh == -1) THEN
               FLUIDBOUNDARY = .TRUE.
            ELSE
               NEIGHBORPG = U2D_GRID%CELL_PG(neigh)
               IF (NEIGHBORPG .NE. -1) THEN
                  IF (GRID_BC(NEIGHBORPG)%VOLUME_BC == SOLID) FLUIDBOUNDARY = .TRUE.
               END IF
            END IF




            IF (.NOT. FLUIDBOUNDARY) THEN ! +++++++++ INTERNAL CELL
               C1 = U2D_GRID%CELL_CENTROIDS(:,eleID)
               C2 = U2D_GRID%CELL_CENTROIDS(:,neigh)
               dLR = NORM2(C2-C1)
               Uneigh = U(:,neigh)
               gradUprim_neigh = gradUprim(:,:,neigh)
            else ! +++++++++ BOUNDARY CELL
               dLR = SQRT(Acell)
               FACE_PG = U2D_GRID%CELL_EDGES_PG(intID,eleID)

               if (GRID_BC(FACE_PG)%PARTICLE_BC == STATE) then ! ++++++++ STATE BOUNDARY +++++++++++++++++++
                  Uneigh = GRID_BC(FACE_PG)%U_BOUND
                  gradUprim_neigh = gradUprim(:,:,eleID)

                  !CALL compute_fluxes_diffusive(U(FIRST:LAST,eleID), U_neigh, &
                  !gradUprim(:,FIRST:LAST,eleID), gradUprim_neigh(:,FIRST:LAST), &
                  !nx, ny, F_dot_n_diff(FIRST:LAST), Acell, dLR, I, eleID)

               else if (GRID_BC(FACE_PG)%PARTICLE_BC == NOSLIP) then ! ++++++++ WALL BOUNDARY ++++++++++++++++++++
                  call compute_noslip_state(U(:,eleID), GRID_BC(FACE_PG)%TEMP, Uneigh)
                  gradUprim_neigh = gradUprim(:,:,eleID)
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == SYMMETRY) then ! ++++++++ SYMMETRY BOUNDARY ++++++++++++++++++++
                  call compute_sym_state(U(:,eleID), nx, ny, Uneigh)
                  gradUprim_neigh = gradUprim(:,:,eleID)
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == MOVING) then ! ++++++++ MOVING BOUNDARY ++++++++++++++++++++
                  call compute_moving_state(U(:,eleID), nx, ny, GRID_BC(FACE_PG)%TEMP, &
                  GRID_BC(FACE_PG)%UX, GRID_BC(FACE_PG)%UY, Uneigh)
                  gradUprim_neigh = gradUprim(:,:,eleID)
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == KINETIC) then ! ++++++++ KINETIC BOUNDARY ++++++++++++++++++++
                  C1 = U2D_GRID%CELL_CENTROIDS(:,eleID)
                  THETA = ATAN2(C1(2), C1(1))
                  Twall = GRID_BC(FACE_PG)%TEMP + 200.d0*(SIN(THETA))**2

                  call compute_moving_state(U(:,eleID), nx, ny, Twall, &
                  GRID_BC(FACE_PG)%UX, GRID_BC(FACE_PG)%UY, Uneigh)
                  !Uneigh = U(:,eleID)
                  gradUprim_neigh = gradUprim(:,:,eleID)
               else
                  print*, "ERROR! UNKNOWN BOUNDARY TYPE ", neigh, " for element ", eleID, &
                  " Check the mesh or the pre-processing."
                  print*, "ABORTING!"
                  STOP
               end if
            end if

            ! Compute the inviscid flux
            DO I = 1, N_SPECIES_FLUID
               FIRST = (I-1)*Neq+1
               LAST = I*Neq

               IF (FLUX_FUNCTION == AUSM) THEN
                  call compute_fluxes_AUSMplusup(U(FIRST:LAST,eleID), Uneigh(FIRST:LAST), &
                  nx, ny, F_dot_n_hyper(FIRST:LAST), Acell, I)
               ELSE IF (FLUX_FUNCTION == HLL) THEN
                  call compute_fluxes_HLL(U(FIRST:LAST,eleID), Uneigh(FIRST:LAST), &
                  nx, ny, F_dot_n_hyper(FIRST:LAST), Acell, I)
               ELSE IF (FLUX_FUNCTION == SLAU) THEN
                  call compute_fluxes_SLAU2(U(FIRST:LAST,eleID), Uneigh(FIRST:LAST), &
                  nx, ny, F_dot_n_hyper(FIRST:LAST), Acell, I)
               ELSE IF (FLUX_FUNCTION == CENTRAL) THEN
                  call compute_fluxes_central(U(FIRST:LAST,eleID), Uneigh(FIRST:LAST), &
                  nx, ny, F_dot_n_hyper(FIRST:LAST), Acell, I)
               ELSE
                  CALL ERROR_ABORT('Error! Flux function not available.')
               END IF

               ! call compute_fluxes_diffusive(U(FIRST:LAST,eleID), U_neigh, &
               ! gradUprim(:,FIRST:LAST,eleID), gradUprim_neigh, &
               ! nx, ny, F_dot_n_diff(FIRST:LAST), Acell, dLR, I, eleID)

            END DO

            CALL compute_fluxes_diffusive_binary(U(:,eleID), Uneigh, &
            gradUprim(:,:,eleID), gradUprim_neigh, &
            nx, ny, F_dot_n_diff(:), Acell, dLR, eleID, neigh)

            IF (FLUIDBOUNDARY) THEN
               FACE_PG = U2D_GRID%CELL_EDGES_PG(intID,eleID)
               IF (GRID_BC(FACE_PG)%PARTICLE_BC == KINETIC) THEN ! +++++++++ KINETIC FLUX AT THE WALL
                  F_dot_n_hyper = 0.d0
                  !F_dot_n_diff = 0.d0

                  CALL compute_kinetic_wall_fluxes(U(:,eleID), nx, ny, Twall, &
                  GRID_BC(FACE_PG)%UX, GRID_BC(FACE_PG)%UY, GRID_BC(FACE_PG)%REACT, F_dot_n_wall)
               END IF
            END IF
            
            !WRITE(*,*) F_dot_n_hyper
            !WRITE(*,*) F_dot_n_diff
            !WRITE(*,*) F_dot_n_wall


            ! Update solution
            U_new(:,eleID) = U_new(:,eleID) - dt*(F_dot_n_hyper + F_dot_n_diff + F_dot_n_wall)*Aface/Vcell

            !U_new(:,eleID) = U_new(:,eleID) - dt*(F_dot_n_hyper + F_dot_n_diff)*Aface/Vcell
            !WRITE(*,*) F_dot_n_diff
            !U_new(:,eleID) = U_new(:,eleID) - dt*(F_dot_n_hyper + F_dot_n_wall)*Aface/Vcell
            !U_new(:,eleID) = U_new(:,eleID) - dt*(F_dot_n_hyper)*Aface/Vcell



            ! Check that the solution did not diverge
            do eqID = 1, N_SPECIES_FLUID*Neq
               if (isnan(U_new(eqID,eleID))) then 
                  print*, 'Solution diverged, try with a smaller time step! Aborting.'
                  print*, 'Solution that diverged: ', U_new(:,eleID)
                  print*, 'in cell ID = ', eleID
                  ERROR STOP
               end if
            end do

         end do ! End loop on interfaces

         ! Add source terms
         call compute_source_terms(U(:,eleID), gradUprim(:,:,eleID), S, eleID)

         U_new(:,eleID) = U_new(:,eleID) + dt*S

   
      end do ! End loop on elements
   
      U = U_new ! Save results

      DEALLOCATE(F_dot_n_hyper)
      DEALLOCATE(F_dot_n_diff)
      DEALLOCATE(S)
      DEALLOCATE(Uprim)
      DEALLOCATE(Uneigh)
      DEALLOCATE(gradUprim)


   end subroutine


   subroutine compute_source_terms(U, gradUprim, S, eleID)

      ! This function computes the source terms of the equations
      ! These originate from reactive processes and intra-species collisions
      ! A term is also added when the equations are axisymmetric

      implicit none

      real(kind=8), dimension(:),   intent(in)  :: U
      real(kind=8), dimension(:,:), intent(in)  :: gradUprim
      real(kind=8), dimension(:),   intent(out) :: S
      INTEGER,   intent(in)                     :: eleID
      INTEGER :: I, J, JR, R1_SP_ID, R2_SP_ID, P, P_SP_ID
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: prim
      INTEGER :: FIRST, LAST, RHOI, RHOJ, MOMXI, MOMXJ, MOMYI, MOMYJ, ENEI, ENEJ, RHOP, MOMXP, MOMYP, ENEP
      REAL(KIND=8) :: MI, MJ, MIJ, NI, NJ, UXI, UXJ, UYI, UYJ, TI, TJ, TIJ, VTHIJ, S_MOM_IJ, UIDOTW, UJDOTW, S_ENE_IJ, QIJ
      REAL(KIND=8) :: K_FORWARD, TEMP, RATE_OF_PROGRESS
      REAL(KIND=8) :: MP, NP, UXP, UYP, EP, TP, EI
      REAL(KIND=8) :: QXI, QYI, QXJ, QYJ, S_SORET_IJ

      ALLOCATE(prim(N_SPECIES_FLUID*Neq))
      S = 0.d0
      
      DO I = 1, N_SPECIES_FLUID
         FIRST = (I-1)*Neq+1
         LAST = I*Neq
         call compute_primitive_from_conserved(U(FIRST:LAST), prim(FIRST:LAST), I)
      END DO

      ! Elastic collisions between the different fluids

      IF (.TRUE.) THEN
         DO I = 1, N_SPECIES_FLUID
            DO J = I+1, N_SPECIES_FLUID
               RHOI  = (I-1)*Neq+1
               RHOJ  = (J-1)*Neq+1
               MOMXI = (I-1)*Neq+2
               MOMXJ = (J-1)*Neq+2
               MOMYI = (I-1)*Neq+3
               MOMYJ = (J-1)*Neq+3
               ENEI  = (I-1)*Neq+4
               ENEJ  = (J-1)*Neq+4

               ! Momentum and energy elastic source terms from [Benilov, Phys. Plasmas 4, 521–528 (1997)]
               ! (in the low-Mach number limit)
               
               MI = SPECIES(I)%MOLECULAR_MASS
               MJ = SPECIES(J)%MOLECULAR_MASS
               MIJ = MI*MJ/(MI+MJ)

               NI = prim(RHOI)/MI
               NJ = prim(RHOJ)/MJ
               UXI = prim(MOMXI)
               UXJ = prim(MOMXJ)
               UYI = prim(MOMYI)
               UYJ = prim(MOMYJ)
               TI = prim(ENEI)
               TJ = prim(ENEJ)
               TIJ = (MI*TJ+MJ*TI)/(MI+MJ)

               VTHIJ = SQRT(8*KB*TIJ/(PI*MIJ))

               QIJ = PI*(0.5*(SPECIES(I)%DIAM + SPECIES(J)%DIAM))**2

               QXI = -KAPPA_GRID(I,I,eleID)*gradUprim(1,ENEI)
               QYI = -KAPPA_GRID(I,I,eleID)*gradUprim(2,ENEI)
               QXJ = -KAPPA_GRID(J,J,eleID)*gradUprim(1,ENEJ)
               QYJ = -KAPPA_GRID(J,J,eleID)*gradUprim(2,ENEJ)

               S_MOM_IJ = 4./3.*MIJ*VTHIJ*QIJ*NI*NJ
               S_SORET_IJ = 4./15.*MIJ*VTHIJ*QIJ/(TI*MJ+TJ*MI)/KB
               S(MOMXI) = S(MOMXI) + S_MOM_IJ*( UXJ - UXI ) + S_SORET_IJ*( QXJ*MI*NI - QXI*MJ*NJ )
               S(MOMXJ) = S(MOMXJ) - S_MOM_IJ*( UXJ - UXI ) - S_SORET_IJ*( QXJ*MI*NI - QXI*MJ*NJ )
               S(MOMYI) = S(MOMYI) + S_MOM_IJ*( UYJ - UYI ) + S_SORET_IJ*( QYJ*MI*NI - QYI*MJ*NJ )
               S(MOMYJ) = S(MOMYJ) - S_MOM_IJ*( UYJ - UYI ) - S_SORET_IJ*( QYJ*MI*NI - QYI*MJ*NJ )

               UIDOTW = UXI*(UXI-UXJ) + UYI*(UYI-UYJ)
               UJDOTW = UXJ*(UXI-UXJ) + UYJ*(UYI-UYJ)
               S_ENE_IJ = 4./3.*MIJ/(MI+MJ)*VTHIJ*QIJ*NI*NJ*(3*KB*(TI-TJ) + MI*TJ/TIJ*UIDOTW + MJ*TI/TIJ*UJDOTW)
               S(ENEI) = S(ENEI) - S_ENE_IJ
               S(ENEJ) = S(ENEJ) + S_ENE_IJ

               ! Update stability constraints
               invdt_coll = MAX(invdt_coll, 4./3.*NI*QIJ*VTHIJ, 4./3.*NJ*QIJ*VTHIJ)
            END DO
         END DO
      END IF

      ! Reactions with the background
      IF (BOOL_BG_FILE) THEN
         DO I = 1, N_SPECIES_FLUID ! The first reactant must be among the fluid species.

            RHOI  = (I-1)*Neq+1
            MOMXI = (I-1)*Neq+2
            MOMYI = (I-1)*Neq+3
            ENEI  = (I-1)*Neq+4

            MI = SPECIES(I)%MOLECULAR_MASS

            NI = prim(RHOI)/MI
            UXI = prim(MOMXI)
            UYI = prim(MOMYI)
            EI = U(ENEI)/prim(RHOI)
            TI = prim(ENEI)

            DO J = N_SPECIES_FLUID + 1, N_SPECIES ! The second reactant must be in the background.

               MJ = SPECIES(J)%MOLECULAR_MASS

               NJ = BG_CELL_NRHO(J,eleID)
               UXJ = BG_CELL_VX(J,eleID)
               UYJ = BG_CELL_VY(J,eleID)
               TJ = BG_CELL_TEMP(J,eleID)

               MIJ = MI*MJ/(MI+MJ)
               TIJ = (MI*TJ+MJ*TI)/(MI+MJ)
               VTHIJ = SQRT(8*KB*TIJ/(PI*MIJ))

               DO JR = 1, N_REACTIONS
                  R1_SP_ID = REACTIONS(JR)%R1_SP_ID
                  R2_SP_ID = REACTIONS(JR)%R2_SP_ID

                  IF ((R1_SP_ID .NE. I .OR. R2_SP_ID .NE. J) .AND. &
                     (R1_SP_ID .NE. J .OR. R2_SP_ID .NE. I)) CYCLE

                  TEMP = TIJ ! Temperature switch maybe necessary

                  IF (REACTIONS(JR)%TYPE == FIXED_RATE) THEN
                     K_FORWARD = REACTIONS(JR)%CONSTANT_RATE
                  ELSE IF (REACTIONS(JR)%TYPE == ARRHENIUS) THEN
                     K_FORWARD = REACTIONS(JR)%A * TEMP**REACTIONS(JR)%N * EXP(-REACTIONS(JR)%TA / TEMP)
                  ELSE IF (REACTIONS(JR)%TYPE == HARD_SPHERE) THEN
                     K_FORWARD = PI*REACTIONS(JR)%DIAM**2 * SQRT(8*KB*TEMP/(PI*MIJ)) * EXP(-REACTIONS(JR)%TA / TEMP) &
                     * (1. + REACTIONS(JR)%TA / TEMP)
                  ELSE IF (REACTIONS(JR)%TYPE == TABULATED) THEN
                     K_FORWARD = INTERP_RATE(TEMP*KB/QE, REACTIONS(JR)%TABLE_TEMP, REACTIONS(JR)%TABLE_RATE)
                  END IF

                  RATE_OF_PROGRESS = K_FORWARD*NI*NJ ! [1/m3/s]

                  ! Subtract reactants
                  S(RHOI) = S(RHOI) - MI*RATE_OF_PROGRESS
                  S(MOMXI) = S(MOMXI) - MI*RATE_OF_PROGRESS*UXI
                  S(MOMYI) = S(MOMYI) - MI*RATE_OF_PROGRESS*UYI
                  S(ENEI) = S(ENEI) - MI*RATE_OF_PROGRESS*EI + RATE_OF_PROGRESS*REACTIONS(JR)%DELTAE(I)
                  !S(ENEJ) = S(ENEJ) - 0.5*MJ*RATE_OF_PROGRESS*(UXJ**2 + UYJ**2) ! Check, this is not the total energy I think.

                  ! Add products (might be fluid or background)

                  DO P = 1, REACTIONS(JR)%N_PROD
                     IF (P == 1) P_SP_ID = REACTIONS(JR)%P1_SP_ID
                     IF (P == 2) P_SP_ID = REACTIONS(JR)%P2_SP_ID
                     IF (P == 3) P_SP_ID = REACTIONS(JR)%P3_SP_ID
                     IF (P == 4) P_SP_ID = REACTIONS(JR)%P4_SP_ID

                     IF (P_SP_ID .LE. N_SPECIES_FLUID) THEN ! Product is fluid
                        RHOP  = (P_SP_ID-1)*Neq+1
                        MOMXP = (P_SP_ID-1)*Neq+2
                        MOMYP = (P_SP_ID-1)*Neq+3
                        ENEP  = (P_SP_ID-1)*Neq+4

                        MP = SPECIES(P_SP_ID)%MOLECULAR_MASS

                        NP = prim(RHOP)/MP
                        UXP = prim(MOMXP)
                        UYP = prim(MOMYP)
                        EP = U(ENEP)/prim(RHOP)
                        TP = prim(ENEP)

                        S(RHOP) = S(RHOP) + MP*RATE_OF_PROGRESS
                        S(MOMXP) = S(MOMXP) + MP*RATE_OF_PROGRESS*UXP
                        S(MOMYP) = S(MOMYP) + MP*RATE_OF_PROGRESS*UYP
                        S(ENEP) = S(ENEP) + MP*RATE_OF_PROGRESS*EP + RATE_OF_PROGRESS*REACTIONS(JR)%DELTAE(P+2)

                     END IF

                  END DO


               END DO

            END DO
         END DO
      END IF


      IF (.FALSE.) THEN
         DO I = 1, N_SPECIES_FLUID
            RHOI  = (I-1)*Neq+1
            MOMYI = (I-1)*Neq+3
            ENEI  = (I-1)*Neq+4

            S(MOMYI) = S(MOMYI) - 9.81*prim(RHOI)
            S(ENEI) = S(ENEI) - 9.81*prim(RHOI)*prim(MOMYI)
         END DO
      END IF


      IF (AXI) THEN
         DO I = 1, N_SPECIES_FLUID
            RHOI  = (I-1)*Neq+1
            MOMYI = (I-1)*Neq+3
            ENEI  = (I-1)*Neq+4

            S(MOMYI) = S(MOMYI) + prim(RHOI) / SPECIES(I)%MOLECULAR_MASS*KB*prim(ENEI) &
            / U2D_GRID%CELL_CENTROIDS(2,eleID)
            !U2D_GRID%CELL_CENTROIDS(2,eleID)
            !WRITE(*,*)  U2D_GRID%CELL_VOLUMES(eleID) / (U2D_GRID%CELL_AREAS(eleID) * 2.*PI)
            ! U2D_GRID%CELL_CENTROIDS(2,eleID)
         END DO
      END IF

      DEALLOCATE(prim)

   end subroutine


   SUBROUTINE compute_kinetic_wall_fluxes(U, nx, ny, Twall, uwallx, uwally, react, F_dot_n_wall)

      ! This subroutine computes the flux at the interface with a solid
      ! surface. First, the exiting mass fluxes for eache species are computed,
      ! and then the entering fluxes due to surface catalysis are computed.

      implicit none

      real(kind=8), dimension(:), intent(in)  :: U
      real(kind=8),               intent(in)  :: nx, ny, Twall, uwallx, uwally
      LOGICAL,                    intent(in)  :: react
      real(kind=8), dimension(:), intent(out) :: F_dot_n_wall

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Ndotout, Ndotin
      REAL(KIND=8), DIMENSION(Neq) :: prim
      INTEGER :: JR, I, FIRST, LAST, P1_SP_ID, R_SP_ID
      REAL(KIND=8) :: M, T, rho, ux, uy, udotn, Pflux, SN, SToutX, SToutY, Pout, Pin
      REAL(KIND=8) :: uwalldotn, STinX, STinY, tauout, tauin

      ALLOCATE(Ndotout(N_SPECIES_FLUID))
      ALLOCATE(Ndotin(N_SPECIES_FLUID))

      Ndotin = 0.d0
      
      ! Compute exiting particle flux
      DO I = 1, N_SPECIES_FLUID
         FIRST = (I-1)*Neq+1
         LAST = I*Neq

         CALL compute_primitive_from_conserved(U(FIRST:LAST), prim, I)
         rho = prim(1)
         ux = prim(2)
         uy = prim(3)
         T = prim(4)
         M = SPECIES(I)%MOLECULAR_MASS
         SN = (ux*nx + uy*ny)*SQRT(M/(2.*KB*T))
         Ndotout(I) = rho/M*SQRT(KB*T/(2.*PI*M)) * (EXP(-SN**2) + SQRT(PI)*SN*(1+ERF(SN)))
      END DO

      ! Compute entering particle flux
      IF (react) THEN
         DO JR = 1, N_WALL_REACTIONS
            IF (WALL_REACTIONS(JR)%N_PROD == 1) THEN
               P1_SP_ID = WALL_REACTIONS(JR)%P1_SP_ID
               R_SP_ID = WALL_REACTIONS(JR)%R_SP_ID
               Ndotin(P1_SP_ID) = Ndotin(P1_SP_ID) + WALL_REACTIONS(JR)%PROB * Ndotout(R_SP_ID)
            END IF
         END DO
      ELSE
         Ndotin = Ndotout
      END IF

      ! Compute mass, momentum, and energy flux vector
      DO I = 1, N_SPECIES_FLUID
         FIRST = (I-1)*Neq+1
         LAST = I*Neq
         CALL compute_primitive_from_conserved(U(FIRST:LAST), prim, I)
         rho = prim(1)
         ux = prim(2)
         uy = prim(3)
         T = prim(4)
         M = SPECIES(I)%MOLECULAR_MASS
         udotn = ux*nx + uy*ny
         SN = udotn*SQRT(M/(2.*KB*T))
         SToutX = (ux - udotn*nx)*SQRT(M/(2.*KB*T))
         SToutY = (uy - udotn*ny)*SQRT(M/(2.*KB*T))

         uwalldotn = uwallx*nx + uwally*ny
         STinX = (uwallx - uwalldotn*nx)*SQRT(M/(2.*KB*Twall))
         STinY = (uwally - uwalldotn*ny)*SQRT(M/(2.*KB*Twall))

         ! Mass
         F_dot_n_wall(FIRST) = M * (Ndotout(I) - Ndotin(I))

         ! Momentum
         Pout = rho/M*KB*T * (SN/SQRT(PI)*EXP(-SN**2) + (0.5+SN**2)*(1+ERF(SN)))
         Pin = Ndotin(I)*0.5*SQRT(2.*PI*M*KB*Twall)
         Pflux = Pout + Pin
         tauout = Pout*(EXP(-SN**2)/SQRT(PI) + SN*(1+ERF(SN)))
         tauin = Pin/SQRT(PI)
         F_dot_n_wall(FIRST+1) = nx*Pflux + SToutX*tauout - STinX*tauin
         F_dot_n_wall(FIRST+2) = ny*Pflux + SToutY*tauout - STinY*tauin

         ! Energy
         F_dot_n_wall(FIRST+3) = Ndotout(I)*(0.5*M*(ux*ux+uy*uy) + 2.5*KB*T) &
         - 1./8.*rho/M * SQRT(8.*KB*T/(PI*M)) * KB*T*EXP(-SN**2) &
         - Ndotin(I)*(0.5*M*((uwallx - uwalldotn*nx)**2 + (uwally - uwalldotn*ny)**2) + 2.*KB*Twall)

      END DO

      !WRITE(*,*) F_dot_n_wall

   END SUBROUTINE





   SUBROUTINE compute_wall_fluxes_old(U, nx, ny, F_dot_n_wall)

      ! This subroutine computes the flux at the interface with a solid
      ! surface. First, the exiting mass fluxes for eache species are computed,
      ! and then the entering fluxes due to surface catalysis are computed.

      implicit none

      real(kind=8), dimension(:),   intent(in) :: U
      real(kind=8),                 intent(in) :: nx, ny
      real(kind=8), dimension(:), intent(out)  :: F_dot_n_wall

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Ndotout, Ndotin
      REAL(KIND=8), DIMENSION(Neq) :: prim
      INTEGER :: JR, I, FIRST, LAST, P1_SP_ID, R_SP_ID
      REAL(KIND=8) :: M, T, rho, ux, uy, udotn, uparsq, Pflux

      ALLOCATE(Ndotout(N_SPECIES_FLUID))
      ALLOCATE(Ndotin(N_SPECIES_FLUID))
      
      ! Compute exiting particle flux
      DO I = 1, N_SPECIES_FLUID
         FIRST = (I-1)*Neq+1
         LAST = I*Neq

         CALL compute_primitive_from_conserved(U(FIRST:LAST), prim, I)
         rho = prim(1)
         ux = prim(2)
         uy = prim(3)
         T = prim(4)
         M = SPECIES(I)%MOLECULAR_MASS
         Ndotout(I) = rho/M*SQRT(KB*T/(2.*PI*M))
      END DO

      ! Compute entering particle flux
      DO JR = 1, N_WALL_REACTIONS
         IF (WALL_REACTIONS(JR)%N_PROD == 1) THEN
            P1_SP_ID = WALL_REACTIONS(JR)%P1_SP_ID
            R_SP_ID = WALL_REACTIONS(JR)%R_SP_ID
            Ndotin(P1_SP_ID) = Ndotin(P1_SP_ID) + WALL_REACTIONS(JR)%PROB * Ndotout(R_SP_ID)
         END IF
      END DO

      ! Compute mass, momentum, and energy flux vector
      DO I = 1, N_SPECIES_FLUID
         FIRST = (I-1)*Neq+1
         LAST = I*Neq
         CALL compute_primitive_from_conserved(U(FIRST:LAST), prim, I)
         ux = prim(2)
         uy = prim(3)
         T = prim(4)
         M = SPECIES(I)%MOLECULAR_MASS

         ! Mass
         F_dot_n_wall(FIRST) = M * (Ndotout(I) - Ndotin(I))

         ! Momentum
         Pflux =  Ndotout(I)*0.5*SQRT(2.*PI*M*KB*T) + Ndotin(I)*0.5*SQRT(2.*PI*M*KB*Tw)
         F_dot_n_wall(FIRST+1) = nx*Pflux
         F_dot_n_wall(FIRST+2) = ny*Pflux

         ! Energy
         udotn = ux*nx + uy*ny
         uparsq = (ux - udotn*nx)**2 + (uy - udotn*ny)**2
         F_dot_n_wall(FIRST+3) = Ndotout(I)*(0.5*M*Uparsq + 2.*KB*T) - Ndotin(I)*(2.*KB*Tw)

      END DO

      !WRITE(*,*) F_dot_n_wall

   END SUBROUTINE


   subroutine compute_fluxes_HLL(U_L, U_R, nx, ny, F_dot_n, A_ele, SP_ID)

      ! Computes HLL numerical fluxes among the cell eleID and the neighbor cell neigh
      ! The element area Aele is also passed, for the sake of computing the CFL number.

      implicit none

      real(kind=8), dimension(:),   intent(in)  :: U_L, U_R
      real(kind=8),                 intent(in)  :: nx, ny, A_ele
      real(kind=8), dimension(Neq), intent(out) :: F_dot_n
      INTEGER,                      intent(in)  :: SP_ID

      real(kind=8), dimension(Neq) :: F_L, F_R

      ! Wave speeds
      real(kind=8) :: ws_min_L, ws_max_L, ws_min_R, ws_max_R, ws_min, ws_max

      call compute_flux_ws(U_L, F_L, nx, ny, ws_max_L, ws_min_L, SP_ID)
      call compute_flux_ws(U_R, F_R, nx, ny, ws_max_R, ws_min_R, SP_ID)

      ws_min = MIN(ws_min_L, ws_min_R)
      ws_max = MAX(ws_max_L, ws_max_R)

      ! HLL fluxes
      if (ws_min .ge. 0.0) then
         F_dot_n = F_L
      else if (ws_max .lt. 0.0) then
         F_dot_n = F_R
      else
         F_dot_n = (ws_min*ws_max*(U_R - U_L) + ws_max*F_L - ws_min*F_R)/(ws_max - ws_min)
      end if

      ! Update global maximum wave speed (used for setting the time step)
      ws_max = abs(ws_max)
      invdt_adv = MAX(invdt_adv, ws_max/sqrt(A_ele))

   end subroutine



   subroutine compute_fluxes_central(U_L, U_R, nx, ny, F_dot_n, A_ele, SP_ID)

      ! Computes HLL numerical fluxes among the cell eleID and the neighbor cell neigh
      ! The element area Aele is also passed, for the sake of computing the CFL number.

      implicit none

      real(kind=8), dimension(:),   intent(in)  :: U_L, U_R
      real(kind=8),                 intent(in)  :: nx, ny, A_ele
      real(kind=8), dimension(Neq), intent(out) :: F_dot_n
      INTEGER,                      intent(in)  :: SP_ID

      real(kind=8), dimension(Neq) :: F_L, F_R

      ! Wave speeds
      real(kind=8) :: ws_min_L, ws_max_L, ws_min_R, ws_max_R, ws_max

      call compute_flux_ws(U_L, F_L, nx, ny, ws_max_L, ws_min_L, SP_ID)
      call compute_flux_ws(U_R, F_R, nx, ny, ws_max_R, ws_min_R, SP_ID)

      F_dot_n = 0.5*(F_R + F_L)

      ! Update global maximum wave speed (used for setting the time step)
      ws_max  = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))
      invdt_adv = MAX(invdt_adv, ws_max/sqrt(A_ele))

   end subroutine

   subroutine compute_fluxes_AUSMplus(U_L, U_R, nx, ny, flux, A_ele, SP_ID)

      implicit none
    
      real(8), intent(in)  :: U_L(4), U_R(4)
      real(8), intent(in)  :: nx, ny
      real(8), intent(out) :: flux(4)

      INTEGER,                      intent(in)  :: SP_ID
      real(kind=8),                 intent(in)  :: A_ele
    
      real(8) :: rhoL, uL, vL, pL, HL, aL
      real(8) :: rhoR, uR, vR, pR, HR, aR
      real(8) :: unL, unR
      real(8) :: ML, MR
      real(8) :: Mbar, abar
      real(8) :: Mp, Mm
      real(8) :: Pp, Pm
      real(8) :: mdot
      real(8) :: pflux
      real(8) :: beta, kp, ku, sigma
      real(8) :: EL, ER
      real(8) :: gamma
    
      beta  = 0.125d0
      kp    = 0.25d0
      ku    = 0.75d0
      sigma = 1.0d0
    

      gamma = SPECIES(SP_ID)%GAMMA
      !------------------------
      ! Left state
      !------------------------
    
      rhoL = U_L(1)
      uL   = U_L(2)/rhoL
      vL   = U_L(3)/rhoL
      EL   = U_L(4)
    
      pL = (gamma-1d0)*(EL - 0.5d0*rhoL*(uL*uL + vL*vL))
    
      HL = (EL + pL)/rhoL
      aL = sqrt(gamma*pL/rhoL)
    
      !------------------------
      ! Right state
      !------------------------
    
      rhoR = U_R(1)
      uR   = U_R(2)/rhoR
      vR   = U_R(3)/rhoR
      ER   = U_R(4)
    
      pR = (gamma-1d0)*(ER - 0.5d0*rhoR*(uR*uR + vR*vR))
    
      HR = (ER + pR)/rhoR
      aR = sqrt(gamma*pR/rhoR)
    
      !------------------------
      ! Normal velocities
      !------------------------
    
      unL = uL*nx + vL*ny
      unR = uR*nx + vR*ny
    
      !------------------------
      ! Reference speed of sound
      !------------------------
    
      abar = 0.5d0*(aL + aR)
    
      ! Mach numbers
      ML = unL / abar
      MR = unR / abar
    
      !------------------------
      ! Mach splitting (AUSM+)
      !------------------------
    
      if (abs(ML) >= 1d0) then
          Mp = 0.5d0*(ML + abs(ML))
      else
          Mp = 0.25d0*(ML + 1d0)**2 + beta*(ML**2 - 1d0)**2
      end if
    
      if (abs(MR) >= 1d0) then
          Mm = 0.5d0*(MR - abs(MR))
      else
          Mm = -0.25d0*(MR - 1d0)**2 - beta*(MR**2 - 1d0)**2
      end if
    
      Mbar = Mp + Mm
    
      !------------------------
      ! Pressure splitting
      !------------------------
    
      if (abs(ML) >= 1d0) then
          Pp = 0.5d0*(1d0 + sign(1d0,ML))
      else
          Pp = 0.25d0*(ML + 1d0)**2*(2d0 - ML) &
               + kp*(ML**2 - 1d0)**2
      end if
    
      if (abs(MR) >= 1d0) then
          Pm = 0.5d0*(1d0 - sign(1d0,MR))
      else
          Pm = 0.25d0*(MR - 1d0)**2*(2d0 + MR) &
               - kp*(MR**2 - 1d0)**2
      end if
    
      pflux = Pp*pL + Pm*pR
    
      !------------------------
      ! Mass flux
      !------------------------
    
      mdot = abar * Mbar
    
      !------------------------
      ! Upwind selection
      !------------------------
    
      if (mdot >= 0d0) then
    
          flux(1) = mdot * rhoL
          flux(2) = mdot * rhoL * uL + pflux*nx
          flux(3) = mdot * rhoL * vL + pflux*ny
          flux(4) = mdot * rhoL * HL
    
      else
    
          flux(1) = mdot * rhoR
          flux(2) = mdot * rhoR * uR + pflux*nx
          flux(3) = mdot * rhoR * vR + pflux*ny
          flux(4) = mdot * rhoR * HR
    
      end if

      invdt_adv = MAX(invdt_adv, 2.*aL/sqrt(A_ele), 2.*aR/sqrt(A_ele))

    
   end subroutine


   subroutine compute_fluxes_AUSMplusup(U_L, U_R, nx, ny, flux, A_ele, SP_ID)

      ! Computes fluxes at the interface between states U_L and U_R using the
      ! AUSM+-up for all speeds method of Liou [https://doi.org/10.1016/j.jcp.2005.09.020]
      ! The reference Mach number is fixed at 0.3

      implicit none
      
      real(8), intent(in)  :: U_L(4), U_R(4)
      real(8), intent(in)  :: nx, ny
      real(8), intent(out) :: flux(4)

      INTEGER,                      intent(in)  :: SP_ID
      real(kind=8),                 intent(in)  :: A_ele
      
      real(8) :: rhoL,uL,vL,pL,HL,aL
      real(8) :: rhoR,uR,vR,pR,HR,aR
      real(8) :: unL,unR
      real(8) :: ML,MR
      real(8) :: aface
      real(8) :: mdot
      real(8) :: rhoface
      real(8) :: beta,Kp,Ku,sigma
      real(8) :: EL,ER
      real(8) :: gamma
      real(8) :: astarL, astarR, ahatL, ahatR, Mfsq, Mrefsq, fa, alpha, p1, p2
      real(8) :: MLP, MRM, betaLP, betaRM, pu, Mface, pflux, Minf, MP
      
      beta = 0.125d0
      Kp   = 0.25d0
      Ku   = 0.75d0
      sigma = 1.0d0
      
      gamma = SPECIES(SP_ID)%GAMMA

      Minf = 0.3

      !------------------------
      ! Left state
      !------------------------
      
      rhoL = U_L(1)
      uL   = U_L(2)/rhoL
      vL   = U_L(3)/rhoL
      EL   = U_L(4)
      
      pL = (gamma-1.d0)*(EL - 0.5d0*rhoL*(uL*uL + vL*vL))
      HL = (EL + pL)/rhoL
      aL = sqrt(gamma*pL/rhoL)
      
      !------------------------
      ! Right state
      !------------------------
      
      rhoR = U_R(1)
      uR   = U_R(2)/rhoR
      vR   = U_R(3)/rhoR
      ER   = U_R(4)
      
      pR = (gamma-1.d0)*(ER - 0.5d0*rhoR*(uR*uR + vR*vR))
      HR = (ER + pR)/rhoR
      aR = sqrt(gamma*pR/rhoR)
      
      !------------------------
      ! Normal velocity
      !------------------------
      
      unL = uL*nx + vL*ny
      unR = uR*nx + vR*ny
      
      !------------------------
      ! Interface sound speed
      !------------------------
      
      astarL = SQRT(2.*(gamma-1)/(gamma+1)*HL)
      astarR = SQRT(2.*(gamma-1)/(gamma+1)*HR)

      ahatL = astarL*astarL/MAX(astarL,  unL)
      ahatR = astarR*astarR/MAX(astarR, -unR)

      aface = MIN(ahatL,ahatR)
            
      ML = unL / aface
      MR = unR / aface

      Mfsq = 0.5*(ML*ML+MR*MR)
      Mrefsq = MIN(1., MAX(Mfsq, Minf*Minf))

      fa = 2.*SQRT(Mrefsq) - Mrefsq

      alpha = 3.0/16.0*(-4.0+5.0*fa*fa)

      IF (ABS(ML) <= 1.0) THEN
         p1 = 0.25*(ML+1.0)**2
         p2 = (ML*ML-1.0)**2

         MLP = p1 + beta*p2
         betaLP = p1*(2.0-ML) + alpha*ML*p2
      ELSE
         MLP = 0.5*(ML+ABS(ML))
         betaLP = MLP/ML
      END IF

      IF (ABS(MR) <= 1.0) THEN
         p1 = 0.25*(mR-1.0)*(mR-1.0)
         p2 = (mR*mR-1.0)*(mR*mR-1.0)

         MRM = -p1 - beta*p2
         betaRM = p1*(2.0+MR) - alpha*MR*p2
      ELSE
         MRM = 0.5*(MR-ABS(MR))
         betaRM = MRM/MR
      END IF




      rhoface = 0.5*(rhoL + rhoR)
      MP = -(Kp/fa)*MAX((1.0-sigma*Mfsq),0.0)*(pR-pL)/(rhoface*aface*aface)

      Pu = -Ku*fa*betaLP*betaRM*2.0*rhoface*aface*(unR-unL)


      Mface = MLP + MRM !+ MP

      pflux = betaLP*pL + betaRM*pR !+ Pu


      !------------------------
      ! Flux
      !------------------------

      if (Mface >= 0.d0) then
         mdot = Mface * aface * rhoL
         flux(1) = mdot
         flux(2) = mdot*uL + pflux*nx
         flux(3) = mdot*vL + pflux*ny
         flux(4) = mdot*HL
      else
         mdot = Mface * aface * rhoR
         flux(1) = mdot
         flux(2) = mdot*uR + pflux*nx
         flux(3) = mdot*vR + pflux*ny
         flux(4) = mdot*HR
      endif

      invdt_adv = MAX(invdt_adv, 2.*aL/sqrt(A_ele), 2.*aR/sqrt(A_ele))

      
   end subroutine




   subroutine compute_fluxes_SLAU2(U_L, U_R, nx, ny, flux, A_ele, SP_ID)

      ! Computes fluxes at the interface between states U_L and U_R using the
      ! SLAU2 method of Kitamura and Shima [https://doi.org/10.1016/j.jcp.2013.02.046]

      implicit none
      
      real(8), intent(in)  :: U_L(4), U_R(4)
      real(8), intent(in)  :: nx, ny
      real(8), intent(out) :: flux(4)

      INTEGER,                      intent(in)  :: SP_ID
      real(kind=8),                 intent(in)  :: A_ele
      
      real(8) :: rhoL,uL,vL,pL,HL,aL
      real(8) :: rhoR,uR,vR,pR,HR,aR
      real(8) :: unL,unR
      real(8) :: ML,MR
      real(8) :: mdot

      real(8) :: EL,ER
      real(8) :: gamma
      real(8) :: abar, unbar, Mhat, chi, Mbar, pbar, rhobar, betaL, betaR, deltap, deltarho, deltabeta
      real(8) :: unbarL, unbarR, g, pflux
      

      gamma = SPECIES(SP_ID)%GAMMA

      !------------------------
      ! Left state
      !------------------------
      
      rhoL = U_L(1)
      uL   = U_L(2)/rhoL
      vL   = U_L(3)/rhoL
      EL   = U_L(4)
      
      pL = (gamma-1.d0)*(EL - 0.5d0*rhoL*(uL*uL + vL*vL))
      HL = (EL + pL)/rhoL
      aL = sqrt(gamma*pL/rhoL)
      
      !------------------------
      ! Right state
      !------------------------
      
      rhoR = U_R(1)
      uR   = U_R(2)/rhoR
      vR   = U_R(3)/rhoR
      ER   = U_R(4)
      
      pR = (gamma-1.d0)*(ER - 0.5d0*rhoR*(uR*uR + vR*vR))
      HR = (ER + pR)/rhoR
      aR = sqrt(gamma*pR/rhoR)
      
      !------------------------
      ! Normal velocity
      !------------------------
      
      unL = uL*nx + vL*ny
      unR = uR*nx + vR*ny
      
      !------------------------
      ! Interface sound speed
      !------------------------

      abar = 0.5*(aL + aR)
      
      unbar = (rhoL*ABS(unL) + rhoR*ABS(unR))/(rhoL+rhoR)

      Mhat = MIN(1.0, SQRT(0.5*(uL*uL + vL*vL + uR*uR + vR*vR)) / abar)

      chi = (1.0 - Mhat)**2

      Mbar = unbar/abar
      pbar = 0.5*(pL+pR)
      rhobar = 0.5*(rhoL+rhoR)

      ML = unL/abar
      MR = unR/abar

      IF (ABS(ML) < 1.0) THEN
         betaL = 0.25*(2.0-ML)*(ML+1.0)**2
      ELSE
         betaL = 0.5*(1.0+SIGN(1.d0,ML))
      ENDIF

      IF (ABS(MR) < 1.0) THEN
         betaR = 0.25*(2.0+MR)*(MR-1.0)**2
      ELSE
         betaR = 0.5*(1.0+SIGN(1.d0,-MR))
      ENDIF

      deltap = pR-pL
      deltarho = rhoR-rhoL
      deltabeta = betaL-betaR

      g = -MAX(MIN(ML, 0.0), -1.0) * MIN(MAX(MR, 0.0), 1.0)
      
      !------------------------
      ! Flux
      !------------------------

      pflux = pbar + (0.5*deltabeta*(pL-pR)) + SQRT(0.5*(uL*uL + vL*vL + uR*uR + vR*vR)) * &
                     (betaL + betaR - 1.0) * rhobar * abar
      
      unbarL = (1.-g)*unbar + g*ABS(unL)
      unbarR = (1.-g)*unbar + g*ABS(unR)

      mdot = 0.5*(rhoL*(unL+unbarL) + rhoR*(unR-unbarR) - chi/abar*deltap)

      if (mdot >= 0.d0) then
         flux(1) = mdot
         flux(2) = mdot*uL + pflux*nx
         flux(3) = mdot*vL + pflux*ny
         flux(4) = mdot*HL
      else
         flux(1) = mdot
         flux(2) = mdot*uR + pflux*nx
         flux(3) = mdot*vR + pflux*ny
         flux(4) = mdot*HR
      endif

      invdt_adv = MAX(invdt_adv, 2.*aL/sqrt(A_ele), 2.*aR/sqrt(A_ele))

      
   end subroutine


   subroutine compute_transport_simple(U)

      implicit none

      real(kind=8), dimension(:,:), intent(in)  :: U
      INTEGER :: eleID, IS, JS

      DO eleID = 1, NCELLS

         DO IS = 1, N_SPECIES_FLUID
            MU_GRID(IS, eleID) = SPECIES(IS)%MU
            DO JS = 1, N_SPECIES_FLUID
               KAPPA_GRID(IS, JS, eleID) = SPECIES(IS)%KAPPA
            END DO
         END DO

      END DO

   end subroutine



   subroutine compute_transport(U)

      implicit none

      real(kind=8), dimension(:,:), intent(in)  :: U
      real(kind=8), dimension(Neq) :: primI2, primI
      REAL(KIND=8) :: T_I, T_I2, X_I, X_I2, n_I, n_I2
      INTEGER :: eleID

      DO eleID = 1, NCELLS

         ! Skip cells that are not fluid
         IF (U2D_GRID%CELL_PG(eleID) .NE. -1) THEN
            IF (GRID_BC(U2D_GRID%CELL_PG(eleID))%VOLUME_BC == SOLID) CYCLE
         END IF


         call compute_primitive_from_conserved(U(1:4, eleID), primI2, 1)
         call compute_primitive_from_conserved(U(5:8, eleID), primI, 2)

         T_I2 = primI2(4)
         T_I = primI(4)

         n_I2 = primI2(1)/SPECIES(1)%MOLECULAR_MASS
         n_I = primI(1)/SPECIES(2)%MOLECULAR_MASS
         
         X_I2 = n_I2/(n_I2 + n_I)
         X_I = n_I/(n_I2 + n_I)

         ! MU_GRID(1, eleID) = SPECIES(1)%MU
         ! MU_GRID(2, eleID) = SPECIES(2)%MU
         ! KAPPA_GRID(1,1, eleID) = KAPPA_I2_I2(T_I2, 0.d0)
         ! KAPPA_GRID(1,2, eleID) = 0.d0
         ! KAPPA_GRID(2,1, eleID) = 0.d0
         ! KAPPA_GRID(2,2, eleID) = SPECIES(2)%KAPPA

         ! MU_GRID(1, eleID) = MU_I2(T_I2, 0.d0)
         ! MU_GRID(2, eleID) = MU_I(T_I, 1.d0)
         ! KAPPA_GRID(1,1, eleID) = KAPPA_I2_I2(T_I2, 0.d0)
         ! KAPPA_GRID(1,2, eleID) = 0.d0
         ! KAPPA_GRID(2,1, eleID) = 0.d0
         ! KAPPA_GRID(2,2, eleID) = KAPPA_I_I(T_I, 1.d0)

         MU_GRID(1, eleID) = MU_I2(T_I2, X_I)
         MU_GRID(2, eleID) = MU_I(T_I, X_I)
         KAPPA_GRID(1,1, eleID) = KAPPA_I2_I2(T_I2, X_I)
         KAPPA_GRID(1,2, eleID) = KAPPA_I2_I(T_I2, X_I)
         KAPPA_GRID(2,1, eleID) = KAPPA_I_I2(T_I, X_I)
         KAPPA_GRID(2,2, eleID) = KAPPA_I_I(T_I, X_I)

      END DO

   end subroutine


   subroutine compute_fluxes_diffusive(U_L, U_R, gradU_L, gradU_R, nx, ny, F_dot_n, A_ele, dLR, SP_ID, eleID)

      ! Computes HLL numerical fluxes among the cell eleID and the neighbor cell neigh
      ! The cell state is U_L, and the reighboring cell state is U_R.
      ! The corresponding cell-centered gradients for the primitive variables 
      ! gradU_L and gradU_R are also required.
      ! The element area Aele is also passed, for the sake of computing the CFL number.

      implicit none

      real(kind=8), dimension(:),   intent(in)  :: U_L, U_R
      real(kind=8), dimension(:,:), intent(in)  :: gradU_L, gradU_R
      real(kind=8),                 intent(in)  :: nx, ny, A_ele, dLR
      real(kind=8), dimension(Neq), intent(out) :: F_dot_n
      INTEGER,                      intent(in)  :: SP_ID, eleID

      real(kind=8), dimension(Neq) :: prim

      real(kind=8) :: ux_L, ux_R, uy_L, uy_R, T_L, T_R, rho_L
      real(kind=8) :: DUXDX, DUYDX, DUXDY, DUYDY, TAUXX, TAUXY, TAUYY, DTDX, DTDY, UX, UY
      REAL(KIND=8) :: MU, KAPPA, CP
      REAL(KIND=8) :: MU_LIM, KAPPA_LIM


      F_dot_n = 0.d0

      rho_L = U_L(1)

      !MU    = SPECIES(SP_ID)%MU
      !KAPPA = SPECIES(SP_ID)%KAPPA
      MU = MU_GRID(SP_ID, eleID)
      KAPPA = KAPPA_GRID(SP_ID,SP_ID, eleID)
      CP    = SPECIES(SP_ID)%CP

      MU_LIM = CFL_target*(A_ele*rho_L)/(6.*dt_target)
      KAPPA_LIM = CFL_target*(A_ele*rho_L*CP)/(6.*dt_target)

      !MU = MIN(MU, MU_LIM)
      !KAPPA = MIN(KAPPA, KAPPA_LIM)

      DUXDX = 0.5*(gradU_L(1,2) + gradU_R(1,2))
      DUYDX = 0.5*(gradU_L(1,3) + gradU_R(1,3))
      DUXDY = 0.5*(gradU_L(2,2) + gradU_R(2,2))
      DUYDY = 0.5*(gradU_L(2,3) + gradU_R(2,3))

      
      TAUXX = 2./3.*MU*(2.*DUXDX-DUYDY)
      TAUYY = 2./3.*MU*(2.*DUYDY-DUXDX)
      TAUXY = MU*(DUXDY + DUYDX)

      F_dot_n(2) = - (NX*TAUXX + NY*TAUXY)
      F_dot_n(3) = - (NX*TAUXY + NY*TAUYY)



      DTDX = 0.5*(gradU_L(1,4) + gradU_R(1,4))
      DTDY = 0.5*(gradU_L(2,4) + gradU_R(2,4))

      ux_L = U_L(2)/U_L(1)
      ux_R = U_R(2)/U_R(1)

      uy_L = U_L(3)/U_L(1)
      uy_R = U_R(3)/U_R(1)

      UX = 0.5*(ux_L + ux_R)
      UY = 0.5*(uy_L + uy_R)

      call compute_primitive_from_conserved(U_L, prim, SP_ID)
      T_L = prim(4)
      rho_L = prim(1)
      call compute_primitive_from_conserved(U_R, prim, SP_ID)
      T_R = prim(4)

      F_dot_n(4) = - KAPPA*(T_R-T_L)/dLR

      !F_dot_n(4) = - (NX*DTDX + NY*DTDY)*KAPPA &
      !             - (UX*TAUXX + UY*TAUXY)*NX &
      !             - (UX*TAUXY + UY*TAUYY)*NY

      F_dot_n(4) = F_dot_n(4) - (UX*TAUXX + UY*TAUXY)*NX &
                              - (UX*TAUXY + UY*TAUYY)*NY

      invdt_cond = MAX(invdt_cond, 6.*KAPPA/(A_ele*rho_L*CP))
      invdt_diff = MAX(invdt_diff, 6.*MU/(A_ele*rho_L))



   end subroutine



   subroutine compute_fluxes_diffusive_binary(U_L, U_R, gradU_L, gradU_R, nx, ny, F_dot_n, A_ele, dLR, eleID, neigh)

      implicit none

      real(kind=8), dimension(:),   intent(in)  :: U_L, U_R
      real(kind=8), dimension(:,:), intent(in)  :: gradU_L, gradU_R
      real(kind=8),                 intent(in)  :: nx, ny, A_ele, dLR
      real(kind=8), dimension(:),   intent(out) :: F_dot_n
      INTEGER,                      intent(in)  :: eleID, neigh

      real(kind=8), dimension(Neq) :: prim

      real(kind=8) :: ux_L, ux_R, uy_L, uy_R, T_L, T_R, rho_L
      real(kind=8) :: DUXDX, DUYDX, DUXDY, DUYDY, TAUXX, TAUXY, TAUYY, DTDX, DTDY, UX, UY
      REAL(KIND=8) :: MU, KAPPA, CP, KAPPA_L, KAPPA_R, MU_L, MU_R

      INTEGER :: I, J, FIRST, LAST, FIRSTJ, LASTJ

      F_dot_n = 0.d0

      DO I = 1, N_SPECIES_FLUID ! The species to which the fluxes are applied
         FIRST = (I-1)*Neq+1
         LAST = I*Neq

         !MU    = SPECIES(SP_ID)%MU
         !KAPPA = SPECIES(SP_ID)%KAPPA
         !MU = 0.5*(MU_GRID(I, eleID) + MU_GRID(I, neigh))
         MU_L = MU_GRID(I, eleID)
         IF (neigh == -1) THEN
            MU_R = MU_L
         ELSE
            MU_R = MU_GRID(I, neigh)
         END IF
         !MU = 2.*MU_L*MU_R/(MU_L+MU_R)
         MU = 0.5*(MU_L + MU_R)


         !MU_LIM = CFL_target*(A_ele*rho_L)/(6.*dt_target)
         !KAPPA_LIM = CFL_target*(A_ele*rho_L*CP)/(6.*dt_target)

         !MU = MIN(MU, MU_LIM)
         !KAPPA = MIN(KAPPA, KAPPA_LIM)

         DUXDX = 0.5*(gradU_L(1,FIRST+1) + gradU_R(1,FIRST+1))
         DUYDX = 0.5*(gradU_L(1,FIRST+2) + gradU_R(1,FIRST+2))
         DUXDY = 0.5*(gradU_L(2,FIRST+1) + gradU_R(2,FIRST+1))
         DUYDY = 0.5*(gradU_L(2,FIRST+2) + gradU_R(2,FIRST+2))

         
         TAUXX = 2./3.*MU*(2.*DUXDX-DUYDY)
         TAUYY = 2./3.*MU*(2.*DUYDY-DUXDX)
         TAUXY = MU*(DUXDY + DUYDX)

         F_dot_n(FIRST+1) = - (NX*TAUXX + NY*TAUXY)
         F_dot_n(FIRST+2) = - (NX*TAUXY + NY*TAUYY)



         DTDX = 0.5*(gradU_L(1,FIRST+3) + gradU_R(1,FIRST+3))
         DTDY = 0.5*(gradU_L(2,FIRST+3) + gradU_R(2,FIRST+3))

         ux_L = U_L(2)/U_L(1)
         ux_R = U_R(2)/U_R(1)

         uy_L = U_L(3)/U_L(1)
         uy_R = U_R(3)/U_R(1)

         UX = 0.5*(ux_L + ux_R)
         UY = 0.5*(uy_L + uy_R)


         rho_L = U_L(FIRST)
  
         F_dot_n(FIRST+3) = F_dot_n(FIRST+3) - (UX*TAUXX + UY*TAUXY)*NX &
                                             - (UX*TAUXY + UY*TAUYY)*NY


         invdt_diff = MAX(invdt_diff, 6.*MU/(A_ele*rho_L))

         DO J = 1, N_SPECIES_FLUID ! The species that appears in the gradient
            FIRSTJ = (J-1)*Neq+1
            LASTJ = J*Neq

            KAPPA_L = KAPPA_GRID(I,J, eleID)
            IF (neigh == -1) THEN
               KAPPA_R = KAPPA_L
            ELSE
               KAPPA_R = KAPPA_GRID(I,J, neigh)
            END IF
            KAPPA = 0.5*(KAPPA_L + KAPPA_R)
            CP    = SPECIES(I)%CP

            rho_L = U_L(FIRST)
            call compute_primitive_from_conserved(U_L(FIRSTJ:LASTJ), prim, J)
            T_L = prim(4)
            call compute_primitive_from_conserved(U_R(FIRSTJ:LASTJ), prim, J)
            T_R = prim(4)

            F_dot_n(FIRST+3) = F_dot_n(FIRST+3) - KAPPA*(T_R-T_L)/dLR !KAPPA*(NX*DTDX + NY*DTDY) ! 

            invdt_cond = MAX(invdt_cond, 6.*KAPPA/(A_ele*rho_L*CP))
         END DO
      END DO

   end subroutine


end module
