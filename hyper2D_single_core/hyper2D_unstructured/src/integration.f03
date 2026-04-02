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
      
      real(kind=8), dimension(Neq) :: U_sym, U_wall, U_neigh
      real(kind=8), dimension(2,Neq) :: gradUprim_neigh
      real(kind=8), dimension(:), allocatable :: F_dot_n_hyper, F_dot_n_diff, F_dot_n_wall, S
      real(kind=8) :: nx, ny, Aface, Acell, Vcell, dLR
      integer      :: neigh, FACE_PG, I, FIRST, LAST, NEIGHBORPG
      real(kind=8), dimension(3) :: C1, C2

      real(kind=8), dimension(:,:), ALLOCATABLE :: Uprim
      real(kind=8), dimension(:,:,:), ALLOCATABLE :: gradUprim

      LOGICAL :: FLUIDBOUNDARY

      ALLOCATE(F_dot_n_hyper(N_SPECIES_FLUID*Neq))
      ALLOCATE(F_dot_n_diff(N_SPECIES_FLUID*Neq))

      ALLOCATE(F_dot_n_wall(N_SPECIES_FLUID*Neq))
      ALLOCATE(S(N_SPECIES_FLUID*Neq))
      ALLOCATE(Uprim(N_SPECIES_FLUID*Neq,NCELLS))
      ALLOCATE(gradUprim(2,N_SPECIES_FLUID*Neq,NCELLS))

      DO eleID = 1, NCELLS
         ! Skip cells that are not fluid
         IF (U2D_GRID%CELL_PG(eleID) .NE. -1) THEN
            IF (GRID_BC(U2D_GRID%CELL_PG(eleID))%VOLUME_BC == SOLID) CYCLE
         END IF

         DO I = 1, N_SPECIES_FLUID
            FIRST = (I-1)*Neq+1
            LAST = I*Neq+1
            call compute_primitive_from_conserved(U(FIRST:LAST,eleID), Uprim(FIRST:LAST,eleID), I)
         END DO
      END DO

      call compute_cell_centered_gradients_weighted_least_squares(Uprim, gradUprim)

      do eleID = 1, NCELLS
         ! Skip cells that are not fluid
         IF (U2D_GRID%CELL_PG(eleID) .NE. -1) THEN
            IF (GRID_BC(U2D_GRID%CELL_PG(eleID))%VOLUME_BC == SOLID) CYCLE
         END IF

         Acell = U2D_GRID%CELL_AREAS(eleID)
         Vcell = U2D_GRID%CELL_VOLUMES(eleID)
         U_new(:,eleID) = U(:,eleID) ! Init

         do intID = 1, 3 ! Only quad element supported

            F_dot_n_hyper = 0.0 ! Init
            F_dot_n_diff = 0.0 ! Init

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


            ! Check what neighbor is it
            DO I = 1, N_SPECIES_FLUID
               FIRST = (I-1)*Neq+1
               LAST = I*Neq+1

               IF (.NOT. FLUIDBOUNDARY) THEN ! +++++++++ INTERNAL CELL

                  C1 = U2D_GRID%CELL_CENTROIDS(:,eleID)
                  C2 = U2D_GRID%CELL_CENTROIDS(:,neigh)
                  dLR = NORM2(C2-C1)
                  U_neigh = U(FIRST:LAST,neigh)
                  gradUprim_neigh = gradUprim(:,FIRST:LAST,neigh)
               else
                  dLR = SQRT(Acell)
                  FACE_PG = U2D_GRID%CELL_EDGES_PG(intID,eleID)
                  if (GRID_BC(FACE_PG)%PARTICLE_BC == STATE) then ! ++++++++ GENERIC BOUNDARY +++++++++++++++++++
                     U_neigh = GRID_BC(FACE_PG)%U_BOUND(FIRST:LAST)
                     gradUprim_neigh = gradUprim(:,FIRST:LAST,eleID)
                  else if (GRID_BC(FACE_PG)%PARTICLE_BC == WALL) then ! ++++++++ WALL BOUNDARY ++++++++++++++++++++
                     call compute_wall_state(U(FIRST:LAST,eleID), &
                     nx, ny, U_wall, I)
                     U_neigh = U_wall
                     gradUprim_neigh = gradUprim(:,FIRST:LAST,eleID)
                  else if (GRID_BC(FACE_PG)%PARTICLE_BC == SYMMETRY) then ! ++++++++ SYM BOUNDARY ++++++++++++++++++++
                     call compute_sym_state(U(FIRST:LAST,eleID), &
                     nx, ny, U_sym, I)
                     U_neigh = U_sym
                     gradUprim_neigh = gradUprim(:,FIRST:LAST,eleID)
                  else
                     print*, "ERROR! UNKNOWN BOUNDARY TYPE ", neigh, " for element ", eleID, &
                     " Check the mesh or the pre-processing."
                     print*, "ABORTING!"
                     STOP
                  end if
               end if

               IF (FLUX_FUNCTION == AUSM) THEN
                  call compute_fluxes_AUSMplusup(U(FIRST:LAST,eleID), U_neigh, &
                  nx, ny, F_dot_n_hyper(FIRST:LAST), Acell, I)
               ELSE IF (FLUX_FUNCTION == HLL) THEN
                  call compute_fluxes_HLL(U(FIRST:LAST,eleID), U_neigh, &
                  nx, ny, F_dot_n_hyper(FIRST:LAST), Acell, I)
               ELSE IF (FLUX_FUNCTION == SLAU) THEN
                  call compute_fluxes_SLAU2(U(FIRST:LAST,eleID), U_neigh, &
                  nx, ny, F_dot_n_hyper(FIRST:LAST), Acell, I)
               ELSE
                  CALL ERROR_ABORT('Error! Flux function not available.')
               END IF

               call compute_fluxes_diffusive(U(FIRST:LAST,eleID), U_neigh, &
               gradUprim(:,FIRST:LAST,eleID), gradUprim_neigh, &
               nx, ny, F_dot_n_diff(FIRST:LAST), Acell, dLR, I)
               !F_dot_n_diff = 0.d0

            END DO


            F_dot_n_wall = 0.d0
            IF (neigh == -1) THEN ! +++++++++ BOUNDARY CELL
               dLR = SQRT(Acell)
               FACE_PG = U2D_GRID%CELL_EDGES_PG(intID,eleID)
               IF (GRID_BC(FACE_PG)%REACT) THEN
                  !WRITE(*,*) F_dot_n_hyper
                  F_dot_n_hyper = 0.d0
                  !WRITE(*,*) 'Set zero flux for cell ', eleID, ' face ', intID
                  CALL compute_wall_fluxes(U(:,eleID), nx, ny, F_dot_n_wall)
                  !WRITE(*,*) F_dot_n_hyper
               END IF
            END IF
            !F_dot_n_wall = 0.d0


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
                  stop
               end if
            end do

         end do ! End loop on interfaces

         ! Add source terms
         call compute_source_terms(U(:,eleID), S, eleID)

         U_new(:,eleID) = U_new(:,eleID) + dt*S

   
      end do ! End loop on elements
   
      U = U_new ! Save results

      DEALLOCATE(F_dot_n_hyper)
      DEALLOCATE(F_dot_n_diff)
      DEALLOCATE(S)
      DEALLOCATE(Uprim)
      DEALLOCATE(gradUprim)


   end subroutine


   subroutine compute_source_terms(U, S, eleID)

      ! This function computes the source terms of the equations
      ! These originate from reactive processes and intra-species collisions
      ! A term is also added when the equations are axisymmetric

      implicit none

      real(kind=8), dimension(:),   intent(in)  :: U
      real(kind=8), dimension(:),   intent(out) :: S
      INTEGER,   intent(in)                     :: eleID
      INTEGER :: I, J, JR, R1_SP_ID, R2_SP_ID, P, P_SP_ID
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: prim
      INTEGER :: FIRST, LAST, RHOI, RHOJ, MOMXI, MOMXJ, MOMYI, MOMYJ, ENEI, ENEJ, RHOP, MOMXP, MOMYP, ENEP
      REAL(KIND=8) :: MI, MJ, MIJ, NI, NJ, UXI, UXJ, UYI, UYJ, TI, TJ, TIJ, VTHIJ, S_MOM_IJ, UIDOTW, UJDOTW, S_ENE_IJ, QIJ
      REAL(KIND=8) :: K_FORWARD, TEMP, RATE_OF_PROGRESS
      REAL(KIND=8) :: MP, NP, UXP, UYP, EP, TP, EI

      ALLOCATE(prim(N_SPECIES_FLUID*Neq))
      S = 0.d0
      
      DO I = 1, N_SPECIES_FLUID
         FIRST = (I-1)*Neq+1
         LAST = I*Neq+1
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

               S_MOM_IJ = 4./3.*MIJ*VTHIJ*QIJ*NI*NJ
               S(MOMXI) = S(MOMXI) + S_MOM_IJ*( UXJ - UXI )
               S(MOMXJ) = S(MOMXJ) - S_MOM_IJ*( UXJ - UXI )
               S(MOMYI) = S(MOMYI) + S_MOM_IJ*( UYJ - UYI )
               S(MOMYJ) = S(MOMYJ) - S_MOM_IJ*( UYJ - UYI )

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


   SUBROUTINE compute_wall_fluxes(U, nx, ny, F_dot_n_wall)

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
         LAST = I*Neq+1

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
         LAST = I*Neq+1
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
      real(kind=8) :: ws_min_L, ws_max_L, ws_min_R, ws_max_R, ws_max

      call compute_flux_ws(U_L, F_L, nx, ny, ws_max_L, ws_min_L, SP_ID)
      call compute_flux_ws(U_R, F_R, nx, ny, ws_max_R, ws_min_R, SP_ID)

      ! ws_min = MIN(ws_min_L, ws_min_R)
      ! ws_max = MAX(ws_max_L, ws_max_R)

      ! HLL fluxes, check AUSM+-UP/preconditioning
      ! if (ws_min .ge. 0.0) then
      !    F_dot_n = F_L
      ! else if (ws_max .lt. 0.0) then
      !    F_dot_n = F_R
      ! else
      !    F_dot_n = (ws_min*ws_max*(U_R - U_L) + ws_max*F_L - ws_min*F_R)/(ws_max - ws_min)
      ! end if


      ws_max  = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))

      F_dot_n = 0.5*(F_R + F_L) - ws_max/2.0*(U_R-U_L) ! Rusanov flux
   
   

      ! Update global maximum wave speed (used for setting the time step)
      ws_max = abs(ws_max)
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


      Mface = MLP + MRM + MP

      pflux = betaLP*pL + betaRM*pR + Pu

      

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

      pflux = pbar !+ (0.5*deltabeta*(pL-pR)) + SQRT(0.5*(uL*uL + vL*vL + uR*uR + vR*vR)) * &
                    !    (betaL + betaR - 1.0) * rhobar * abar
      
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


   subroutine compute_fluxes_diffusive(U_L, U_R, gradU_L, gradU_R, nx, ny, F_dot_n, A_ele, dLR, SP_ID)

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
      INTEGER,                      intent(in)  :: SP_ID

      real(kind=8), dimension(Neq) :: prim

      real(kind=8) :: ux_L, ux_R, uy_L, uy_R, T_L, T_R, rho_L
      real(kind=8) :: DUXDX, DUYDX, DUXDY, DUYDY, TAUXX, TAUXY, TAUYY, MU, DTDX, DTDY, UX, UY


      F_dot_n = 0.d0

      rho_L = U_L(1)

      DUXDX = 0.5*(gradU_L(1,2) + gradU_R(1,2))
      DUYDX = 0.5*(gradU_L(1,3) + gradU_R(1,3))
      DUXDY = 0.5*(gradU_L(2,2) + gradU_R(2,2))
      DUYDY = 0.5*(gradU_L(2,3) + gradU_R(2,3))

      MU = SPECIES(SP_ID)%MU
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

      F_dot_n(4) = - SPECIES(SP_ID)%KAPPA*(T_R-T_L)/dLR

      !F_dot_n(4) = - (NX*DTDX + NY*DTDY)*SPECIES(SP_ID)%KAPPA &
      !             - (UX*TAUXX + UY*TAUXY)*NX &
      !             - (UX*TAUXY + UY*TAUYY)*NY

      F_dot_n(4) = F_dot_n(4) - (UX*TAUXX + UY*TAUXY)*NX &
                              - (UX*TAUXY + UY*TAUYY)*NY

      invdt_cond = MAX(invdt_cond, 6.*SPECIES(SP_ID)%KAPPA/(A_ele*rho_L*SPECIES(SP_ID)%CP))
      invdt_diff = MAX(invdt_diff, 6.*SPECIES(SP_ID)%MU/(A_ele*rho_L))



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
               if (GRID_BC(FACE_PG)%PARTICLE_BC == STATE) then ! ++++++++ GENERIC BOUNDARY +++++++++++++++++++
                  DO SP_ID = 1, N_SPECIES_FLUID
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq+1
                     CALL compute_primitive_from_conserved(GRID_BC(FACE_PG)%U_BOUND(FIRST:LAST), U_adj(FIRST:LAST), SP_ID)
                  END DO
                  Uface = 0.5*(U(:,I) + U_adj)
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == WALL) then ! ++++++++ WALL BOUNDARY ++++++++++++++++++++
                  U_adj = U(:,I)
                  DO SP_ID = 1, N_SPECIES_FLUID
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq+1
                     U_adj(FIRST+1) = 0.d0
                     U_adj(FIRST+2) = 0.d0
                     U_adj(FIRST+3) = Tw
                  END DO
                  Uface = 0.5*(U(:,I) + U_adj)
                  !Uface = U(:,I)
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == SYMMETRY) then ! ++++++++ SYM BOUNDARY ++++++++++++++++++++
                  ONAXIS = .TRUE.
                  U_adj = U(:,I)
                  DO SP_ID = 1, N_SPECIES_FLUID
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq+1
                     U_adj(FIRST+1) = U(FIRST+1,I) - 2.0*(U(FIRST+1,I)*nx + U(FIRST+2,I)*ny)*nx
                     U_adj(FIRST+2) = U(FIRST+2,I) - 2.0*(U(FIRST+1,I)*nx + U(FIRST+2,I)*ny)*ny
                  END DO
                  Uface = 0.5*(U(:,I) + U_adj)
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
               if (GRID_BC(FACE_PG)%PARTICLE_BC == STATE) then ! ++++++++ GENERIC BOUNDARY +++++++++++++++++++
                  DO SP_ID = 1, N_SPECIES_FLUID
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq+1
                     CALL compute_primitive_from_conserved(GRID_BC(FACE_PG)%U_BOUND(FIRST:LAST), U_adj(FIRST:LAST), SP_ID)
                  END DO
                  Uface = 0.5*(U(:,I) + U_adj)
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == WALL) then ! ++++++++ WALL NO-SLIP BOUNDARY ++++++++++++++++++++
                  U_adj = U(:,I)
                  DO SP_ID = 1, N_SPECIES_FLUID
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq+1
                     U_adj(FIRST+1) = 0.d0
                     U_adj(FIRST+2) = 0.d0
                     U_adj(FIRST+3) = Tw
                  END DO
                  Uface = 0.5*(U(:,I) + U_adj)
                  !Uface = U(:,I)
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == SYMMETRY) then ! ++++++++ SYM BOUNDARY ++++++++++++++++++++
                  ONAXIS = .TRUE.
                  U_adj = U(:,I)
                  DO SP_ID = 1, N_SPECIES_FLUID
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq+1
                     U_adj(FIRST+1) = U(FIRST+1,I) - 2.0*(U(FIRST+1,I)*nx + U(FIRST+2,I)*ny)*nx
                     U_adj(FIRST+2) = U(FIRST+2,I) - 2.0*(U(FIRST+1,I)*nx + U(FIRST+2,I)*ny)*ny
                  END DO
                  Uface = 0.5*(U(:,I) + U_adj)
               else
                  print*, "ERROR! UNKNOWN BOUNDARY TYPE ", neigh, " for element ", I, &
                  " Check the mesh or the pre-processing."
                  print*, "ABORTING!"
                  STOP
               end if
            end if

            gradU(1,:,I) = gradU(1,:,I) + (U_adj - U(:,I))*U2D_GRID%LSTSQ_COEFFS(1,J,I)
            gradU(2,:,I) = gradU(2,:,I) + (U_adj - U(:,I))*U2D_GRID%LSTSQ_COEFFS(2,J,I)

         END DO

         !IF (ONAXIS) gradU(2,:,I) = 0.d0

      END DO

      DEALLOCATE(Uface)
      DEALLOCATE(U_adj)

   end subroutine

end module
