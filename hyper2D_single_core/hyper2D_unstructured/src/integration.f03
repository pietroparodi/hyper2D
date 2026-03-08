module integration

   use pde
   use grid
   use global_module

   implicit none
   
   contains

   ! ======================================================================== 
   
   subroutine forward_Euler_step(U, U_new, dt)

      ! This function performs one step of the Forward Euler explicit time integrator 

      implicit none
   
      real(kind=8), dimension(:,:), intent(inout)  :: U, U_new

      real(kind=8), intent(in) :: dt
      integer :: eleID, eqID, intID
      
      real(kind=8), dimension(Neq) :: U_sym, U_wall, U_neigh
      real(kind=8), dimension(:), allocatable :: F_dot_n_hyper, F_dot_n_diff, S
      real(kind=8) :: nx, ny, Lint, Aele, dLR
      integer      :: neigh, FACE_PG, I, FIRST, LAST
      real(kind=8), dimension(3) :: C1, C2

      real(kind=8), dimension(:,:), ALLOCATABLE :: Uprim
      real(kind=8), dimension(:,:,:), ALLOCATABLE :: gradUprim

      ALLOCATE(F_dot_n_hyper(NSPECIES*Neq))
      ALLOCATE(F_dot_n_diff(NSPECIES*Neq))
      ALLOCATE(S(NSPECIES*Neq))
      ALLOCATE(Uprim(NSPECIES*Neq,NCELLS))
      ALLOCATE(gradUprim(2,NSPECIES*Neq,NCELLS))

      DO eleID = 1, NCELLS
         DO I = 1, NSPECIES
            FIRST = (I-1)*Neq+1
            LAST = I*Neq+1
            call compute_primitive_from_conserved(U(FIRST:LAST,eleID), Uprim(FIRST:LAST,eleID), I)
         END DO
      END DO

      call compute_cell_centered_gradients_green_gauss(Uprim, gradUprim)

      do eleID = 1, NCELLS

         Aele = U2D_GRID%CELL_AREAS(eleID)
         U_new(:,eleID) = U(:,eleID) ! Init

         do intID = 1, 3 ! Only quad element supported

            F_dot_n_hyper = 0.0 ! Init
            F_dot_n_diff = 0.0 ! Init

            ! Extract data
            Lint = U2D_GRID%CELL_EDGES_LEN(intID,eleID)
            nx = U2D_GRID%EDGE_NORMAL(1,intID,eleID)
            ny = U2D_GRID%EDGE_NORMAL(2,intID,eleID)

            neigh = U2D_GRID%CELL_NEIGHBORS(intID,eleID)

            ! Check what neighbor is it
            DO I = 1, NSPECIES
               FIRST = (I-1)*Neq+1
               LAST = I*Neq+1

               if (neigh .NE. -1) then ! +++++++++ INTERNAL CELL

                  C1 = U2D_GRID%CELL_CENTROIDS(:,eleID)
                  C2 = U2D_GRID%CELL_CENTROIDS(:,neigh)
                  dLR = NORM2(C2-C1)
                  U_neigh = U(FIRST:LAST,neigh)
               else
                  dLR = SQRT(Aele)
                  FACE_PG = U2D_GRID%CELL_EDGES_PG(intID,eleID)
                  if (GRID_BC(FACE_PG)%PARTICLE_BC == INLET) then ! ++++++++ INLET BOUNDARY +++++++++++++++++++
                     U_neigh = U_inlet(FIRST:LAST)
                  else if (GRID_BC(FACE_PG)%PARTICLE_BC == OUTLET) then ! ++++++++ OUTLET BOUNDARY +++++++++++++++++++
                     U_neigh = U_outlet(FIRST:LAST)
                  else if (GRID_BC(FACE_PG)%PARTICLE_BC == WALL) then ! ++++++++ WALL BOUNDARY ++++++++++++++++++++
                     call compute_wall_state(U(FIRST:LAST,eleID), &
                     nx, ny, U_wall, I)
                     U_neigh = U_wall
                  else if (GRID_BC(FACE_PG)%PARTICLE_BC == SYMMETRY) then ! ++++++++ SYM BOUNDARY ++++++++++++++++++++
                     call compute_sym_state(U(FIRST:LAST,eleID), &
                     nx, ny, U_sym, I)
                     U_neigh = U_sym
                  else
                     print*, "ERROR! UNKNOWN BOUNDARY TYPE ", neigh, " for element ", eleID, &
                     " Check the mesh or the pre-processing."
                     print*, "ABORTING!"
                     STOP
                  end if
               end if

               call compute_fluxes_HLL(U(FIRST:LAST,eleID), U_neigh, &
               nx, ny, F_dot_n_hyper(FIRST:LAST), Aele, I)
               call compute_fluxes_diffusive(U(FIRST:LAST,eleID), U_neigh, &
               gradUprim(:,FIRST:LAST,eleID), gradUprim(:,FIRST:LAST,eleID), &
               nx, ny, F_dot_n_diff(FIRST:LAST), Aele, dLR, I)

            END DO
            ! Update solution
            !U_new(:,eleID) = U_new(:,eleID) - dt*(F_dot_n_hyper + F_dot_n_diff)*Lint/Aele
            U_new(:,eleID) = U_new(:,eleID) - dt*(F_dot_n_hyper)*Lint/Aele

            ! Check that the solution did not diverge
            do eqID = 1, NSPECIES*Neq
               if (isnan(U_new(eqID,eleID))) then 
                  print*, 'Solution diverged, try with a smaller time step! Aborting.'
                  print*, 'Solution that diverged: ', U_new(:,eleID)
                  print*, 'in cell ID = ', eleID
                  stop
               end if
            end do

         end do ! End loop on interfaces

         ! Add source terms
         !call compute_source_terms(U(:,eleID), Aele, S)
         !WRITE(*,*) S
         !WRITE(*,*) ' ' 
         !U_new(:,eleID) = U_new(:,eleID) + dt*S
   
      end do ! End loop on elements
   
      U = U_new ! Save results

      DEALLOCATE(F_dot_n_hyper)
      DEALLOCATE(F_dot_n_diff)
      DEALLOCATE(S)
      DEALLOCATE(Uprim)
      DEALLOCATE(gradUprim)

   end subroutine




   ! ======================================================================== 

   subroutine compute_source_terms(U, S)

      implicit none

      real(kind=8), dimension(:),   intent(in)  :: U
      real(kind=8), dimension(:),   intent(out) :: S
      INTEGER :: I, J
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: prim
      INTEGER :: FIRST, LAST, RHOI, RHOJ, MOMXI, MOMXJ, MOMYI, MOMYJ, ENEI, ENEJ
      REAL(KIND=8) :: MI, MJ, MIJ, NI, NJ, UXI, UXJ, UYI, UYJ, TI, TJ, TIJ, VTHIJ, S_MOM_IJ, UIDOTW, UJDOTW, S_ENE_IJ, QIJ

      ALLOCATE(prim(NSPECIES*Neq))
      S = 0.d0
      
      DO I = 1, NSPECIES
         FIRST = (I-1)*Neq+1
         LAST = I*Neq+1
         call compute_primitive_from_conserved(U(FIRST:LAST), prim(FIRST:LAST), I)
      END DO

      DO I = 1, NSPECIES
         DO J = I+1, NSPECIES
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
            !S(MOMXI) = S(MOMXI) + S_MOM_IJ*( UXJ - UXI )
            !S(MOMXJ) = S(MOMXJ) - S_MOM_IJ*( UXJ - UXI )
            !S(MOMYI) = S(MOMYI) + S_MOM_IJ*( UYJ - UYI )
            !S(MOMYJ) = S(MOMYJ) - S_MOM_IJ*( UYJ - UYI )

            UIDOTW = UXI*(UXI-UXJ) + UYI*(UYI-UYJ)
            UJDOTW = UXJ*(UXI-UXJ) + UYJ*(UYI-UYJ)
            S_ENE_IJ = 4./3.*MIJ/(MI+MJ)*VTHIJ*QIJ*NI*NJ*(3*KB*(TI-TJ) + MI*TJ/TIJ*UIDOTW + MJ*TI/TIJ*UJDOTW)
            !S(ENEI) = S(ENEI) - S_ENE_IJ
            !S(ENEJ) = S(ENEJ) + S_ENE_IJ

            ! Update stability constraints
            !ws_over_sqrtA_maxabs = MAX(ws_over_sqrtA_maxabs, 4./3.*NI*QIJ*VTHIJ, 4./3.*NJ*QIJ*VTHIJ)
         END DO
      END DO

      DEALLOCATE(prim)

   end subroutine


   ! ======================================================================== 

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
      real(kind=8) :: ws_min_L, ws_max_L, ws_min_R, ws_max_R, ws_max, ws_min

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
      ws_over_sqrtA_maxabs = MAX(ws_over_sqrtA_maxabs, ws_max/sqrt(A_ele))

   end subroutine



   subroutine compute_fluxes_AUSMplusup(U_L, U_R, nx, ny, flux, A_ele, SP_ID) !(U_L, U_R, nx, ny, F_dot_n, A_ele, dLR, SP_ID)
      
      implicit none


      real(kind=8), dimension(:),   intent(in)  :: U_L, U_R
      real(kind=8),                 intent(in)  :: nx, ny
      real(kind=8), dimension(Neq), intent(out) :: flux

      INTEGER,                      intent(in)  :: SP_ID
      real(kind=8),                 intent(in)  :: A_ele


      real(kind=8) :: rhoL,rhoR,uL,vL,uR,vR
      real(kind=8) :: pL,pR,EL,ER,HL,HR
      real(kind=8) :: unL,unR
      real(kind=8) :: aL,aR,a
      real(kind=8) :: ML,MR
      real(kind=8) :: MpL,MmR
      real(kind=8) :: ppL,pmR
      real(kind=8) :: mdot,pface
      real(kind=8) :: Mbar2
      real(kind=8) :: alpha,beta,Kp,Ku
      real(kind=8) :: u,v,H
      real(kind=8) :: pav
      real(kind=8) :: gamma

      alpha = 3.0d0/16.0d0
      beta  = 1.0d0/8.0d0
      Kp    = 0.25d0
      Ku    = 0.75d0

      gamma = SPECIES(SP_ID)%GAMMA
      !--------------------------------
      ! Left primitive variables
      !--------------------------------

      rhoL = U_L(1)
      uL = U_L(2)/rhoL
      vL = U_L(3)/rhoL
      EL = U_L(4)/rhoL

      pL = (gamma-1.d0)*rhoL*(EL - 0.5d0*(uL*uL+vL*vL))
      HL = EL + pL/rhoL

      !--------------------------------
      ! Right primitive variables
      !--------------------------------

      rhoR = U_R(1)
      uR = U_R(2)/rhoR
      vR = U_R(3)/rhoR
      ER = U_R(4)/rhoR

      pR = (gamma-1.d0)*rhoR*(ER - 0.5d0*(uR*uR+vR*vR))
      HR = ER + pR/rhoR

      !--------------------------------
      ! Normal velocities
      !--------------------------------

      unL = uL*nx + vL*ny
      unR = uR*nx + vR*ny

      !--------------------------------
      ! Sound speeds
      !--------------------------------

      aL = sqrt(gamma*pL/rhoL)
      aR = sqrt(gamma*pR/rhoR)

      a = max(aL,aR)

      ML = unL/a
      MR = unR/a

      !--------------------------------
      ! Mach splitting
      !--------------------------------

      if (abs(ML) >= 1.d0) then
         MpL = 0.5d0*(ML + abs(ML))
      else
         MpL = 0.25d0*(ML+1.d0)**2 + beta*(ML*ML-1.d0)**2
      end if

      if (abs(MR) >= 1.d0) then
         MmR = 0.5d0*(MR - abs(MR))
      else
         MmR = -0.25d0*(MR-1.d0)**2 - beta*(MR*MR-1.d0)**2
      end if

      !--------------------------------
      ! Pressure splitting
      !--------------------------------

      if (abs(ML) >= 1.d0) then
         ppL = 0.5d0*(1.d0 + sign(1.d0,ML))
      else
         ppL = 0.25d0*(ML+1.d0)**2*(2.d0-ML) + alpha*ML*(ML*ML-1.d0)**2
      end if

      if (abs(MR) >= 1.d0) then
         pmR = 0.5d0*(1.d0 - sign(1.d0,MR))
      else
         pmR = 0.25d0*(MR-1.d0)**2*(2.d0+MR) - alpha*MR*(MR*MR-1.d0)**2
      end if

      !--------------------------------
      ! Mass flux
      !--------------------------------

      Mbar2 = 0.5d0*(ML*ML + MR*MR)

      mdot = a*(MpL*rhoL + MmR*rhoR) + &
            Kp*max(1.d0-Mbar2,0.d0)*(pR-pL)/a

      !--------------------------------
      ! Interface pressure
      !--------------------------------

      pav = 0.5d0*(pL+pR)

      pface = ppL*pL + pmR*pR + Ku*pav*(MR-ML)

      !--------------------------------
      ! Upwind state
      !--------------------------------

      if (mdot > 0.d0) then
         u = uL
         v = vL
         H = HL
      else
         u = uR
         v = vR
         H = HR
      end if

      !--------------------------------
      ! Flux vector
      !--------------------------------

      flux(1) = mdot
      flux(2) = mdot*u + pface*nx
      flux(3) = mdot*v + pface*ny
      flux(4) = mdot*H

      ! Update global maximum wave speed (used for setting the time step)
      ws_over_sqrtA_maxabs = MAX(ws_over_sqrtA_maxabs, a/sqrt(A_ele))

   end subroutine

   subroutine compute_fluxes_diffusive(U_L, U_R, gradU_L, gradU_R, nx, ny, F_dot_n, A_ele, dLR, SP_ID)

      ! Computes HLL numerical fluxes among the cell eleID and the neighbor cell neigh
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

      ! Update global maximum wave speed (used for setting the time step)
      ws_over_sqrtA_maxabs = MAX(ws_over_sqrtA_maxabs, &
      SPECIES(SP_ID)%KAPPA/(A_ele*rho_L*SPECIES(SP_ID)%CP))

      DUXDX = 0.5*(gradU_L(1,2) + gradU_R(1,2))
      DUYDX = 0.5*(gradU_L(1,3) + gradU_R(1,3))
      DUXDY = 0.5*(gradU_L(2,2) + gradU_R(2,2))
      DUYDY = 0.5*(gradU_L(2,3) + gradU_R(2,3))

      MU = SPECIES(SP_ID)%MU
      TAUXX = 2./3.*MU*(2.*DUXDX-DUYDY)
      TAUYY = 2./3.*MU*(2.*DUYDY-DUXDX)
      TAUXY = MU*(DUXDY + DUYDX)

      F_dot_n(2) = F_dot_n(2) - (NX*TAUXX + NY*TAUXY)
      F_dot_n(3) = F_dot_n(3) - (NX*TAUXY + NY*TAUYY)



      DTDX = 0.5*(gradU_L(1,4) + gradU_R(1,4))
      DTDY = 0.5*(gradU_L(2,4) + gradU_R(2,4))

      UX = 0.5*(ux_L + ux_R)
      UY = 0.5*(ux_L + ux_R)

      F_dot_n(4) = F_dot_n(4) - (NX*DTDX + NY*DTDY)*SPECIES(SP_ID)%KAPPA &
                 - (UX*TAUXX + UY*TAUXY)*NX &
                 - (UX*TAUXY + UY*TAUYY)*NY

   end subroutine



   subroutine compute_cell_centered_gradients(U, gradU)

      implicit none
   
      real(kind=8), dimension(:,:), intent(in)  :: U
      real(kind=8), dimension(:,:,:), intent(inout)  :: gradU

      real(kind=8), dimension(:,:), allocatable  :: Unodes
      real(kind=8), dimension(:), allocatable  :: TOTALW

      INTEGER :: I, J, INODE, IDIR
      REAL(KIND=8) :: WEIGHT

      ALLOCATE(Unodes(NSPECIES*Neq, NNODES))
      ALLOCATE(TOTALW(NNODES))
      
      Unodes = 0.d0
      TOTALW = 0.d0

      ! Interpolate from cell centers to cell nodes

      DO I = 1, NCELLS
         DO J = 1, 3
            INODE = U2D_GRID%CELL_NODES(J,I)
            WEIGHT = 1./NORM2(U2D_GRID%NODE_COORDS(:,INODE) - U2D_GRID%CELL_CENTROIDS(:,I))
            Unodes(:,INODE) = Unodes(:,INODE) + U(:,I)*WEIGHT
            TOTALW(INODE) = TOTALW(INODE) + WEIGHT
         END DO
      END DO

      DO I = 1, NNODES
         Unodes(:,I) = Unodes(:,I)/TOTALW(I)
      END DO


      ! ! Test test
      ! DO I = 1, NCELLS
      !    DO J = 1, 3
      !       IF (U2D_GRID%CELL_EDGES_PG(J,I) .NE. -1) THEN
      !          IF (GRID_BC(U2D_GRID%CELL_EDGES_PG(J,I))%PARTICLE_BC == WALL) THEN
      !             IF (J == 1) THEN
      !                INODE = U2D_GRID%CELL_NODES(1,I)
      !                Unodes(2,INODE) = 0.d0
      !                Unodes(3,INODE) = 0.d0
      !                INODE = U2D_GRID%CELL_NODES(2,I)
      !                Unodes(2,INODE) = 0.d0
      !                Unodes(3,INODE) = 0.d0
      !             ELSE IF (J == 2) THEN
      !                INODE = U2D_GRID%CELL_NODES(2,I)
      !                Unodes(2,INODE) = 0.d0
      !                Unodes(3,INODE) = 0.d0
      !                INODE = U2D_GRID%CELL_NODES(3,I)
      !                Unodes(2,INODE) = 0.d0
      !                Unodes(3,INODE) = 0.d0
      !             ELSE IF (J == 3) THEN
      !                INODE = U2D_GRID%CELL_NODES(3,I)
      !                Unodes(2,INODE) = 0.d0
      !                Unodes(3,INODE) = 0.d0
      !                INODE = U2D_GRID%CELL_NODES(1,I)
      !                Unodes(2,INODE) = 0.d0
      !                Unodes(3,INODE) = 0.d0
      !             END IF
      !          END IF
      !       END IF
      !    END DO
      ! END DO
      ! ! End test test

      DEALLOCATE(TOTALW)

      ! Compute gradient in each cell from nodal values
      gradU = 0.d0

      DO I = 1, NCELLS
         DO IDIR = 1, 2 ! The coordinate direction (x, y)
            DO J = 1, 3 ! The cell node
               INODE = U2D_GRID%CELL_NODES(J,I)
               gradU(IDIR,:,I) = gradU(IDIR,:,I) + Unodes(:,INODE) * U2D_GRID%BASIS_COEFFS(IDIR,J,I)
            END DO
         END DO
      END DO

      DEALLOCATE(Unodes)

   end subroutine





   subroutine compute_cell_centered_gradients_green_gauss(U, gradU)

      implicit none
   
      real(kind=8), dimension(:,:), intent(in)  :: U
      real(kind=8), dimension(:,:,:), intent(inout)  :: gradU

      INTEGER :: I, J, neigh, FACE_PG, SP_ID, FIRST, LAST
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Uface
      REAL(KIND=8) :: Lint, nx, ny, Aele
      REAL(KIND=8), DIMENSION(Neq) :: prim
      REAL(KIND=8), DIMENSION(:), allocatable :: U_adj



      ! Compute gradient in each cell from nodal values
      gradU = 0.d0
      ALLOCATE(Uface(NSPECIES*Neq))
      ALLOCATE(U_adj(NSPECIES*Neq))

      DO I = 1, NCELLS
         DO J = 1, 3 ! The cell face
            ! Extract data
            Lint = U2D_GRID%CELL_EDGES_LEN(J,I)
            nx = U2D_GRID%EDGE_NORMAL(1,J,I)
            ny = U2D_GRID%EDGE_NORMAL(2,J,I)
            Aele = U2D_GRID%CELL_AREAS(I)

            neigh = U2D_GRID%CELL_NEIGHBORS(J,I)
            IF (neigh .NE. -1) THEN
               Uface = 0.5*(U(:,I) + U(:,neigh))

            else
               FACE_PG = U2D_GRID%CELL_EDGES_PG(J,I)
               if (GRID_BC(FACE_PG)%PARTICLE_BC == INLET) then ! ++++++++ INLET BOUNDARY +++++++++++++++++++
                  DO SP_ID = 1, NSPECIES
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq+1
                     call compute_primitive_from_conserved(U_inlet(FIRST:LAST), U_adj(FIRST:LAST), SP_ID)
                  END DO
                  Uface = 0.5*(U(:,I) + U_adj)
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == OUTLET) then ! ++++++++ OUTLET BOUNDARY +++++++++++++++++++
                  DO SP_ID = 1, NSPECIES
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq+1
                     call compute_primitive_from_conserved(U_outlet(FIRST:LAST), U_adj(FIRST:LAST), SP_ID)
                  END DO
                  Uface = 0.5*(U(:,I) + U_adj)
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == WALL) then ! ++++++++ WALL BOUNDARY ++++++++++++++++++++
                  U_adj = U(:,I)
                  DO SP_ID = 1, NSPECIES
                     FIRST = (SP_ID-1)*Neq+1
                     LAST = SP_ID*Neq+1
                     U_adj(FIRST+1) = 0.d0
                     U_adj(FIRST+2) = 0.d0
                     U_adj(FIRST+3) = Tw
                  END DO
                  Uface = 0.5*(U(:,I) + U_adj)
                  !Uface = U(:,I)
               else if (GRID_BC(FACE_PG)%PARTICLE_BC == SYMMETRY) then ! ++++++++ SYM BOUNDARY ++++++++++++++++++++
                  U_adj = U(:,I)
                  DO SP_ID = 1, NSPECIES
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

            gradU(1,:,I) = gradU(1,:,I) + Uface*nx*Lint/Aele
            gradU(2,:,I) = gradU(2,:,I) + Uface*ny*Lint/Aele

         END DO

      END DO

      DEALLOCATE(Uface)
      DEALLOCATE(U_adj)

   end subroutine

end module
