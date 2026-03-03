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
      
      real(kind=8), dimension(Neq) :: U_sym
      real(kind=8), dimension(:), allocatable :: F_dot_n, S
      real(kind=8) :: nx, ny, Lint, Aele, dLR
      integer      :: neigh, FACE_PG, I, FIRST, LAST
      real(kind=8), dimension(3) :: C1, C2

      ALLOCATE(F_dot_n(NSPECIES*Neq))
      ALLOCATE(S(NSPECIES*Neq))

      do eleID = 1, NCELLS


         Aele = U2D_GRID%CELL_AREAS(eleID)
         U_new(:,eleID) = U(:,eleID) ! Init

         do intID = 1, 3 ! Only quad element supported

            F_dot_n = 0.0 ! Init

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

                  call compute_fluxes_HLL_cond(U(FIRST:LAST,eleID), U(FIRST:LAST,neigh), &
                  nx, ny, F_dot_n(FIRST:LAST), Aele, dLR, I)

               else
                  dLR = SQRT(Aele)
                  FACE_PG = U2D_GRID%CELL_EDGES_PG(intID,eleID)
                  if (GRID_BC(FACE_PG)%PARTICLE_BC == INLET) then ! ++++++++ INLET BOUNDARY +++++++++++++++++++
                     call compute_fluxes_HLL_cond(U(FIRST:LAST,eleID), U_inlet(FIRST:LAST), &
                     nx, ny, F_dot_n(FIRST:LAST), Aele, dLR, I)
                  else if (GRID_BC(FACE_PG)%PARTICLE_BC == OUTLET) then ! ++++++++ OUTLET BOUNDARY +++++++++++++++++++
                     call compute_fluxes_HLL_cond(U(FIRST:LAST,eleID), U_outlet(FIRST:LAST), &
                     nx, ny, F_dot_n(FIRST:LAST), Aele, dLR, I)
                  else if (GRID_BC(FACE_PG)%PARTICLE_BC == WALL) then ! ++++++++ WALL BOUNDARY ++++++++++++++++++++
                     call compute_wall_flux(U(FIRST:LAST,eleID), &
                     nx, ny, F_dot_n(FIRST:LAST), Aele, dLR, Tw, I)
                  else if (GRID_BC(FACE_PG)%PARTICLE_BC == SYMMETRY) then ! ++++++++ SYM BOUNDARY ++++++++++++++++++++
                     call compute_sym_state(U(FIRST:LAST,eleID), &
                     nx, ny, U_sym, I)
                     call compute_fluxes_HLL_cond(U(FIRST:LAST,eleID), U_sym, &
                     nx, ny, F_dot_n(FIRST:LAST), Aele, dLR, I)
                  else
                     print*, "ERROR! UNKNOWN BOUNDARY TYPE ", neigh, " for element ", eleID, &
                     " Check the mesh or the pre-processing."
                     print*, "ABORTING!"
                     STOP
                  end if
               end if

            END DO
            ! Update solution
            U_new(:,eleID) = U_new(:,eleID) - dt*F_dot_n*Lint/Aele

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
         call compute_source_terms(U(:,eleID), Aele, S)
         U_new(:,eleID) = U_new(:,eleID) + dt*S
   
      end do ! End loop on elements
   
      U = U_new ! Save results

      DEALLOCATE(F_dot_n)
      DEALLOCATE(S)

   end subroutine




   ! ======================================================================== 

   subroutine compute_source_terms(U, A_ele, S)

      implicit none

      real(kind=8), dimension(:),   intent(in)  :: U
      real(kind=8),                 intent(in)  :: A_ele
      real(kind=8), dimension(:),   intent(out) :: S
      INTEGER :: I, J

      DO I = 1, NSPECIES
         DO J = I+1, NSPECIES
            MI = SPECIES(I)%MOLECULAR_MASS
            MJ = SPECIES(J)%MOLECULAR_MASS
            MIJ = MI*MJ/(MI+MJ)
            S_mom_ij = 16/3*MIJ/(MI+MJ)*OMEGA11(I,J)*NI*NJ
            S((I-1)*Neq+2) = 
         END DO
      END DO

      S = 0.d0


   end subroutine


   ! ======================================================================== 

   subroutine compute_fluxes_HLL_cond(U_L, U_R, nx, ny, F_dot_n, A_ele, dLR, SP_ID)

      ! Computes HLL numerical fluxes among the cell eleID and the neighbor cell neigh
      ! The element area Aele is also passed, for the sake of computing the CFL number.

      implicit none

      real(kind=8), dimension(:),   intent(in)  :: U_L, U_R
      real(kind=8),                 intent(in)  :: nx, ny, A_ele, dLR
      real(kind=8), dimension(Neq), intent(out) :: F_dot_n
      INTEGER,                      intent(in)  :: SP_ID

      real(kind=8), dimension(Neq) :: F_L, F_R, prim

      real(kind=8) :: T_L, T_R, rho

      ! Wave speeds
      real(kind=8) :: ws_min_L, ws_max_L, ws_min_R, ws_max_R, ws_max, ws_min

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
         F_dot_n = (ws_min*ws_max*(U_R - U_L) + ws_max*F_L - ws_min*F_R)/(ws_max - ws_min);
      end if


    
      ! Compute heat diffusion
      call compute_primitive_from_conserved(U_L, prim, SP_ID)
      T_L   = prim(4)
      rho = prim(1)
      call compute_primitive_from_conserved(U_R, prim, SP_ID)
      T_R   = prim(4)

      F_dot_n(4) = F_dot_n(4) + SPECIES(SP_ID)%KAPPA*(T_R-T_L)/dLR
      

      ! Update global maximum wave speed (used for setting the time step)
      ws_max = abs(ws_max)
      ws_over_sqrtA_maxabs = MAX(ws_over_sqrtA_maxabs, ws_max/sqrt(A_ele), &
      SPECIES(SP_ID)%KAPPA/(2*A_ele*rho*SPECIES(SP_ID)%CP))

   end subroutine



end module
