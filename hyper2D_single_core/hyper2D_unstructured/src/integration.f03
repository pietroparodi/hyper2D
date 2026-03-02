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
      
      real(kind=8), dimension(Neq) :: F_dot_n, U_sym
      real(kind=8) :: nx, ny, Lint, Aele
      integer      :: neigh, FACE_PG

      do eleID = 1, NCELLS

         !Aele = ele_geom(eleID, 1)
         Aele = U2D_GRID%CELL_AREAS(eleID)
         U_new(:,eleID) = U(:,eleID) ! Init

         do intID = 1, 3 ! Only quad element supported

            F_dot_n = 0.0 ! Init

            ! Extract data
            !Lint = ele_int_len(eleID,intID)
            Lint = U2D_GRID%CELL_EDGES_LEN(intID,eleID)
            !nx   = ele_int_nx(eleID,intID)
            nx = U2D_GRID%EDGE_NORMAL(1,intID,eleID)
            !ny   = ele_int_ny(eleID,intID)
            ny = U2D_GRID%EDGE_NORMAL(2,intID,eleID)

            !neigh = ele_neigh(eleID,intID)
            neigh = U2D_GRID%CELL_NEIGHBORS(intID,eleID)

            ! Check what neighbor is it
            if (neigh .NE. -1) then ! +++++++++ INTERNAL CELL

               call compute_fluxes_HLL(U(:,eleID), U(:,neigh), nx, ny, F_dot_n, Aele)

            else
               FACE_PG = U2D_GRID%CELL_EDGES_PG(intID,eleID)
               if (GRID_BC(FACE_PG)%PARTICLE_BC == INLET) then ! ++++++++ INLET BOUNDARY +++++++++++++++++++

                  call compute_fluxes_HLL(U(:,eleID), U_inlet, nx, ny, F_dot_n, Aele)

               else if (GRID_BC(FACE_PG)%PARTICLE_BC == OUTLET) then ! ++++++++ OUTLET BOUNDARY +++++++++++++++++++

                  call compute_fluxes_HLL(U(:,eleID), U_outlet, nx, ny, F_dot_n, Aele)

               else if (GRID_BC(FACE_PG)%PARTICLE_BC == WALL) then ! ++++++++ WALL BOUNDARY ++++++++++++++++++++

                  call compute_wall_flux(U(:,eleID), nx, ny, F_dot_n, Aele)

               else if (GRID_BC(FACE_PG)%PARTICLE_BC == SYMMETRY) then ! ++++++++ SYM BOUNDARY ++++++++++++++++++++

                  call compute_sym_state(U(:,eleID), nx, ny, U_sym)
                  call compute_fluxes_HLL(U(:,eleID), U_sym, nx, ny, F_dot_n, Aele)

               else
                  print*, "ERROR! UNKNOWN BOUNDARY TYPE ", neigh, " for element ", eleID, " Check the mesh or the pre-processing."
                  print*, "ABORTING!"
                  STOP
               end if
            end if

            ! Update solution
            U_new(:,eleID) = U_new(:,eleID) - dt*F_dot_n*Lint/Aele

            ! Check that the solution did not diverge
            do eqID = 1, Neq
               if (isnan(U_new(eqID,eleID))) then 
                  print*, 'Solution diverged, try with a smaller time step! Aborting.'
                  print*, 'Solution that diverged: ', U_new(:,eleID)
                  print*, 'in cell ID = ', eleID
                  stop
               end if
            end do

         end do ! End loop on interfaces
   
      end do ! End loop on elements
   
      U = U_new ! Save results

   end subroutine

   ! ======================================================================== 

   subroutine compute_fluxes_HLL(U_L, U_R, nx, ny, F_dot_n, A_ele)

      ! Computes HLL numerical fluxes among the cell eleID and the neighbor cell neigh
      ! The element area Aele is also passed, for the sake of computing the CFL number.

      implicit none

      real(kind=8), dimension(:),   intent(in)  :: U_L, U_R
      real(kind=8),                 intent(in)  :: nx, ny, A_ele
      real(kind=8), dimension(Neq), intent(out) :: F_dot_n 

      real(kind=8), dimension(Neq) :: F_L, F_R

      ! Wave speeds
      real(kind=8) :: ws_min_L, ws_max_L, ws_min_R, ws_max_R, ws_max, ws_min

      call compute_flux_ws(U_L, F_L, nx, ny, ws_max_L, ws_min_L)
      call compute_flux_ws(U_R, F_R, nx, ny, ws_max_R, ws_min_R)

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

      ! Update global maximum wave speed (used for setting the time step)
      ws_max = abs(ws_max)
      ws_over_sqrtA_maxabs = MAX(ws_over_sqrtA_maxabs, ws_max/sqrt(A_ele))

   end subroutine


   ! ======================================================================== 

   subroutine compute_fluxes_Rusanov(U_L, U_R, nx, ny, F_dot_n, A_ele)

      ! Computes Rusanov numerical fluxes among the cell eleID and the neighbor cell neigh
      ! The element area A_ele is also passed, for the sake of computing the Courant number

      implicit none

      real(kind=8), dimension(:),   intent(in)  :: U_L, U_R
      real(kind=8),                 intent(in)  :: nx, ny, A_ele
      real(kind=8), dimension(Neq), intent(out) :: F_dot_n 

      real(kind=8), dimension(Neq) :: F_L, F_R

      ! Wave speeds
      real(kind=8) :: ws_min_L, ws_max_L, ws_min_R, ws_max_R, ws_max

      call compute_flux_ws(U_L, F_L, nx, ny, ws_max_L, ws_min_L)
      call compute_flux_ws(U_R, F_R, nx, ny, ws_max_R, ws_min_R)

      ws_max  = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))

      F_dot_n = 0.5*(F_R + F_L) - ws_max/2.0*(U_R-U_L) ! Rusanov flux
   
      ! Update global maximum wave speed (used for setting the time step)
      ws_over_sqrtA_maxabs = MAX(ws_over_sqrtA_maxabs, ws_max/sqrt(A_ele))

   end subroutine

end module
