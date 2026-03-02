module tools

   use global_module
   use grid
   use pde 

   implicit none

   contains





   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE GRID_SAVE -> Saves cumulated average !!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   SUBROUTINE GRID_SAVE(t_ID, U)

      ! ---------------------------------------------------------------------------------------
      ! This subroutine exports data in the legacy VTK format, that can be visualized using
      ! ParaView.
      ! This helps for checking that all is good.
      ! The output file is put in the "dumps" directory and is named "sol_..." where ... is the
      ! time ID "t_ID" given as an input to the subroutine.
      ! U: conserved variables
      ! 
      ! Notes: 
      !  - All cells are printed, including the two layers of ghost cells for each boundary.
      !    This helps for debugging. 
      !  - The data is exported for each cell center. In ParaView, you may visualize the data 
      !    using the option "Surface with Edges": this is misleading, and these edges do not
      !    represent the cells boundaries. Instead, each point is a cell center.
      ! ---------------------------------------------------------------------------------------
      IMPLICIT NONE



      integer, intent(in) :: t_ID
      real(kind=8), dimension(:,:), intent(in) :: U

      real(kind=8), dimension(:,:), allocatable :: prim ! Primitive variables
      
      integer      :: eqID

      CHARACTER*256                      :: string, file_name
      CHARACTER*7, DIMENSION(33)         :: MOMENT_STRING 


      INTEGER                            :: I

      REAL(KIND=8) :: SPWT


      MOMENT_STRING = ['rho_   ', &
                        'Ux_    ', 'Uy_    ', 'Uz_    ', &
                        'Pxx_   ', 'Pxy_   ', 'Pxz_   ', 'Pyy_   ', 'Pyz_   ', 'Pzz_   ', &
                        'qx_    ', 'qy_    ', 'qz_    ', &
                        'Qxxx_  ', 'Qxxy_  ', 'Qxyy_  ', 'Qyyy_  ', 'Qyyz_  ', &
                        'Qyzz_  ', 'Qzzz_  ', 'Qxxz_  ', 'Qxzz_  ', 'Qxyz_  ', &
                        'Riijj_ ', 'Rxxjj_ ', 'Rxyjj_ ', 'Rxzjj_ ', 'Ryyjj_ ', 'Ryzjj_ ', 'Rzzjj_ ', &
                        'Sxiijj_', 'Syiijj_', 'Sziijj_']


      ! ----- Compute primitive variables on the grid ------

      allocate(prim(Neq,NCELLS))

      prim = 0.0 ! Init
      do i = 1, NCELLS
         call compute_primitive_from_conserved(U(:,i), prim(:,i))
      end do

      ! ------- Write VTK file -------

      WRITE(file_name,'(A, A, I0, A)') TRIM(ADJUSTL(FLOWFIELD_SAVE_PATH)), 'dsmc_flowfield_', t_ID, '.vtk'

      IF (BOOL_BINARY_OUTPUT) THEN
         OPEN(54321, FILE=file_name, ACCESS='STREAM', FORM='UNFORMATTED', STATUS='NEW', CONVERT='BIG_ENDIAN')

         WRITE(54321) '# vtk DataFile Version 3.0'//ACHAR(10)
         WRITE(54321) 'vtk output'//ACHAR(10)
         WRITE(54321) 'BINARY'//ACHAR(10)


         WRITE(54321) 'DATASET UNSTRUCTURED_GRID'//ACHAR(10)

         WRITE(54321) 'POINTS '//ITOA(U2D_GRID%NUM_NODES)//' double'//ACHAR(10)
         DO I = 1, U2D_GRID%NUM_NODES
            WRITE(54321) U2D_GRID%NODE_COORDS(:,I)
         END DO

         WRITE(54321) 'CELLS '//ITOA(U2D_GRID%NUM_CELLS)//' '//ITOA(4*U2D_GRID%NUM_CELLS)//ACHAR(10)
         DO I = 1, U2D_GRID%NUM_CELLS
            WRITE(54321) 3, (U2D_GRID%CELL_NODES(:,I) - 1)
         END DO

         WRITE(54321) 'CELL_TYPES '//ITOA(U2D_GRID%NUM_CELLS)//ACHAR(10)
         DO I = 1, U2D_GRID%NUM_CELLS
            WRITE(54321) 5
         END DO
        

         WRITE(54321) 'CELL_DATA '//ITOA(NCELLS)//ACHAR(10)
         IF (BOOL_DUMP_MOMENTS) THEN
            WRITE(54321) 'FIELD FieldData '//ITOA( 4 )//ACHAR(10)
         ELSE
            WRITE(54321) 'FIELD FieldData '//ITOA( 4 )//ACHAR(10)
         END IF


         ! Write per-cell value
         DO eqID = 1, Neq

            WRITE(54321) prim_names(eqID)//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
            WRITE(54321) prim(eqID,:), ACHAR(10)

         END DO

         CLOSE(54321)
         
      ELSE  ! Write ASCII output.
         OPEN(54321, FILE=file_name, ACCESS='SEQUENTIAL', FORM='FORMATTED', STATUS='NEW')

         WRITE(54321,'(A)') '# vtk DataFile Version 3.0'
         WRITE(54321,'(A)') 'vtk output'
         WRITE(54321,'(A)') 'ASCII'


         WRITE(54321,'(A)') 'DATASET UNSTRUCTURED_GRID'
         
         WRITE(54321,'(A,I10,A7)') 'POINTS', U2D_GRID%NUM_NODES, 'double'
         DO I = 1, U2D_GRID%NUM_NODES
            WRITE(54321,*) U2D_GRID%NODE_COORDS(:,I)
         END DO

         WRITE(54321,'(A,I10,I10)') 'CELLS', U2D_GRID%NUM_CELLS, 4*U2D_GRID%NUM_CELLS 
         DO I = 1, U2D_GRID%NUM_CELLS
            WRITE(54321,*) 3, (U2D_GRID%CELL_NODES(:,I) - 1)
         END DO

         WRITE(54321,'(A,I10)') 'CELL_TYPES', U2D_GRID%NUM_CELLS
         DO I = 1, U2D_GRID%NUM_CELLS
            WRITE(54321,*) 5
         END DO
 
         
         WRITE(54321,'(A,I10)') 'CELL_DATA', NCELLS
         IF (BOOL_DUMP_MOMENTS) THEN
            WRITE(54321,'(A,I10)') 'FIELD FieldData', 4
         ELSE
            WRITE(54321,'(A,I10)') 'FIELD FieldData', 4
         END IF


         ! Write per-cell value
         DO eqID = 1, Neq
            
            WRITE(54321,'(A,I10,I10,A8)') prim_names(eqID), 1, NCELLS, 'integer'
            WRITE(54321) prim(eqID,:)
            
         END DO

         CLOSE(54321)
      
      END IF


   END SUBROUTINE GRID_SAVE


   FUNCTION ITOA(I) RESULT(RES)
      CHARACTER(:), ALLOCATABLE :: RES
      INTEGER, INTENT(IN) :: I
      CHARACTER(RANGE(I)+2) :: TMP
      WRITE(TMP,'(I0)') I
      RES = TRIM(TMP)
   END FUNCTION


end module 
