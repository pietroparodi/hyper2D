module tools

   ! Most of the code in this module comes from the Pantera PIC-DSMC code,
   ! (https://github.com/vonkarmaninstitute/pantera-pic-dsmc) for reasons of compatibility.
   ! Pantera is a free software by the von Karman Institute for Fluid Dynamics (VKI), 
   ! distributed under a GNU GPLv3 license.

   use global_module
   use grid
   use pde

   implicit none

   contains

   SUBROUTINE GRID_SAVE(t_ID, TIME, U)

      ! This subroutine export data in the VTK (legacy) format.
      ! Data can be saved both in ASCII and in binary format (preferrable).
      ! The files can be read directly in Paraview.
      ! A field named "TIME" contains the physical time of the snapshot.
      ! This is not directly used by Paraview but can be used in post-processing. 

      IMPLICIT NONE

      integer, intent(in) :: t_ID
      real(kind=8), intent(in) :: TIME
      real(kind=8), dimension(:,:), intent(in) :: U

      real(kind=8), dimension(:,:), allocatable :: prim, Uprim ! Primitive variables
      real(kind=8), dimension(:,:,:), allocatable :: gradUprim
   
      integer      :: eqID
      CHARACTER*256                      :: file_name
      INTEGER                            :: I, J, IC, FIRST, LAST

      character(len=20), dimension(:), allocatable :: prim_names, transport_names, grad_names

      ALLOCATE(prim_names(N_SPECIES_FLUID*5))
      ALLOCATE(transport_names(N_SPECIES_FLUID*3))
      ALLOCATE(grad_names(N_SPECIES_FLUID*4))


      DO I = 1, N_SPECIES_FLUID
         J = 5*(I-1)
         prim_names(J+1) = 'nrho_mean_'//TRIM(SPECIES(I)%NAME)
         prim_names(J+2) = 'rho_mean_'//TRIM(SPECIES(I)%NAME)
         prim_names(J+3) = 'vx_mean_'//TRIM(SPECIES(I)%NAME)
         prim_names(J+4) = 'vy_mean_'//TRIM(SPECIES(I)%NAME)
         prim_names(J+5) = 'Ttr_mean_'//TRIM(SPECIES(I)%NAME)
      END DO

      DO I = 1, N_SPECIES_FLUID
         J = 4*(I-1)
         grad_names(J+1) = 'grad_rho_'//TRIM(SPECIES(I)%NAME)
         grad_names(J+2) = 'grad_vx_'//TRIM(SPECIES(I)%NAME)
         grad_names(J+3) = 'grad_vy_'//TRIM(SPECIES(I)%NAME)
         grad_names(J+4) = 'grad_T_'//TRIM(SPECIES(I)%NAME)
      END DO

      DO I = 1, N_SPECIES_FLUID
         J = 3*(I-1)
         transport_names(J+1) = 'mu_'//TRIM(SPECIES(I)%NAME)
         transport_names(J+2) = 'kappa_'//TRIM(SPECIES(I)%NAME)//'I'
         transport_names(J+3) = 'kappa_'//TRIM(SPECIES(I)%NAME)//'I2'
      END DO

      ! ----- Compute primitive variables on the grid ------

      allocate(prim(N_SPECIES_FLUID*5,NCELLS))
      allocate(Uprim(N_SPECIES_FLUID*4,NCELLS))
      allocate(gradUprim(2,N_SPECIES_FLUID*4,NCELLS))

      prim = 0.0 ! Init
      do IC = 1, NCELLS
         ! Skip cells that are not fluid
         IF (U2D_GRID%CELL_PG(IC) .NE. -1) THEN
            IF (GRID_BC(U2D_GRID%CELL_PG(IC))%VOLUME_BC == SOLID) CYCLE
         END IF
         
         DO I = 1, N_SPECIES_FLUID
            call compute_primitive_from_conserved(U((I-1)*Neq+1:(I-1)*Neq+4,IC), prim((I-1)*5+2:(I-1)*5+5,IC), I)
            prim((I-1)*5+1,IC) = prim((I-1)*5+2,IC) / SPECIES(I)%MOLECULAR_MASS
         END DO
      end do

      DO IC = 1, NCELLS
         ! Skip cells that are not fluid
         IF (U2D_GRID%CELL_PG(IC) .NE. -1) THEN
            IF (GRID_BC(U2D_GRID%CELL_PG(IC))%VOLUME_BC == SOLID) CYCLE
         END IF

         DO I = 1, N_SPECIES_FLUID
            FIRST = (I-1)*Neq+1
            LAST = I*Neq
            call compute_primitive_from_conserved(U(FIRST:LAST,IC), Uprim(FIRST:LAST,IC), I)
         END DO
      END DO

      CALL compute_cell_centered_gradients_weighted_least_squares(Uprim, gradUprim)
      !CALL compute_cell_centered_gradients_green_gauss(Uprim, gradUprim)

      ! ------- Write VTK file -------

      WRITE(file_name,'(A, A, I0, A)') TRIM(ADJUSTL(FLOWFIELD_SAVE_PATH)), 'fluid_flowfield_', t_ID, '.vtk'

      IF (BOOL_BINARY_OUTPUT) THEN
         OPEN(54321, FILE=file_name, ACCESS='STREAM', FORM='UNFORMATTED', STATUS='NEW', CONVERT='BIG_ENDIAN')

         WRITE(54321) '# vtk DataFile Version 3.0'//ACHAR(10)
         WRITE(54321) 'vtk output'//ACHAR(10)
         WRITE(54321) 'BINARY'//ACHAR(10)


         WRITE(54321) 'DATASET UNSTRUCTURED_GRID'//ACHAR(10)

         WRITE(54321) 'FIELD FieldData 1'//ACHAR(10)
         WRITE(54321) 'TIME 1 1 double'//ACHAR(10)
         WRITE(54321) TIME


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

         WRITE(54321) 'FIELD FieldData '//ITOA( N_SPECIES_FLUID*(5+2*4) + 6 )//ACHAR(10)


         ! Write per-cell value
         DO eqID = 1, N_SPECIES_FLUID*5

            WRITE(54321) prim_names(eqID)//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
            WRITE(54321) prim(eqID,:), ACHAR(10)

         END DO

         ! Write per-cell value
         DO eqID = 1, N_SPECIES_FLUID*4

            WRITE(54321) TRIM(grad_names(eqID))//'_x '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
            WRITE(54321) gradUprim(1,eqID,:), ACHAR(10)
            
            WRITE(54321) TRIM(grad_names(eqID))//'_y '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
            WRITE(54321) gradUprim(2,eqID,:), ACHAR(10)

         END DO


         WRITE(54321) 'mu_I2 '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
         WRITE(54321) MU_GRID(1,:), ACHAR(10)
         WRITE(54321) 'mu_I '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
         WRITE(54321) MU_GRID(2,:), ACHAR(10)
         
         WRITE(54321) 'kappa_I2_I2 '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
         WRITE(54321) KAPPA_GRID(1,1,:), ACHAR(10)

         WRITE(54321) 'kappa_I2_I '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
         WRITE(54321) KAPPA_GRID(2,1,:), ACHAR(10)

         WRITE(54321) 'kappa_I_I2 '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
         WRITE(54321) KAPPA_GRID(1,2,:), ACHAR(10)

         WRITE(54321) 'kappa_I_I '//ITOA(1)//' '//ITOA(NCELLS)//' double'//ACHAR(10)
         WRITE(54321) KAPPA_GRID(2,2,:), ACHAR(10)

         CLOSE(54321)
         
      ELSE  ! Write ASCII output.
         OPEN(54321, FILE=file_name, ACCESS='SEQUENTIAL', FORM='FORMATTED', STATUS='NEW')

         WRITE(54321,'(A)') '# vtk DataFile Version 3.0'
         WRITE(54321,'(A)') 'vtk output'
         WRITE(54321,'(A)') 'ASCII'


         WRITE(54321,'(A)') 'DATASET UNSTRUCTURED_GRID'

         WRITE(54321,'(A)') 'FIELD FieldData 1'
         WRITE(54321,'(A)') 'TIME 1 1 double'
         WRITE(54321,*) TIME
         
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
         WRITE(54321,'(A,I10)') 'FIELD FieldData', N_SPECIES_FLUID*5



         ! Write per-cell value
         DO eqID = 1, N_SPECIES_FLUID*5
            
            WRITE(54321,'(A,I10,I10,A8)') prim_names(eqID), 1, NCELLS, 'integer'
            WRITE(54321) prim(eqID,:)
            
         END DO

         CLOSE(54321)
      
      END IF

      DEALLOCATE(prim_names)
      DEALLOCATE(grad_names)
      DEALLOCATE(prim)
      DEALLOCATE(Uprim)
      DEALLOCATE(gradUprim)


   END SUBROUTINE GRID_SAVE


   FUNCTION ITOA(I) RESULT(RES)

      ! Function to turn an integer into a string
      ! Useful when writing files in binary mode.

      CHARACTER(:), ALLOCATABLE :: RES
      INTEGER, INTENT(IN) :: I
      CHARACTER(RANGE(I)+2) :: TMP
      WRITE(TMP,'(I0)') I
      RES = TRIM(TMP)
   END FUNCTION




   INTEGER FUNCTION SPECIES_NAME_TO_ID(NAME)

      ! Function that returns the ID of a species given its name string,
      ! as defined in the species definition file.
      ! This is primarily used to parse the input file.

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)  :: NAME
      INTEGER                       :: INDEX, MATCH
      MATCH = -1
      DO INDEX = 1, N_SPECIES
         IF (SPECIES(INDEX)%NAME == NAME) MATCH = INDEX
      END DO


      SPECIES_NAME_TO_ID = MATCH

   END FUNCTION SPECIES_NAME_TO_ID


   INTEGER FUNCTION MIXTURE_NAME_TO_ID(NAME)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)  :: NAME
      INTEGER                       :: INDEX, MATCH
      MATCH = -1
      DO INDEX = 1, N_MIXTURES
         IF (MIXTURES(INDEX)%NAME == NAME) MATCH = INDEX
      END DO

      IF (MATCH .EQ. -1) THEN
         WRITE(*,*) 'Error! Mixture name not found.'
      END IF

      MIXTURE_NAME_TO_ID = MATCH

   END FUNCTION MIXTURE_NAME_TO_ID



   SUBROUTINE STRIP_COMMENTS(str,c)

      ! This subroutine checks if there is any character of type 'c' in the string str,
      ! and if so, it removes all the remaining of the string, tagged as a comment.
      ! Also, white spaces are removed by a call to the trim() function.

      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(INOUT) :: str
      CHARACTER(LEN=1),INTENT(IN)    :: c !comment character

      CHARACTER(LEN=LEN(str)) :: str_tmp
      INTEGER :: i
      
      ! Check if there is any comment to trim
      i = INDEX(str,c)
      IF (i .GT. 0) THEN
         str_tmp = str(1:i-1)
      ELSE
         str_tmp = str
      END IF
      
      ! Assign str, removing trailing blank spaces if any
      str = TRIM(str_tmp)

   END SUBROUTINE STRIP_COMMENTS



   SUBROUTINE SKIP_TO(UNIT, STR, STAT)

      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: STR
      INTEGER, INTENT(IN) :: UNIT
      INTEGER, INTENT(OUT) :: STAT
      CHARACTER :: CH
      INTEGER :: IO
    
      DO
         READ(UNIT, IOSTAT=IO) CH

         IF (IO/=0) THEN
            STAT = 1
            RETURN
         END IF
    
         IF (CH==STR(1:1)) THEN
            IF (LEN(STR) == 1) THEN
               STAT = 0
               RETURN
            END IF
            CALL CHECK(UNIT, STR(2:), STAT)
            IF (STAT == 0) RETURN
         END IF
    
      END DO
   END SUBROUTINE

    
   SUBROUTINE CHECK(UNIT, STR, STAT)
      CHARACTER(*), INTENT(IN) :: STR
      INTEGER, INTENT(IN) :: UNIT
      INTEGER, INTENT(OUT) :: STAT
      CHARACTER :: CH
      INTEGER :: I, IO

      STAT = 1
      I = 0

      DO
         I = I + 1

         READ(UNIT, IOSTAT=IO) CH

         IF (IO/=0) RETURN

         IF (CH/=STR(I:I)) RETURN

         IF (I==LEN(STR)) THEN
            STAT = 0
            RETURN
         END IF
      END DO
   END SUBROUTINE CHECK



   FUNCTION INTERP_RATE(VALUE_TEMP, TABLE_TEMP, TABLE_RATE) RESULT(VALUE_RATE)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: VALUE_TEMP
      INTEGER :: L, R
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: TABLE_TEMP, TABLE_RATE
      INTEGER :: INDEX
      REAL(KIND=8) :: VALUE_RATE

      L = LBOUND(TABLE_TEMP, DIM=1)
      R = UBOUND(TABLE_TEMP, DIM=1)

      INDEX = -1
      IF (VALUE_TEMP .LT. TABLE_TEMP(L)) THEN
         ! Lower than lower energy value
         !VALUE_RATE = TABLE_RATE(L)
         VALUE_RATE = 0
         RETURN
      ELSE IF (VALUE_TEMP .GT. TABLE_TEMP(R)) THEN
         ! Higher than highest energy value
         !VALUE_RATE = TABLE_RATE(R)
         VALUE_RATE = 0
         RETURN
      ELSE IF (R == L+1) THEN
         ! Only two values in the table
         VALUE_RATE = TABLE_RATE(L) + (TABLE_RATE(R)-TABLE_RATE(L))*(VALUE_TEMP-TABLE_TEMP(L))/(TABLE_TEMP(R)-TABLE_TEMP(L))
         RETURN
      ELSE
         DO
            INDEX = (L+R)/2
            IF (TABLE_TEMP(INDEX) .LE. VALUE_TEMP) THEN
               IF (TABLE_TEMP(INDEX+1) .GT. VALUE_TEMP) EXIT
               L = INDEX
            ELSE
               IF (TABLE_TEMP(INDEX-1) .LE. VALUE_TEMP) THEN
                  INDEX = INDEX-1
                  EXIT
               END IF
               R = INDEX
            END IF
         END DO
         ! The value we are looking for is between INDEX and INDEX+1.
         L = INDEX
         R = INDEX+1
         VALUE_RATE = TABLE_RATE(L) + (TABLE_RATE(R)-TABLE_RATE(L))*(VALUE_TEMP-TABLE_TEMP(L))/(TABLE_TEMP(R)-TABLE_TEMP(L))
         RETURN
      END IF

   END FUNCTION INTERP_RATE



   SUBROUTINE WRITE_RESTART(TIMESTEP, U)

      IMPLICIT NONE

      INTEGER, intent(in) :: TIMESTEP
      real(kind=8), dimension(:,:), intent(in) :: U
      INTEGER :: ios

      character(len=100) :: iomsg

      CHARACTER(LEN=512)  :: filename

      WRITE(filename, "(A,A8,I0.8)") TRIM(ADJUSTL(FLOWFIELD_SAVE_PATH)), "restart_", TIMESTEP ! Compose filename

      ! Open file for writing
      IF (BOOL_BINARY_OUTPUT) THEN
         OPEN(1010, FILE=filename, ACCESS='STREAM', FORM='UNFORMATTED', STATUS='NEW', &
         CONVERT='BIG_ENDIAN', IOSTAT=ios, IOMSG=iomsg)

         IF (ios .NE. 0) THEN
            CALL ERROR_ABORT('Attention, could not open restart file for writing! ABORTING.')
         ENDIF

         WRITE(1010, IOSTAT=ios) U
         !IF (ios < 0) EXIT

         CLOSE(1010)
      ELSE
         OPEN(1010, FILE=filename, STATUS='OLD', IOSTAT=ios)
         
         IF (ios .NE. 0) THEN
            CALL ERROR_ABORT('Attention, could not open restart file for writing! ABORTING.')
         ENDIF

         WRITE(1010,*,IOSTAT=ios) U
         !IF (ios < 0) EXIT

         CLOSE(1010)
      END IF

   END SUBROUTINE



   SUBROUTINE CHECKS(TIME, U_CONS)

      IMPLICIT NONE

      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), INTENT(IN) :: U_CONS
      REAL(KIND=8), INTENT(IN) :: TIME

      INTEGER                            :: JS, JC, FIRST
   
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: TOT_MASS, TOT_ENERGY
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: TOT_MOMENTUM
   
      REAL(KIND=8) :: VOL
      CHARACTER*256                      :: file_name
      CHARACTER*2048                     :: HEADER_STRING
      LOGICAL                            :: FILE_EXISTS

      ALLOCATE(TOT_MASS(N_SPECIES_FLUID))
      ALLOCATE(TOT_ENERGY(N_SPECIES_FLUID))
      ALLOCATE(TOT_MOMENTUM(3,N_SPECIES_FLUID))

      TOT_MASS = 0.d0
      TOT_MOMENTUM = 0.d0
      TOT_ENERGY = 0.d0

      DO JC = 1, NCELLS
         DO JS = 1, N_SPECIES_FLUID
            FIRST = (JS-1)*Neq+1
            VOL = U2D_GRID%CELL_VOLUMES(JC)
            TOT_MASS(JS) = TOT_MASS(JS) + U_CONS(FIRST, JC)*VOL
            TOT_MOMENTUM(1, JS) = TOT_MOMENTUM(1, JS) + U_CONS(FIRST+1, JC)*VOL
            TOT_MOMENTUM(2, JS) = TOT_MOMENTUM(2, JS) + U_CONS(FIRST+2, JC)*VOL
            TOT_MOMENTUM(3, JS) = 0.d0
            TOT_ENERGY(JS) = TOT_ENERGY(JS) + U_CONS(FIRST+3, JC)*VOL
         END DO
      END DO

      WRITE(file_name,'(A, A)') TRIM(ADJUSTL(FLOWFIELD_SAVE_PATH)), 'conservation_checks'

      INQUIRE(FILE=file_name, EXIST=FILE_EXISTS)

      OPEN(54331, FILE=file_name, POSITION='append', STATUS='unknown', ACTION='write')
      IF (.NOT. FILE_EXISTS) THEN
         HEADER_STRING = ''
         HEADER_STRING = TRIM(HEADER_STRING) // 'time'
         DO JS = 1, N_SPECIES
            HEADER_STRING = TRIM(HEADER_STRING) // ' mass_' // TRIM(SPECIES(JS)%NAME)
         END DO
         HEADER_STRING = TRIM(HEADER_STRING) // ' totmass'
         DO JS = 1, N_SPECIES
            HEADER_STRING = TRIM(HEADER_STRING) // ' xmom_' // TRIM(SPECIES(JS)%NAME) // ' ' &
                                          // 'ymom_' // TRIM(SPECIES(JS)%NAME) // ' ' &
                                          // 'zmom_' // TRIM(SPECIES(JS)%NAME)
         END DO
         HEADER_STRING = TRIM(HEADER_STRING) // ' totxmom totymom totzmom'
         DO JS = 1, N_SPECIES
            HEADER_STRING = TRIM(HEADER_STRING) // ' ene_' // TRIM(SPECIES(JS)%NAME)
         END DO
         HEADER_STRING = TRIM(HEADER_STRING) // ' totene'

         WRITE(54331,*) TRIM(HEADER_STRING)
      END IF

      WRITE(54331,*) TIME, TOT_MASS, SUM(TOT_MASS), TOT_MOMENTUM, SUM(TOT_MOMENTUM, DIM=2), &
      TOT_ENERGY, SUM(TOT_ENERGY)
      CLOSE(54331)


   END SUBROUTINE CHECKS



   FUNCTION MU_I(T, X_I) RESULT(MU)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: T, X_I
      REAL(KIND=8) :: MU

      IF (.NOT. TABLE_MU_I%INITIALIZED) THEN
         CALL read_2d_table('mu_I_vss.txt', TABLE_MU_I%NROWS, TABLE_MU_I%NCOLS, &
         TABLE_MU_I%ROW_IDX, TABLE_MU_I%COL_IDX, TABLE_MU_I%DATA)
         TABLE_MU_I%INITIALIZED = .TRUE.
      END IF
      MU = interpolate_2d_table(TABLE_MU_I%ROW_IDX, TABLE_MU_I%COL_IDX, &
      TABLE_MU_I%DATA, T, X_I)

   END FUNCTION

   FUNCTION MU_I2(T, X_I) RESULT(MU)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: T, X_I
      REAL(KIND=8) :: MU

      IF (.NOT. TABLE_MU_I2%INITIALIZED) THEN
         CALL read_2d_table('mu_I2_vss.txt', TABLE_MU_I2%NROWS, TABLE_MU_I2%NCOLS, &
         TABLE_MU_I2%ROW_IDX, TABLE_MU_I2%COL_IDX, TABLE_MU_I2%DATA)
         TABLE_MU_I2%INITIALIZED = .TRUE.
      END IF
      MU = interpolate_2d_table(TABLE_MU_I2%ROW_IDX, TABLE_MU_I2%COL_IDX, &
      TABLE_MU_I2%DATA, T, X_I)

   END FUNCTION

   FUNCTION KAPPA_I_I(T, X_I) RESULT(KAPPA)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: T, X_I
      REAL(KIND=8) :: KAPPA

      IF (.NOT. TABLE_KAPPA_I_I%INITIALIZED) THEN
         CALL read_2d_table('kappa_I_I_vss.txt', TABLE_KAPPA_I_I%NROWS, TABLE_KAPPA_I_I%NCOLS, &
         TABLE_KAPPA_I_I%ROW_IDX, TABLE_KAPPA_I_I%COL_IDX, TABLE_KAPPA_I_I%DATA)
         TABLE_KAPPA_I_I%INITIALIZED = .TRUE.
      END IF
      KAPPA = interpolate_2d_table(TABLE_KAPPA_I_I%ROW_IDX, TABLE_KAPPA_I_I%COL_IDX, &
      TABLE_KAPPA_I_I%DATA, T, X_I)

   END FUNCTION

   FUNCTION KAPPA_I_I2(T, X_I) RESULT(KAPPA)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: T, X_I
      REAL(KIND=8) :: KAPPA

      IF (.NOT. TABLE_KAPPA_I_I2%INITIALIZED) THEN
         CALL read_2d_table('kappa_I_I2_vss.txt', TABLE_KAPPA_I_I2%NROWS, TABLE_KAPPA_I_I2%NCOLS, &
         TABLE_KAPPA_I_I2%ROW_IDX, TABLE_KAPPA_I_I2%COL_IDX, TABLE_KAPPA_I_I2%DATA)
         TABLE_KAPPA_I_I2%INITIALIZED = .TRUE.
      END IF
      KAPPA = interpolate_2d_table(TABLE_KAPPA_I_I2%ROW_IDX, TABLE_KAPPA_I_I2%COL_IDX, &
      TABLE_KAPPA_I_I2%DATA, T, X_I)

   END FUNCTION

   FUNCTION KAPPA_I2_I(T, X_I) RESULT(KAPPA)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: T, X_I
      REAL(KIND=8) :: KAPPA

      IF (.NOT. TABLE_KAPPA_I2_I%INITIALIZED) THEN
         CALL read_2d_table('kappa_I2_I_vss.txt', TABLE_KAPPA_I2_I%NROWS, TABLE_KAPPA_I2_I%NCOLS, &
         TABLE_KAPPA_I2_I%ROW_IDX, TABLE_KAPPA_I2_I%COL_IDX, TABLE_KAPPA_I2_I%DATA)
         TABLE_KAPPA_I2_I%INITIALIZED = .TRUE.
      END IF
      KAPPA = interpolate_2d_table(TABLE_KAPPA_I2_I%ROW_IDX, TABLE_KAPPA_I2_I%COL_IDX, &
      TABLE_KAPPA_I2_I%DATA, T, X_I)

   END FUNCTION

   FUNCTION KAPPA_I2_I2(T, X_I) RESULT(KAPPA)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: T, X_I
      REAL(KIND=8) :: KAPPA

      IF (.NOT. TABLE_KAPPA_I2_I2%INITIALIZED) THEN
         CALL read_2d_table('kappa_I2_I2_vss.txt', TABLE_KAPPA_I2_I2%NROWS, TABLE_KAPPA_I2_I2%NCOLS, &
         TABLE_KAPPA_I2_I2%ROW_IDX, TABLE_KAPPA_I2_I2%COL_IDX, TABLE_KAPPA_I2_I2%DATA)
         TABLE_KAPPA_I2_I2%INITIALIZED = .TRUE.
      END IF
      KAPPA = interpolate_2d_table(TABLE_KAPPA_I2_I2%ROW_IDX, TABLE_KAPPA_I2_I2%COL_IDX, &
      TABLE_KAPPA_I2_I2%DATA, T, X_I)

   END FUNCTION


   function interpolate_2d_table(row_indices, col_indices, table_data, row_val, col_val) result(interp_val)
      implicit none

      real(KIND=8), intent(in) :: row_indices(:), col_indices(:), table_data(:,:)
      real(KIND=8), intent(in) :: row_val, col_val
      integer :: ierr
      real(KIND=8) :: interp_val

      integer :: nrows, ncols, i, j
      real(KIND=8) :: x1, x2, y1, y2, f11, f12, f21, f22, t, u

      ierr = 0
      nrows = size(row_indices)
      ncols = size(col_indices)

      ! Check if row_val and col_val are within the bounds of the table
      if (row_val < row_indices(1) .or. row_val > row_indices(nrows) .or. &
         col_val < col_indices(1) .or. col_val > col_indices(ncols)) then
         ierr = 1
         print *, "Error: Target row or column value is out of bounds for interpolation."
         interp_val = 0.0
         return
      end if

      ! Find the interval for row_val
      do i = 1, nrows - 1
         if (row_val >= row_indices(i) .and. row_val <= row_indices(i+1)) exit
      end do

      ! Find the interval for col_val
      do j = 1, ncols - 1
         if (col_val >= col_indices(j) .and. col_val <= col_indices(j+1)) exit
      end do

      ! Set the four surrounding points
      x1 = row_indices(i)
      x2 = row_indices(i+1)
      y1 = col_indices(j)
      y2 = col_indices(j+1)

      f11 = table_data(i, j)
      f12 = table_data(i, j+1)
      f21 = table_data(i+1, j)
      f22 = table_data(i+1, j+1)

      ! Calculate the interpolation weights
      t = (row_val - x1) / (x2 - x1)
      u = (col_val - y1) / (y2 - y1)

      ! Perform bilinear interpolation
      interp_val = (1 - t) * (1 - u) * f11 + t * (1 - u) * f21 + &
                  (1 - t) * u * f12 + t * u * f22

   end function interpolate_2d_table


   subroutine read_2d_table(filename, nrows, ncols, row_indices, col_indices, table_data)
      implicit none

      character(len=*), intent(in) :: filename
      integer, intent(out) :: nrows, ncols
      real(KIND=8), allocatable, intent(out) :: row_indices(:), col_indices(:), table_data(:,:)

      integer :: i, j, io_status, ierr

      ! Initialize error flag
      ierr = 0

      ! Open the file
      open(unit=10, file=filename, status='old', iostat=io_status)
      if (io_status /= 0) then
         ierr = 1
         print *, "Error: Could not open file ", filename
         return
      end if

      ! Read number of rows and columns
      read(10, *, iostat=io_status) nrows
      read(10, *, iostat=io_status) ncols
      if (io_status /= 0) then
         ierr = 2
         print *, "Error: Could not read nrows and ncols"
         close(10)
         return
      end if

      ! Allocate arrays for row and column indices
      allocate(row_indices(nrows), stat=io_status)
      if (io_status /= 0) then
         ierr = 3
         print *, "Error: Could not allocate row_indices"
         close(10)
         return
      end if

      allocate(col_indices(ncols), stat=io_status)
      if (io_status /= 0) then
         ierr = 4
         print *, "Error: Could not allocate col_indices"
         deallocate(row_indices)
         close(10)
         return
      end if

      ! Allocate array for table data
      allocate(table_data(nrows, ncols), stat=io_status)
      if (io_status /= 0) then
         ierr = 5
         print *, "Error: Could not allocate table_data"
         deallocate(row_indices, col_indices)
         close(10)
         return
      end if

      ! Read row indices
      read(10, *, iostat=io_status) (row_indices(i), i=1, nrows)
      if (io_status /= 0) then
         ierr = 6
         print *, "Error: Could not read row indices"
         deallocate(row_indices, col_indices, table_data)
         close(10)
         return
      end if

      ! Read column indices
      read(10, *, iostat=io_status) (col_indices(j), j=1, ncols)
      if (io_status /= 0) then
         ierr = 7
         print *, "Error: Could not read column indices"
         deallocate(row_indices, col_indices, table_data)
         close(10)
         return
      end if

      ! Read table data in 1.8e format
      do i = 1, nrows
         read(10, *, iostat=io_status) (table_data(i, j), j=1, ncols)
         if (io_status /= 0) then
               ierr = 8
               print *, "Error: Could not read table data for row ", i
               deallocate(row_indices, col_indices, table_data)
               close(10)
               return
         end if
      end do

      ! Close the file
      close(10)

   end subroutine read_2d_table


end module 
