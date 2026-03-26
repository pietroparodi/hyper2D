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

      real(kind=8), dimension(:,:), allocatable :: prim ! Primitive variables
   
      integer      :: eqID
      CHARACTER*256                      :: file_name
      INTEGER                            :: I, J, IC

      character(len=20), dimension(:), allocatable :: prim_names

      ALLOCATE(prim_names(N_SPECIES_FLUID*Neq))


      DO I = 1, N_SPECIES_FLUID
         J = 4*(I-1)
         prim_names(J+1) = 'rho_'//TRIM(SPECIES(I)%NAME)
         prim_names(J+2) = 'ux_'//TRIM(SPECIES(I)%NAME)
         prim_names(J+3) = 'uy_'//TRIM(SPECIES(I)%NAME)
         prim_names(J+4) = 'T_'//TRIM(SPECIES(I)%NAME)
      END DO

      ! ----- Compute primitive variables on the grid ------

      allocate(prim(N_SPECIES_FLUID*Neq,NCELLS))

      prim = 0.0 ! Init
      do IC = 1, NCELLS
         DO I = 1, N_SPECIES_FLUID
            call compute_primitive_from_conserved(U((I-1)*Neq+1:I*Neq+1,IC), prim((I-1)*Neq+1:I*Neq+1,IC), I)
         END DO
      end do

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

         WRITE(54321) 'FIELD FieldData '//ITOA( N_SPECIES_FLUID*Neq )//ACHAR(10)


         ! Write per-cell value
         DO eqID = 1, N_SPECIES_FLUID*Neq

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
         WRITE(54321,'(A,I10)') 'FIELD FieldData', N_SPECIES_FLUID*Neq



         ! Write per-cell value
         DO eqID = 1, N_SPECIES_FLUID*Neq
            
            WRITE(54321,'(A,I10,I10,A8)') prim_names(eqID), 1, NCELLS, 'integer'
            WRITE(54321) prim(eqID,:)
            
         END DO

         CLOSE(54321)
      
      END IF

      DEALLOCATE(prim_names)

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


end module 
