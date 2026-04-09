MODULE initialization

   ! Most of the code in this module comes from the Pantera PIC-DSMC code,
   ! (https://github.com/vonkarmaninstitute/pantera-pic-dsmc) for reasons of compatibility.
   ! Pantera is a free software by the von Karman Institute for Fluid Dynamics (VKI), 
   ! distributed under a GNU GPLv3 license.

   USE global_module
   USE grid
   USE tools

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE READINPUT

      IMPLICIT NONE

      INTEGER, PARAMETER :: in1 = 444
      INTEGER            :: ios
      CHARACTER*512      :: line
      INTEGER            :: ReasonEOF

      CHARACTER*512      :: MIXTURE_DEFINITION, VSS_PARAMS_FILENAME, BG_FILENAME
      CHARACTER*512      :: BC_DEFINITION, SPECIES_FILENAME
      CHARACTER*64       :: BG_MIX_NAME, FLUX_FUNCTION_STRING

      ! Open input file for reading
      OPEN(UNIT=in1,FILE='input.fluid', STATUS='old',IOSTAT=ios)

      IF (ios.NE.0) THEN
         CALL ERROR_ABORT('Attention, "input" file not found! ABORTING.')
      ENDIF

      line = '' ! Init empty

      ! +++++++ Read until the end of file ++++++++
      DO

         READ(in1,'(A)', IOSTAT=ReasonEOF) line ! Read line
         CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

         IF (ReasonEOF < 0) EXIT ! End of file reached

         ! ~~~~~~~~~~~~~  Geometry and computational domain  ~~~~~~~~~~~~~~~~~

         IF (line=='Axisymmetric:')            READ(in1,*) AXI

         IF (line=='Wall_reactions_file:') THEN
            READ(in1,'(A)') WALL_REACTIONS_FILENAME
            CALL READ_WALL_REACTIONS(WALL_REACTIONS_FILENAME)
         END IF

         IF (line=='Mesh_file_SU2:') THEN
            READ(in1,'(A)') MESH_FILENAME
            CALL READ_2D_UNSTRUCTURED_GRID_SU2(MESH_FILENAME)
         END IF

         IF (line=='Boundary_condition:') THEN
            READ(in1,'(A)') BC_DEFINITION
            CALL DEF_BOUNDARY_CONDITION(BC_DEFINITION)
         END IF

         IF (line=='Domain_type:') THEN
            READ(in1,'(A)') BC_DEFINITION
            CALL DEF_DOMAIN_TYPE(BC_DEFINITION)
         END IF


         ! ~~~~~~~~~~~~~  Numerical settings  ~~~~~~~~~~~~~~~~~
         IF (line=='Timestep:')                READ(in1,*) dt_save
         IF (line=='Number_of_timesteps:')     READ(in1,*) Nt
         IF (line=='Max_timestep:')            READ(in1,*) dtmax
         IF (line=='Target_CFL:')              READ(in1,*) CFL_target
         IF (line=='Stats_every:')             READ(in1,*) STATS_EVERY
         IF (line=='Number_fluid_species:') THEN
            READ(in1,*) N_SPECIES_FLUID
            IF (N_SPECIES_FLUID > N_SPECIES) CALL ERROR_ABORT( &
            'Error! The requested number of fluid species is larger than the number of loaded species.')
         END IF
         
         ! ~~~~~~~~~~~~~  File output ~~~~~~~~~~~~~~~

         IF (line=='Output_path:') THEN
            READ(in1,'(A)') FLOWFIELD_SAVE_PATH
         END IF
         IF (line=='Binary_output:')           READ(in1,*) BOOL_BINARY_OUTPUT
         !IF (line=='Dump_grid_every:')         READ(in1,*) DUMP_GRID_EVERY
         !IF (line=='Dump_grid_start:')         READ(in1,*) DUMP_GRID_START

         IF (line=='Dump_state_every:')         READ(in1,*) DUMP_STATE_EVERY
         !IF (line=='Dump_state_start:')         READ(in1,*) DUMP_STATE_START


         ! ~~~~~~~~~~~~~  Multifluid ~~~~~~~~~~~~~~~
         IF (line=='Species_file:') THEN
            READ(in1,'(A)') SPECIES_FILENAME
            CALL READ_SPECIES(SPECIES_FILENAME)
         END IF

         IF (line=='Def_mixture:') THEN
            READ(in1,'(A)') MIXTURE_DEFINITION
            CALL DEF_MIXTURE(MIXTURE_DEFINITION)
         END IF

         IF (line=='Background_dens:')     READ(in1,*) BG_DENS
         IF (line=='Background_temp:')     READ(in1,*) BG_TEMP
         IF (line=='Background_mixture:') THEN
            READ(in1,*) BG_MIX_NAME
            BG_MIX = MIXTURE_NAME_TO_ID(BG_MIX_NAME)
         END IF
         IF (line=='Background_file:') THEN
            READ(in1,'(A)') BG_FILENAME
            CALL READ_BACKGROUND_FILE(BG_FILENAME)
         END IF

         IF (line=='VSS_parameters_file:')     THEN
            READ(in1,'(A)') VSS_PARAMS_FILENAME
            CALL READ_VSS(VSS_PARAMS_FILENAME)
         END IF

         IF (line=='VSS_parameters_binary_file:')     THEN
            READ(in1,'(A)') VSS_PARAMS_FILENAME
            CALL READ_VSS_BINARY(VSS_PARAMS_FILENAME)
         END IF

        ! ~~~~~~~~~~~~~  Reactions ~~~~~~~~~~~~~~~
         IF (line=='Reactions_file:') THEN
            READ(in1,'(A)') REACTIONS_FILENAME
            CALL READ_REACTIONS(REACTIONS_FILENAME)
         END IF

         ! ~~~~~~~~~~~~~  Initial particles seeding  ~~~~~~~~~~~~~~~~~
         IF (line=='Initial_state:') THEN
            READ(in1,'(A)') BC_DEFINITION
            CALL DEF_INITIAL_STATE(BC_DEFINITION)
         END IF


         IF (line=='Flux_function:') THEN
            READ(in1,*) FLUX_FUNCTION_STRING
            IF (FLUX_FUNCTION_STRING == "AUSM") THEN
               FLUX_FUNCTION = AUSM
            ELSE IF (FLUX_FUNCTION_STRING == "HLL") THEN
               FLUX_FUNCTION = HLL
            ELSE IF (FLUX_FUNCTION_STRING == "SLAU") THEN
               FLUX_FUNCTION = SLAU
            ELSE
               CALL ERROR_ABORT('Specified flux function in input file does not exist.')
            END IF
         END IF


         IF (line=='Read_restart:') THEN
            READ(in1,'(A)') RESTART_FILENAME
            BOOL_RESTART = .TRUE.
         END IF

      END DO ! Loop for reading input file

      CLOSE(in1) ! Close input file

      
   END SUBROUTINE READINPUT




   SUBROUTINE READ_WALL_REACTIONS(FILENAME)

      IMPLICIT NONE

      CHARACTER*64, INTENT(IN) :: FILENAME
      
      INTEGER, PARAMETER :: in4 = 5657
      INTEGER            :: ios
      CHARACTER*512      :: line
      INTEGER            :: ReasonEOF

      CHARACTER(LEN=80)   :: DEFINITION
      INTEGER :: N_STR, I
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)
      TYPE(WALL_REACTIONS_DATA_STRUCTURE) :: NEW_REACTION
      
      ! Open input file for reading
      OPEN(UNIT=in4,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios.NE.0) THEN
         CALL ERROR_ABORT('Attention, wall reactions definition file not found! ABORTING.')
      ENDIF

      line = '' ! Init empty

      ! +++++++ Read until the end of file ++++++++
      DO

         READ(in4,'(A)', IOSTAT=ReasonEOF) DEFINITION ! Read reaction components line         
         CALL STRIP_COMMENTS(DEFINITION, '!')         ! Remove comments from line

         IF (ReasonEOF < 0) EXIT ! End of file reached

         ! ~~~~~~~~~~~~~  Geometry and computational domain  ~~~~~~~~~~~~~~~~~
         
         CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)

         IF (N_STR == 0) CYCLE ! This is an empty line
         IF (STRARRAY(2) .NE. '-->' ) THEN
            CALL ERROR_ABORT('Attention, format is not respected in wall reactions file.')
         END IF

         NEW_REACTION%R_SP_ID = SPECIES_NAME_TO_ID(STRARRAY(1))
         IF (STRARRAY(3) .EQ. 'none' ) THEN
            NEW_REACTION%N_PROD = 0
            IF (N_STR .NE. 3) THEN
               CALL ERROR_ABORT('Attention, format is not respected in wall reactions file.')
            END IF
         ELSE
            NEW_REACTION%N_PROD = (N_STR-1)/2
            DO I = 1, NEW_REACTION%N_PROD
               IF (I == 1) THEN
                  NEW_REACTION%P1_SP_ID = SPECIES_NAME_TO_ID( STRARRAY(3) )
               ELSE IF (I == 2) THEN
                  NEW_REACTION%P2_SP_ID = SPECIES_NAME_TO_ID(STRARRAY(5))
                  IF (STRARRAY(4) .NE. '+' ) THEN
                     CALL ERROR_ABORT('Attention, format is not respected in wall reactions file.')
                  END IF
               ELSE
                  ! Only acceps up to two products.
                  CALL ERROR_ABORT('Attention, format is not respected in wall reactions file.')
               END IF
            END DO
         END IF
      
         READ(in4,'(A)', IOSTAT=ReasonEOF) DEFINITION ! Read reaction parameters line
         CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)
         READ(STRARRAY(1), '(ES14.0)') NEW_REACTION%PROB

         IF (ReasonEOF < 0) EXIT ! End of file reached
         
         IF (ALLOCATED(WALL_REACTIONS)) THEN
            ALLOCATE(TEMP_WALL_REACTIONS(N_WALL_REACTIONS+1)) ! Append the mixture to the list
            TEMP_WALL_REACTIONS(1:N_WALL_REACTIONS) = WALL_REACTIONS(1:N_WALL_REACTIONS)
            CALL MOVE_ALLOC(TEMP_WALL_REACTIONS, WALL_REACTIONS)
         ELSE
            ALLOCATE(WALL_REACTIONS(1))
         END IF
         N_WALL_REACTIONS = N_WALL_REACTIONS + 1
         WALL_REACTIONS(N_WALL_REACTIONS) = NEW_REACTION

         DEALLOCATE(STRARRAY)

      END DO
      
      CLOSE(in4) ! Close input file

      DO I = 1, N_WALL_REACTIONS
         WRITE(*,*) 'Wall reaction ', I, 'Has 1 reactant with id:', WALL_REACTIONS(I)%R_SP_ID
         WRITE(*,*) WALL_REACTIONS(I)%N_PROD, 'products with ids:',  WALL_REACTIONS(I)%P1_SP_ID, ' and ', WALL_REACTIONS(I)%P2_SP_ID
         WRITE(*,*) 'Reaction probability:', WALL_REACTIONS(I)%PROB
      END DO

   END SUBROUTINE READ_WALL_REACTIONS



   SUBROUTINE READ_REACTIONS(FILENAME)

      IMPLICIT NONE

      CHARACTER*256, INTENT(IN) :: FILENAME
      
      INTEGER, PARAMETER :: in3 = 457
      INTEGER, PARAMETER :: in4 = 754
      INTEGER            :: ios
      CHARACTER*512      :: line
      INTEGER            :: ReasonEOF, ReasonEOFCS

      CHARACTER(LEN=80)   :: DEFINITION
      INTEGER :: N_STR
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)
      TYPE(REACTIONS_DATA_STRUCTURE) :: NEW_REACTION
      INTEGER :: I, NROWS, INDEX
      CHARACTER*512      :: LINECS
      CHARACTER*256      :: REACTION_FILENAME
      
      !CHARACTER*64 :: SP_NAME
      !INTEGER      :: SP_ID
      !REAL(KIND=8) :: DIAM
      !REAL(KIND=8) :: OMEGA
      !REAL(KIND=8) :: TREF
      !REAL(KIND=8) :: ALPHA

      !INTEGER      :: IS, JS
      !REAL(KIND=8) :: M1, M2, MRED


      ! Open input file for reading
      OPEN(UNIT=in3,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios.NE.0) THEN
         CALL ERROR_ABORT('Attention, reactions definition file not found! ABORTING.')
      ENDIF

      line = '' ! Init empty

      ! +++++++ Read until the end of file ++++++++
      DO

         READ(in3,'(A)', IOSTAT=ReasonEOF) DEFINITION ! Read reaction components line         
         CALL STRIP_COMMENTS(DEFINITION, '!')         ! Remove comments from line

         IF (ReasonEOF < 0) EXIT

         CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)

         IF (N_STR == 0) CYCLE ! This is an empty line
         IF (STRARRAY(2) .NE. '+' .OR. (STRARRAY(4) .NE. '-->') .OR. STRARRAY(6) .NE. '+') THEN
            WRITE(*,*) 'Line in reactions file:', DEFINITION
            CALL ERROR_ABORT('Attention, format is not respected in reactions file.')
         END IF

         IF (N_STR .GE. 7) THEN
            NEW_REACTION%N_PROD = 2
            NEW_REACTION%R1_SP_ID = SPECIES_NAME_TO_ID(STRARRAY(1))
            NEW_REACTION%R2_SP_ID = SPECIES_NAME_TO_ID(STRARRAY(3))
            NEW_REACTION%P1_SP_ID = SPECIES_NAME_TO_ID(STRARRAY(5))
            NEW_REACTION%P2_SP_ID = SPECIES_NAME_TO_ID(STRARRAY(7))
         END IF

         IF (N_STR .GE. 9) THEN
            IF (STRARRAY(8) .NE. '+') THEN
               WRITE(*,*) 'Line in reactions file:', DEFINITION
               CALL ERROR_ABORT('Attention, format is not respected in reactions file.')
            END IF
            NEW_REACTION%N_PROD = 3
            NEW_REACTION%P3_SP_ID = SPECIES_NAME_TO_ID(STRARRAY(9))
         END IF

         IF (N_STR .GE. 11) THEN
            IF (STRARRAY(10) .NE. '+') THEN
               WRITE(*,*) 'Line in reactions file:', DEFINITION
               CALL ERROR_ABORT('Attention, format is not respected in reactions file.')
            END IF
            NEW_REACTION%N_PROD = 4
            NEW_REACTION%P4_SP_ID = SPECIES_NAME_TO_ID(STRARRAY(11))
         END IF

         READ(in3,'(A)', IOSTAT=ReasonEOF) DEFINITION ! Read reaction parameters line
         IF (ReasonEOF < 0) EXIT
         CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)
         IF (STRARRAY(1) == 'constant') THEN
            NEW_REACTION%TYPE = FIXED_RATE
            READ(STRARRAY(2), *) NEW_REACTION%CONSTANT_RATE
         ELSE IF (STRARRAY(1) == 'hardsphere') THEN
            NEW_REACTION%TYPE = HARD_SPHERE
            READ(STRARRAY(2), *) NEW_REACTION%DIAM
            READ(STRARRAY(3), *) NEW_REACTION%TA
         ELSE IF (STRARRAY(1) == 'arrhenius') THEN
            NEW_REACTION%TYPE = ARRHENIUS
            READ(STRARRAY(2), *) NEW_REACTION%A
            READ(STRARRAY(3), *) NEW_REACTION%N
            READ(STRARRAY(4), *) NEW_REACTION%TA
         ELSE IF (STRARRAY(1) == 'tabulated') THEN
            NEW_REACTION%TYPE = TABULATED
            READ(STRARRAY(2), *) REACTION_FILENAME


            OPEN(UNIT=in4,FILE=REACTION_FILENAME, STATUS='old',IOSTAT=ios)

            IF (ios .NE. 0) THEN
               CALL ERROR_ABORT('Attention, reactions rate file not found! ABORTING.')
            ENDIF
      
            LINECS = '' ! Init empty

            DO
               READ(in4,'(A)', IOSTAT=ReasonEOFCS) LINECS
               IF (ReasonEOFCS < 0) CALL ERROR_ABORT('Attention, reactions rate file format error! ABORTING.')
                  
               IF (LINECS(1:5) == '-----') EXIT
            END DO
            NROWS = 0
            DO
               READ(in4,'(A)', IOSTAT=ReasonEOFCS) LINECS  
               IF (ReasonEOFCS < 0) CALL ERROR_ABORT('Attention, reactions rate file format error! ABORTING.')
                  
               IF (LINECS(1:5) == '-----') EXIT
               NROWS = NROWS + 1
            END DO
            REWIND(in4)
            IF (ALLOCATED(NEW_REACTION%TABLE_TEMP)) DEALLOCATE(NEW_REACTION%TABLE_TEMP)
            IF (ALLOCATED(NEW_REACTION%TABLE_RATE)) DEALLOCATE(NEW_REACTION%TABLE_RATE)
            ALLOCATE(NEW_REACTION%TABLE_TEMP(NROWS))
            ALLOCATE(NEW_REACTION%TABLE_RATE(NROWS))
            DO
               READ(in4,'(A)', IOSTAT=ReasonEOFCS) LINECS
               IF (ReasonEOFCS < 0) CALL ERROR_ABORT('Attention, reactions rate file format error! ABORTING.')
                  
               IF (LINECS(1:5) == '-----') EXIT
            END DO
            DO I = 1, NROWS
               READ(in4,'(A)', IOSTAT=ReasonEOFCS) LINECS
               IF (ReasonEOFCS < 0) CALL ERROR_ABORT('Attention, reactions rate file format error! ABORTING.')
               READ(LINECS, *) NEW_REACTION%TABLE_TEMP(I), NEW_REACTION%TABLE_RATE(I)
            END DO

            CLOSE(in4)

         END IF


         READ(in3,'(A)', IOSTAT=ReasonEOF) DEFINITION ! Read reaction parameters line
         CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)
         DO I = 1, MIN(N_STR, 6)
            READ(STRARRAY(I), *) NEW_REACTION%DELTAE(I)
         END DO


         IF (ReasonEOF < 0) EXIT ! End of file reached
         
         IF (ALLOCATED(REACTIONS)) THEN
            ALLOCATE(TEMP_REACTIONS(N_REACTIONS+1)) ! Append the reaction to the list
            TEMP_REACTIONS(1:N_REACTIONS) = REACTIONS(1:N_REACTIONS)
            CALL MOVE_ALLOC(TEMP_REACTIONS, REACTIONS)
         ELSE
            ALLOCATE(REACTIONS(1))
         END IF
         N_REACTIONS = N_REACTIONS + 1
         REACTIONS(N_REACTIONS) = NEW_REACTION

         DEALLOCATE(STRARRAY)

      END DO
      
      CLOSE(in3) ! Close input file


      DO INDEX = 1, N_REACTIONS
         WRITE(*,*) 'Reaction ', index, 'Has 2 reactants with ids:', REACTIONS(index)%R1_SP_ID, ' and ', &
         REACTIONS(index)%R2_SP_ID
         IF (REACTIONS(index)%N_PROD == 2) THEN
            WRITE(*,*) REACTIONS(index)%N_PROD, 'products with ids:',  REACTIONS(index)%P1_SP_ID, ' and ', &
            REACTIONS(index)%P2_SP_ID
         ELSE IF (REACTIONS(index)%N_PROD == 3) THEN
            WRITE(*,*) REACTIONS(index)%N_PROD, 'products with ids:',  REACTIONS(index)%P1_SP_ID, ', ', &
            REACTIONS(index)%P2_SP_ID, ' and ', REACTIONS(index)%P3_SP_ID
         END IF
         IF (REACTIONS(index)%TYPE == FIXED_RATE) THEN
            WRITE(*,*) 'Fixed rate reaction with rate:', REACTIONS(index)%CONSTANT_RATE
         ELSE IF (REACTIONS(index)%TYPE == HARD_SPHERE) THEN
            WRITE(*,*) 'Hard-sphere type with diameter ', REACTIONS(index)%DIAM, &
            ' and activation temperature ', REACTIONS(index)%TA
         ELSE IF (REACTIONS(index)%TYPE == ARRHENIUS) THEN
            WRITE(*,*) 'Arrhenius type with factor ', REACTIONS(index)%A, &
            ' exponent ', REACTIONS(index)%N, ' and activation  temperature ', REACTIONS(index)%TA
         ELSE IF (REACTIONS(index)%TYPE == TABULATED) THEN
            WRITE(*,*) 'Reaction ', index, ' is from tabulated data.'
            WRITE(*,*) 'The temperature goes from ', MINVAL(REACTIONS(index)%TABLE_TEMP), &
                       ' to ', MAXVAL(REACTIONS(index)%TABLE_TEMP)
            WRITE(*,*) 'The rate goes from ', MINVAL(REACTIONS(index)%TABLE_RATE), &
                       ' to ', MAXVAL(REACTIONS(index)%TABLE_RATE)
         END IF
         WRITE(*,*) 'The energy balance for each component are ', REACTIONS(INDEX)%DELTAE(1:2+REACTIONS(INDEX)%N_PROD)
      END DO


   END SUBROUTINE READ_REACTIONS



   SUBROUTINE DEF_BOUNDARY_CONDITION(DEFINITION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: DEFINITION

      INTEGER :: N_STR, I, IPG, SP_ID, FIRST, LAST
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)
      CHARACTER(LEN=64) :: MIX_NAME
      REAL(KIND=8), DIMENSION(Neq) :: U_BOUND_PRIM
      REAL(KIND=8) :: FRAC


      CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)

      ! phys_group type parameters
      IPG = -1
      DO I = 1, N_GRID_BC
         IF (GRID_BC(I)%PHYSICAL_GROUP_NAME == STRARRAY(1)) IPG = I
      END DO
      IF (IPG == -1) THEN
         WRITE(*,*) 'Group ', STRARRAY(1), ' not found.'
         CALL ERROR_ABORT('Error in boundary condition definition. Group name not found.')
      END IF

      IF (STRARRAY(2) == 'state') THEN
         GRID_BC(IPG)%PARTICLE_BC = STATE
         READ(STRARRAY(3),'(A10)') MIX_NAME
         GRID_BC(IPG)%MIX_ID = MIXTURE_NAME_TO_ID(MIX_NAME)
         READ(STRARRAY(4), '(ES14.0)') GRID_BC(IPG)%NRHO
         READ(STRARRAY(5), '(ES14.0)') GRID_BC(IPG)%UX
         READ(STRARRAY(6), '(ES14.0)') GRID_BC(IPG)%UY
         READ(STRARRAY(7), '(ES14.0)') GRID_BC(IPG)%UZ
         READ(STRARRAY(8), '(ES14.0)') GRID_BC(IPG)%TEMP

         ALLOCATE(GRID_BC(IPG)%U_BOUND(N_SPECIES_FLUID*Neq))
         GRID_BC(IPG)%U_BOUND = 0.d0
         DO I = 1, MIXTURES(GRID_BC(IPG)%MIX_ID)%N_COMPONENTS
            SP_ID = MIXTURES(GRID_BC(IPG)%MIX_ID)%COMPONENTS(I)%ID
            FIRST = (SP_ID-1)*Neq+1
            LAST = SP_ID*Neq+1

            FRAC = MIXTURES(GRID_BC(IPG)%MIX_ID)%COMPONENTS(I)%MOLFRAC
            U_BOUND_PRIM(1) = FRAC*GRID_BC(IPG)%NRHO*SPECIES(SP_ID)%MOLECULAR_MASS
            U_BOUND_PRIM(2) = U_BOUND_PRIM(1)*GRID_BC(IPG)%UX
            U_BOUND_PRIM(3) = U_BOUND_PRIM(1)*GRID_BC(IPG)%UY
            U_BOUND_PRIM(4) = GRID_BC(IPG)%TEMP
            CALL compute_primitive_from_conserved(U_BOUND_PRIM, GRID_BC(IPG)%U_BOUND(FIRST:LAST), SP_ID)
         END DO
      ELSE IF (STRARRAY(2) == 'wall') THEN
         GRID_BC(IPG)%PARTICLE_BC = WALL
      ELSE IF (STRARRAY(2) == 'symmetry') THEN
         GRID_BC(IPG)%PARTICLE_BC = SYMMETRY
      ELSE IF (STRARRAY(2) == 'react') THEN
         GRID_BC(IPG)%REACT = .TRUE.
      ELSE
         CALL ERROR_ABORT('Error in boundary condition definition.')
      END IF

   END SUBROUTINE DEF_BOUNDARY_CONDITION



   SUBROUTINE DEF_DOMAIN_TYPE(DEFINITION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: DEFINITION

      INTEGER :: N_STR, I, IPG
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)


      CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)

      ! phys_group type parameters
      IPG = -1
      DO I = 1, N_GRID_BC
         IF (GRID_BC(I)%PHYSICAL_GROUP_NAME == STRARRAY(1)) IPG = I
      END DO
      IF (IPG == -1) THEN
         WRITE(*,*) 'Group ', STRARRAY(1), ' not found.'
         CALL ERROR_ABORT('Error in domain type definition. Group name not found.')
      END IF
      
      IF (STRARRAY(2) == 'fluid') THEN
         GRID_BC(IPG)%VOLUME_BC = FLUID
      ELSE IF (STRARRAY(2) == 'solid') THEN
         GRID_BC(IPG)%VOLUME_BC = SOLID
      ELSE
         CALL ERROR_ABORT('Error in boundary condition definition.')
      END IF

   END SUBROUTINE DEF_DOMAIN_TYPE



   SUBROUTINE READ_SPECIES(FILENAME)

      IMPLICIT NONE

      CHARACTER*256, INTENT(IN) :: FILENAME

      INTEGER, PARAMETER :: in2 = 445
      INTEGER            :: ios
      CHARACTER*512      :: line
      INTEGER            :: ReasonEOF

      CHARACTER*10 :: NAME
      REAL(KIND=8) :: MOLWT
      REAL(KIND=8) :: MOLECULAR_MASS
      REAL(KIND=8) :: DIAM
      REAL(KIND=8) :: GAMMA
      REAL(KIND=8) :: KAPPA
      REAL(KIND=8) :: MU

   
      ! Open input file for reading
      OPEN(UNIT=in2,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios.NE.0) THEN
         CALL ERROR_ABORT('Attention, species definition file not found! ABORTING.')
      ENDIF

      line = '' ! Init empty

      ! +++++++ Read until the end of file ++++++++
      N_SPECIES = 0
      DO

         READ(in2,'(A)', IOSTAT=ReasonEOF) line ! Read line         
         CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

         IF (ReasonEOF < 0) EXIT ! End of file reached

         ! ~~~~~~~~~~~~~  Geometry and computational domain  ~~~~~~~~~~~~~~~~~
         
         !READ(line,'(A2, ES14.3, ES14.3, I1, ES14.3, I1, ES14.3, ES14.3, ES14.3, ES14.3)') &
         READ(line,*) &
         NAME, MOLWT, MOLECULAR_MASS, DIAM, GAMMA, KAPPA, MU
      
         IF (ALLOCATED(SPECIES)) THEN
            ALLOCATE(TEMP_SPECIES(N_SPECIES+1)) ! Append the species to the list
            TEMP_SPECIES(1:N_SPECIES) = SPECIES(1:N_SPECIES)
            CALL MOVE_ALLOC(TEMP_SPECIES, SPECIES)
         ELSE
            ALLOCATE(SPECIES(1))
         END IF
         N_SPECIES = N_SPECIES + 1

         SPECIES(N_SPECIES)%NAME = NAME
         SPECIES(N_SPECIES)%MOLECULAR_MASS = MOLECULAR_MASS
         SPECIES(N_SPECIES)%DIAM = DIAM
         SPECIES(N_SPECIES)%GAMMA = GAMMA
         SPECIES(N_SPECIES)%KAPPA = KAPPA
         SPECIES(N_SPECIES)%MU = MU
         SPECIES(N_SPECIES)%CP = GAMMA/(GAMMA-1)*KB/MOLECULAR_MASS

      END DO

      CLOSE(in2) ! Close input file

   END SUBROUTINE READ_SPECIES



   SUBROUTINE DEF_MIXTURE(DEFINITION)
      
      IMPLICIT NONE

      CHARACTER*64 :: MIX_NAME
      CHARACTER*64 :: COMP_NAME
      REAL(KIND=8) :: MOLFRAC

      CHARACTER(LEN=*), INTENT(IN) :: DEFINITION
      INTEGER :: N_STR
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)
      INTEGER :: i, N_COMP
      !TYPE(MIXTURE_COMPONENT) NEW_MIXTURE
      TYPE(MIXTURE), DIMENSION(:), ALLOCATABLE :: TEMP_MIXTURES
      TYPE(MIXTURE_COMPONENT), DIMENSION(:), ALLOCATABLE :: TEMP_COMPONENTS


      CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)
      
      READ(STRARRAY(1),'(A10)') MIX_NAME
      N_COMP = (N_STR-1)/2
      ALLOCATE(TEMP_COMPONENTS(N_COMP))
      DO i = 1, N_COMP
         READ(STRARRAY(2*i),'(A10)') COMP_NAME
         READ(STRARRAY(2*i+1),'(ES14.3)') MOLFRAC

         TEMP_COMPONENTS(i)%NAME = COMP_NAME
         TEMP_COMPONENTS(i)%MOLFRAC = MOLFRAC
         TEMP_COMPONENTS(i)%ID = SPECIES_NAME_TO_ID(COMP_NAME)
      END DO


      IF (ALLOCATED(MIXTURES)) THEN
         ALLOCATE(TEMP_MIXTURES(N_MIXTURES+1)) ! Append the mixture to the list
         TEMP_MIXTURES(1:N_MIXTURES) = MIXTURES(1:N_MIXTURES)
         CALL MOVE_ALLOC(TEMP_MIXTURES, MIXTURES)
      ELSE
         ALLOCATE(MIXTURES(1))
      END IF
      N_MIXTURES = N_MIXTURES + 1

      MIXTURES(N_MIXTURES)%NAME = MIX_NAME
      MIXTURES(N_MIXTURES)%N_COMPONENTS = N_COMP

      CALL MOVE_ALLOC(TEMP_COMPONENTS, MIXTURES(N_MIXTURES)%COMPONENTS)


   END SUBROUTINE DEF_MIXTURE


   !=============================================================
   ! Helper routine to read ASCII lines from stream file
   !=============================================================
   logical function read_line(unit, line)
      integer, intent(in) :: unit
      character(len=*), intent(out) :: line
      integer :: ios
      character(len=1) :: ch
      integer :: i

      line = ''
      i = 1

      do
         read(unit, iostat=ios) ch
         if (ios /= 0) then
            read_line = .false.
            return
         end if

         if (ch == new_line('a')) exit
         if (i <= len(line)) then
            line(i:i) = ch
            i = i + 1
         end if
      end do

      read_line = .true.
   end function read_line


   SUBROUTINE READ_BACKGROUND_FILE(FILENAME)

      IMPLICIT NONE

      CHARACTER*512, INTENT(IN) :: FILENAME
      CHARACTER*30       :: line
      INTEGER, PARAMETER :: in5 = 5557
      INTEGER            :: ios
      INTEGER            :: ReasonEOF
      INTEGER            :: N_STR, SP_ID
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)

      ! Open input file for reading
      !OPEN(UNIT=in5,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (BOOL_BINARY_OUTPUT) THEN
         OPEN(UNIT=in5, FILE=FILENAME, STATUS='OLD', IOSTAT=ios, ACCESS='STREAM', CONVERT='BIG_ENDIAN')
      ELSE
         OPEN(UNIT=in5, FILE=FILENAME, STATUS='OLD', IOSTAT=ios, FORM='FORMATTED')
      END IF

      IF (ios.NE.0) THEN
         CALL ERROR_ABORT('Attention, vtk file not found! ABORTING.')
      END IF

      ALLOCATE(BG_CELL_NRHO(N_SPECIES, NCELLS))
      ALLOCATE(BG_CELL_TEMP(N_SPECIES, NCELLS))
      ALLOCATE(BG_CELL_VX(N_SPECIES, NCELLS))
      ALLOCATE(BG_CELL_VY(N_SPECIES, NCELLS))
      ALLOCATE(BG_CELL_VZ(N_SPECIES, NCELLS))

      ! ++++++++ Read until the end of file ++++++++
 
      DO
         IF (BOOL_BINARY_OUTPUT) THEN
            IF (.NOT. READ_LINE(in5, line)) EXIT
         ELSE
            READ(in5, IOSTAT=ReasonEOF) line ! Read line
            IF (ReasonEOF < 0) EXIT ! End of file reached
         END IF
      
         CALL SPLIT_STR(line, ' ', STRARRAY, N_STR)

         IF (STRARRAY(1)(1:10) == 'nrho_mean_') THEN
            SP_ID = SPECIES_NAME_TO_ID(STRARRAY(1)(11:))
            IF (SP_ID == -1) CALL ERROR_ABORT('Error! Species in MCC background vtk file not found.')
            WRITE(*,*) 'Found n data for species ', SP_ID
            READ(in5, IOSTAT=ReasonEOF) BG_CELL_NRHO(SP_ID, :)
            IF (ReasonEOF < 0) EXIT ! End of file reached
         END IF

         IF (STRARRAY(1)(1:9) == 'Ttr_mean_') THEN
            SP_ID = SPECIES_NAME_TO_ID(STRARRAY(1)(10:))
            IF (SP_ID == -1) CALL ERROR_ABORT('Error! Species in MCC background vtk file not found.')
            WRITE(*,*) 'Found T data for species ', SP_ID
            READ(in5, IOSTAT=ReasonEOF) BG_CELL_TEMP(SP_ID, :)
            IF (ReasonEOF < 0) EXIT ! End of file reached
         END IF

         IF (STRARRAY(1)(1:8) == 'vx_mean_') THEN
            SP_ID = SPECIES_NAME_TO_ID(STRARRAY(1)(9:))
            IF (SP_ID == -1) CALL ERROR_ABORT('Error! Species in MCC background vtk file not found.')
            WRITE(*,*) 'Found vx data for species ', SP_ID
            READ(in5, IOSTAT=ReasonEOF) BG_CELL_VX(SP_ID, :)
            IF (ReasonEOF < 0) EXIT ! End of file reached
         END IF

         IF (STRARRAY(1)(1:8) == 'vy_mean_') THEN
            SP_ID = SPECIES_NAME_TO_ID(STRARRAY(1)(9:))
            IF (SP_ID == -1) CALL ERROR_ABORT('Error! Species in MCC background vtk file not found.')
            WRITE(*,*) 'Found vy data for species ', SP_ID
            READ(in5, IOSTAT=ReasonEOF) BG_CELL_VY(SP_ID, :)
            IF (ReasonEOF < 0) EXIT ! End of file reached
         END IF

         IF (STRARRAY(1)(1:8) == 'vz_mean_') THEN
            SP_ID = SPECIES_NAME_TO_ID(STRARRAY(1)(9:))
            IF (SP_ID == -1) CALL ERROR_ABORT('Error! Species in MCC background vtk file not found.')
            WRITE(*,*) 'Found vz data for species ', SP_ID
            READ(in5, IOSTAT=ReasonEOF) BG_CELL_VZ(SP_ID, :)
            IF (ReasonEOF < 0) EXIT ! End of file reached
         END IF
      END DO

      CLOSE(in5) ! Close input file

      BOOL_BG_FILE = .TRUE.

      WRITE(*,*) 'Done reading.'

   END SUBROUTINE READ_BACKGROUND_FILE



   SUBROUTINE READ_VSS(FILENAME)

      IMPLICIT NONE

      CHARACTER*64, INTENT(IN) :: FILENAME
      
      INTEGER, PARAMETER :: in2 = 456
      INTEGER            :: ios
      CHARACTER*512      :: line
      INTEGER            :: ReasonEOF

      CHARACTER*64 :: SP_NAME
      INTEGER      :: SP_ID
      REAL(KIND=8) :: DIAM
      REAL(KIND=8) :: OMEGA
      REAL(KIND=8) :: TREF
      REAL(KIND=8) :: ALPHA

      INTEGER      :: IS, JS
      REAL(KIND=8) :: M1, M2, MRED

      ! Open input file for reading
      OPEN(UNIT=in2,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios.NE.0) THEN
         CALL ERROR_ABORT('Attention, VSS parameters definition file not found! ABORTING.')
      ENDIF

      line = '' ! Init empty

      ! +++++++ Read until the end of file ++++++++
      DO

         READ(in2,'(A)', IOSTAT=ReasonEOF) line ! Read line         
         CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

         IF (ReasonEOF < 0) EXIT ! End of file reached

         ! ~~~~~~~~~~~~~  Geometry and computational domain  ~~~~~~~~~~~~~~~~~
         
         !READ(line,'(A2, ES14.3, ES14.3, I1, ES14.3, I1, ES14.3, ES14.3, ES14.3, ES14.3)') &
         READ(line,*) SP_NAME, DIAM, OMEGA, TREF, ALPHA
      
         SP_ID = SPECIES_NAME_TO_ID(SP_NAME)

         SPECIES(SP_ID)%DIAM  = DIAM
         SPECIES(SP_ID)%OMEGA = OMEGA
         SPECIES(SP_ID)%TREF  = TREF
         SPECIES(SP_ID)%ALPHA = ALPHA

         SPECIES(SP_ID)%SIGMA = PI*DIAM**2
         
      END DO

      ALLOCATE(VSS_GREFS(N_SPECIES, N_SPECIES))
      ALLOCATE(VSS_SIGMAS(N_SPECIES, N_SPECIES))
      ALLOCATE(VSS_ALPHAS(N_SPECIES, N_SPECIES))
      ALLOCATE(VSS_OMEGAS(N_SPECIES, N_SPECIES))


      DO IS = 1, N_SPECIES
         DO JS = 1, N_SPECIES
            M1    = SPECIES(IS)%MOLECULAR_MASS
            M2    = SPECIES(JS)%MOLECULAR_MASS
            MRED  = M1*M2/(M1+M2)
            OMEGA    = 0.5 * (SPECIES(IS)%OMEGA + SPECIES(JS)%OMEGA)
            TREF = 0.5 * (SPECIES(IS)%TREF + SPECIES(JS)%TREF)
            
            IF (ABS(OMEGA-0.5) .LT. 1.d-6) THEN
               VSS_GREFS(IS, JS) = (2.*KB*TREF/MRED)**0.5 * 1.2354
            ELSE
               VSS_GREFS(IS, JS) = (2.*KB*TREF/MRED)**0.5 * (GAMMA(2.5-OMEGA))**(-0.5/(OMEGA-0.5))
            END IF

            VSS_SIGMAS(IS, JS) = PI * (0.5 * (SPECIES(IS)%DIAM + SPECIES(JS)%DIAM))**2
            VSS_ALPHAS(IS, JS) = 0.5 * (SPECIES(IS)%ALPHA + SPECIES(JS)%ALPHA)
            VSS_OMEGAS(IS, JS) = 0.5 * (SPECIES(IS)%OMEGA + SPECIES(JS)%OMEGA)

         END DO
      END DO

      CLOSE(in2) ! Close input file


   END SUBROUTINE READ_VSS




   SUBROUTINE READ_VSS_BINARY(FILENAME)

      IMPLICIT NONE

      CHARACTER*64, INTENT(IN) :: FILENAME
      
      INTEGER, PARAMETER :: in2 = 457
      INTEGER            :: ios
      CHARACTER*512      :: line
      INTEGER            :: ReasonEOF

      INTEGER :: N_STR
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)

      INTEGER, ALLOCATABLE :: SP_IDS(:)
      REAL(KIND=8) :: OMEGA
      REAL(KIND=8) :: TREF
      REAL(KIND=8) :: READ_VALUE

      INTEGER      :: IS, JS
      REAL(KIND=8) :: M1, M2, MRED

      ALLOCATE(VSS_GREFS(N_SPECIES, N_SPECIES))
      ALLOCATE(VSS_SIGMAS(N_SPECIES, N_SPECIES))
      ALLOCATE(VSS_ALPHAS(N_SPECIES, N_SPECIES))
      ALLOCATE(VSS_OMEGAS(N_SPECIES, N_SPECIES))

      ! Open input file for reading
      OPEN(UNIT=in2,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios.NE.0) THEN
         CALL ERROR_ABORT('Attention, VSS parameters definition file not found! ABORTING.')
      ENDIF

      ! Read Tref
      line = '' ! Init empty
      DO WHILE (line == '')
         
         READ(in2,'(A)', IOSTAT=ReasonEOF) line ! Read line         
         CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

         IF (ReasonEOF < 0) THEN
            CALL ERROR_ABORT('Attention, VSS definition reading failed! ABORTING.')
         END IF
         
      END DO

      READ(line,*) TREF
      
      ! Read species ordering
      N_STR = 0
      DO WHILE (N_STR == 0)
         line = '' ! Init empty
         READ(in2,'(A)', IOSTAT=ReasonEOF) line ! Read line         
         CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

         IF (ReasonEOF < 0) THEN
            CALL ERROR_ABORT('Attention, VSS definition reading failed! ABORTING.')
         END IF

         CALL SPLIT_STR(line, ' ', STRARRAY, N_STR)
      END DO

      IF (N_STR .NE. N_SPECIES) CALL ERROR_ABORT('Attention, incorrect number of species in VSS definition! ABORTING.')

      ALLOCATE(SP_IDS(N_SPECIES))
      DO IS = 1, N_SPECIES
         SP_IDS(IS) = SPECIES_NAME_TO_ID(STRARRAY(IS))
      END DO


      ! Read molecular diameters
      DO JS = 1, N_SPECIES
         N_STR = 0
         DO WHILE (N_STR == 0)
            line = '' ! Init empty
            READ(in2,'(A)', IOSTAT=ReasonEOF) line ! Read line         
            CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

            IF (ReasonEOF < 0) THEN
               CALL ERROR_ABORT('Attention, VSS definition reading failed! ABORTING.')
            END IF

            CALL SPLIT_STR(line, ' ', STRARRAY, N_STR)
         END DO
         
         IF (N_STR .NE. N_SPECIES+1-JS) CALL ERROR_ABORT('Attention, VSS params matrix format wrong! ABORTING.')
         

         DO IS = JS, N_SPECIES
            READ(STRARRAY(1+IS-JS),*) READ_VALUE
            VSS_SIGMAS(SP_IDS(IS), SP_IDS(JS)) = PI*READ_VALUE**2
            VSS_SIGMAS(SP_IDS(JS), SP_IDS(IS)) = VSS_SIGMAS(SP_IDS(IS), SP_IDS(JS))
         END DO
      END DO
         

      ! Read omegas
      DO JS = 1, N_SPECIES
         N_STR = 0
         DO WHILE (N_STR == 0)
            line = '' ! Init empty
            READ(in2,'(A)', IOSTAT=ReasonEOF) line ! Read line         
            CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

            IF (ReasonEOF < 0) THEN
               CALL ERROR_ABORT('Attention, VSS definition reading failed! ABORTING.')
            END IF

            CALL SPLIT_STR(line, ' ', STRARRAY, N_STR)
         END DO

         DO IS = JS, N_SPECIES
            READ(STRARRAY(1+IS-JS),*) READ_VALUE
            VSS_OMEGAS(SP_IDS(IS), SP_IDS(JS)) = READ_VALUE
            VSS_OMEGAS(SP_IDS(JS), SP_IDS(IS)) = VSS_OMEGAS(SP_IDS(IS), SP_IDS(JS))
         END DO
      END DO


      ! Read alphas
      DO JS = 1, N_SPECIES
         N_STR = 0
         DO WHILE (N_STR == 0)
            line = '' ! Init empty
            READ(in2,'(A)', IOSTAT=ReasonEOF) line ! Read line         
            CALL STRIP_COMMENTS(line, '!')         ! Remove comments from line

            IF (ReasonEOF < 0) THEN
               CALL ERROR_ABORT('Attention, VSS definition reading failed! ABORTING.')
            END IF

            CALL SPLIT_STR(line, ' ', STRARRAY, N_STR)
         END DO

         DO IS = JS, N_SPECIES
            READ(STRARRAY(1+IS-JS),*) READ_VALUE
            VSS_ALPHAS(SP_IDS(IS), SP_IDS(JS)) = READ_VALUE
            VSS_ALPHAS(SP_IDS(JS), SP_IDS(IS)) = VSS_ALPHAS(SP_IDS(IS), SP_IDS(JS))
         END DO
      END DO
 

      DO IS = 1, N_SPECIES
         DO JS = 1, N_SPECIES
            M1    = SPECIES(IS)%MOLECULAR_MASS
            M2    = SPECIES(JS)%MOLECULAR_MASS
            MRED  = M1*M2/(M1+M2)
            OMEGA = VSS_OMEGAS(IS, JS)
                        
            IF (ABS(OMEGA-0.5) .LT. 1.d-6) THEN
               VSS_GREFS(IS, JS) = (2.*KB*TREF/MRED)**0.5 * 1.2354
            ELSE
               VSS_GREFS(IS, JS) = (2.*KB*TREF/MRED)**0.5 * (GAMMA(2.5-OMEGA))**(-0.5/(OMEGA-0.5))
            END IF
         END DO
      END DO

      CLOSE(in2) ! Close input file

   

      ! WRITE(*,*) 'Tref = ', TREF
      ! WRITE(*,*) 'sigma = '
      ! DO IS = 1, N_SPECIES
      !    WRITE(*,*) VSS_SIGMAS(IS,:)
      ! END DO
      ! WRITE(*,*) 'omega = '
      ! DO IS = 1, N_SPECIES
      !    WRITE(*,*) VSS_OMEGAS(IS,:)
      ! END DO
      ! WRITE(*,*) 'alpha = '
      ! DO IS = 1, N_SPECIES
      !    WRITE(*,*) VSS_ALPHAS(IS,:)
      ! END DO
      ! WRITE(*,*) 'gref = '
      ! DO IS = 1, N_SPECIES
      !    WRITE(*,*) VSS_GREFS(IS,:)
      ! END DO

   END SUBROUTINE READ_VSS_BINARY



   SUBROUTINE DEF_INITIAL_STATE(DEFINITION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: DEFINITION

      INTEGER :: N_STR
      CHARACTER(LEN=80), ALLOCATABLE :: STRARRAY(:)

      CHARACTER*64 :: MIX_NAME
   
      CALL SPLIT_STR(DEFINITION, ' ', STRARRAY, N_STR)


      READ(STRARRAY(1),'(A10)') MIX_NAME
      READ(STRARRAY(2), '(ES14.0)') INITIAL_NRHO
      READ(STRARRAY(3), '(ES14.0)') INITIAL_UX
      READ(STRARRAY(4), '(ES14.0)') INITIAL_UY
      READ(STRARRAY(5), '(ES14.0)') INITIAL_UZ
      READ(STRARRAY(6), '(ES14.0)') INITIAL_TEMP

      INITIAL_MIX_ID = MIXTURE_NAME_TO_ID(MIX_NAME)

      ! WRITE(*,*) 'Read initial particle seed definition. Parameters: ', MIX_NAME, ', ', MIX_ID, ', ', NRHO, ', ',&
      !  UX, ', ', UY, ', ', UZ, ', ', TTRAX, ', ', TTRAY, ', ', TTRAZ, ', ', TROT, ', ', TVIB



   END SUBROUTINE DEF_INITIAL_STATE


   SUBROUTINE READ_RESTART(U)

      IMPLICIT NONE

      real(kind=8), dimension(:,:), intent(inout) :: U
      INTEGER :: ios

      character(len=100) :: iomsg

      ! Open file for reading
      WRITE(*,*) 'Opening file ', RESTART_FILENAME, ' for restart.'

      IF (BOOL_BINARY_OUTPUT) THEN
         OPEN(1025, FILE=RESTART_FILENAME, ACCESS='STREAM', FORM='UNFORMATTED', STATUS='OLD', &
         CONVERT='BIG_ENDIAN', IOSTAT=ios, IOMSG=iomsg)

         IF (ios .NE. 0) THEN
            WRITE(*,*) 'iomsg was: ', iomsg
            CALL ERROR_ABORT('Attention, restart file not found! ABORTING.')
         ENDIF

         READ(1025, IOSTAT=ios) U

         CLOSE(1025)
      ELSE
         OPEN(1025, FILE=RESTART_FILENAME, STATUS='OLD', IOSTAT=ios)
         
         IF (ios .NE. 0) THEN
            CALL ERROR_ABORT('Attention, restart file not found! ABORTING.')
         ENDIF

         READ(1025,*,IOSTAT=ios) U

         CLOSE(1025)
      END IF

   END SUBROUTINE

END MODULE initialization