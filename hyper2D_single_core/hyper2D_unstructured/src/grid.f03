module grid

   ! Most of the code in this module comes from the Pantera PIC-DSMC code,
   ! (https://github.com/vonkarmaninstitute/pantera-pic-dsmc) for reasons of compatibility.
   ! Pantera is a free software by the von Karman Institute for Fluid Dynamics (VKI), 
   ! distributed under a GNU GPLv3 license.

   use global_module

   implicit none

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Geometry, domain and grid !!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! INTEGER      :: NX, NY, NZ
   INTEGER      :: DIMS = 2
   ! REAL(KIND=8) :: CELL_VOL
   LOGICAL      :: AXI = .FALSE.
   ! LOGICAL      :: BOOL_X_PERIODIC = .FALSE. ! Assign default value!
   ! LOGICAL      :: BOOL_Y_PERIODIC = .FALSE.
   ! LOGICAL      :: BOOL_Z_PERIODIC = .FALSE.
   ! LOGICAL, DIMENSION(4) :: BOOL_PERIODIC = .FALSE.
   ! LOGICAL, DIMENSION(4) :: BOOL_SPECULAR = .FALSE.
   ! LOGICAL, DIMENSION(4) :: BOOL_REACT    = .FALSE.
   ! LOGICAL, DIMENSION(4) :: BOOL_DIFFUSE  = .FALSE.
   ! INTEGER, DIMENSION(:), ALLOCATABLE :: BOUNDARY_COLL_COUNT, WALL_COLL_COUNT, LINE_EMIT_COUNT
   ! REAL(KIND=8) :: BOUNDTEMP
   ! CHARACTER(LEN=256) :: GRID_FILENAME
   ! LOGICAL :: BOOL_BINARY_OUTPUT = .TRUE.
   ! REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: XCOORD, YCOORD
   ! REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: XSIZE, YSIZE
   ! INTEGER, DIMENSION(:), ALLOCATABLE :: CELL_PROCS
   INTEGER :: NCELLS, NNODES, NBOUNDCELLS, NBOUNDNODES
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: CELL_VOLUMES

   TYPE UNSTRUCTURED_1D_GRID_DATA_STRUCTURE
      INTEGER :: NUM_NODES, NUM_CELLS
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: NODE_COORDS
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_NODES
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_NEIGHBORS
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: EDGE_NORMAL
      INTEGER, DIMENSION(:), ALLOCATABLE        :: CELL_PG
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_EDGES_PG
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: PG_NODES
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: BASIS_COEFFS
      INTEGER, DIMENSION(:), ALLOCATABLE        :: PERIODIC_RELATED_NODE
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: SEGMENT_NODES_BOUNDARY_INDEX
      INTEGER, DIMENSION(:), ALLOCATABLE        :: NODES_BOUNDARY_INDEX
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: CELL_VOLUMES
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: SEGMENT_AREAS
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: SEGMENT_LENGTHS
   END TYPE UNSTRUCTURED_1D_GRID_DATA_STRUCTURE

   TYPE(UNSTRUCTURED_1D_GRID_DATA_STRUCTURE) :: U1D_GRID


   TYPE UNSTRUCTURED_2D_GRID_DATA_STRUCTURE
      INTEGER :: NUM_NODES, NUM_CELLS, NUM_LINES
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: NODE_COORDS
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_NODES
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: LINE_NODES
      INTEGER, DIMENSION(:), ALLOCATABLE        :: POINT_NODES
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_NEIGHBORS
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: EDGE_NORMAL
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: CELL_EDGES_LEN
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: CELL_FACES_AREA
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_EDGES_PG
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_EDGES_BOUNDARY_INDEX
      INTEGER, DIMENSION(:), ALLOCATABLE        :: NODES_BOUNDARY_INDEX
      INTEGER, DIMENSION(:), ALLOCATABLE        :: CELL_PG
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: PG_NODES
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: BASIS_COEFFS
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: LSTSQ_COEFFS
      INTEGER, DIMENSION(:), ALLOCATABLE        :: PERIODIC_RELATED_NODE
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: CELL_VOLUMES
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: CELL_AREAS
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: CELL_CENTROIDS
   END TYPE UNSTRUCTURED_2D_GRID_DATA_STRUCTURE

   TYPE(UNSTRUCTURED_2D_GRID_DATA_STRUCTURE) :: U2D_GRID



   CHARACTER*256 :: MESH_FILENAME

   INTEGER         :: N_GRID_BC = 0

   ENUM, BIND(C)  
      ENUMERATOR VACUUM, SPECULAR, DIFFUSE, CLL, REACT, AXIS, PERIODIC_MASTER, PERIODIC_SLAVE, EMIT, WB_BC, &
      PISTON, INLET, OUTLET, WALL, SYMMETRY
   END ENUM

   ENUM, BIND(C)
      ENUMERATOR DIRICHLET_BC, NEUMANN_BC, DIELECTRIC_BC, CONDUCTIVE_BC, THIN_DIELECTRIC_LAYER_BC, ROBIN_BC, & 
                 PERIODIC_MASTER_BC, PERIODIC_SLAVE_BC, RF_VOLTAGE_BC, DECOUPLED_RF_VOLTAGE_BC, SPICE_NODE_BC, NO_BC
   END ENUM

   ENUM, BIND(C)
      ENUMERATOR FLUID, SOLID
   END ENUM

   TYPE BOUNDARY_CONDITION_DATA_STRUCTURE
      CHARACTER(LEN=256)       :: PHYSICAL_GROUP_NAME
      INTEGER(KIND(VACUUM)) :: PARTICLE_BC = VACUUM !! DBDBDBD Maximum 100 species
      INTEGER(KIND(DIRICHLET_BC)) :: FIELD_BC = NO_BC
      INTEGER(KIND(FLUID)) :: VOLUME_BC = FLUID

      REAL(KIND=8) :: WALL_TEMP
      REAL(KIND=8) :: WALL_POTENTIAL
      REAL(KIND=8) :: WALL_EFIELD
      REAL(KIND=8) :: ACC_N
      REAL(KIND=8) :: ACC_T

      REAL(KIND=8) :: EPS_REL

      REAL(KIND=8) :: WALL_RF_POTENTIAL
      REAL(KIND=8) :: RF_FREQUENCY
      REAL(KIND=8) :: CAPACITANCE

      REAL(KIND=8) :: LAYER_THICKNESS

      REAL(KIND=8) :: SPICE_NODE_POTENTIAL = 0.d0
      REAL(KIND=8) :: SPICE_NODE_CURRENT

      REAL(KIND=8), DIMENSION(2) :: TRANSLATEVEC

      LOGICAL :: REACT = .FALSE.
      LOGICAL :: DUMP_FLUXES = .FALSE.

      ! Piston boundary condition (moving wall)
      REAL(KIND=8), DIMENSION(3) :: U_PISTON = 0.d0

   END TYPE BOUNDARY_CONDITION_DATA_STRUCTURE

   TYPE(BOUNDARY_CONDITION_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: GRID_BC






   contains 


   SUBROUTINE READ_2D_UNSTRUCTURED_GRID_SU2(FILENAME)

      IMPLICIT NONE

      CHARACTER*256, INTENT(IN) :: FILENAME

      CHARACTER*256 :: LINE, GROUPNAME, DUMMYLINE

      INTEGER, PARAMETER :: in5 = 2385
      INTEGER            :: ios
      INTEGER            :: ReasonEOF

      INTEGER            :: NUM, I, J, FOUND, V1, V2, V3, ELEM_TYPE, NUMELEMS
      INTEGER, DIMENSION(3,2) :: IND
      REAL(KIND=8)       :: X1, X2, X3, Y1, Y2, Y3, LEN, RAD
      REAL(KIND=8), DIMENSION(3) :: XYZ, A, B, C

      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: TEMP_CELL_NEIGHBORS

      INTEGER, DIMENSION(2) :: VLIST2, WHICH1, WHICH2
      INTEGER, DIMENSION(3) :: VLIST3

      INTEGER, DIMENSION(:), ALLOCATABLE :: N_CELLS_WITH_NODE, CELL_WITH_NODE, IOF
      INTEGER :: IDX, JN, JC1, JC2

      LOGICAL, DIMENSION(:), ALLOCATABLE :: NODE_ON_BOUNDARY
      INTEGER :: NUM_BOUNDARY_NODES, NUM_BOUNDARY_ELEM

      REAL(KIND=8) :: DX1, DX2, DX3, DY1, DY2, DY3, XP, YP, W1, W2, W3

      ! Open input file for reading
      OPEN(UNIT=in5,FILE='testgrid.su2', STATUS='old',IOSTAT=ios)

      IF (ios .NE. 0) THEN
         CALL ERROR_ABORT('Attention, mesh file not found! ABORTING.')
      ENDIF

      ! Read the mesh file. SU2 file format (*.su2)

      WRITE(*,*) '==========================================='
      WRITE(*,*) 'Reading grid file in SU2 format.'
      WRITE(*,*) '==========================================='

      
      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE, NUM
         IF (ReasonEOF < 0) EXIT 
         !WRITE(*,*) 'Read line:', LINE, ' number ', NUM
         
         IF (LINE == 'NPOIN=') THEN
            ALLOCATE(U2D_GRID%NODE_COORDS(3,NUM))
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) XYZ(1:2)
               XYZ(3) = 0.d0 ! Stay in the x-y plane.
               U2D_GRID%NODE_COORDS(:,I) = XYZ
            END DO
            U2D_GRID%NUM_NODES = NUM
         ELSE IF (LINE == 'NELEM=') THEN
            ALLOCATE(U2D_GRID%CELL_NODES(3,NUM))

            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) ELEM_TYPE, U2D_GRID%CELL_NODES(:,I)
               !WRITE(*,*) 'I read element ', I, ' has nodes ', U2D_GRID%CELL_NODES(:,I)
               IF (ELEM_TYPE .NE. 5) WRITE(*,*) 'Element type was not 5!'
            END DO
            U2D_GRID%CELL_NODES = U2D_GRID%CELL_NODES + 1 ! Start indexing from 1.

            U2D_GRID%NUM_CELLS = NUM

            ALLOCATE(U2D_GRID%CELL_EDGES_PG(3, U2D_GRID%NUM_CELLS))
            U2D_GRID%CELL_EDGES_PG = -1

            ALLOCATE(U2D_GRID%CELL_PG(U2D_GRID%NUM_CELLS))
            U2D_GRID%CELL_PG = -1
      
         ELSE IF (LINE == 'NMARK=') THEN

            ! Assign physical groups to cell edges.
            DO I = 1, NUM
               
               READ(in5,*, IOSTAT=ReasonEOF) LINE, GROUPNAME
               IF (LINE .NE. 'MARKER_TAG=') THEN
                  WRITE(*,*) 'Error! did not find marker name.'
               ELSE
                  !WRITE(*,*) 'Found marker tag, with groupname: ', GROUPNAME
               END IF
         
               
               READ(in5,*, IOSTAT=ReasonEOF) LINE, NUMELEMS
               IF (LINE .NE. 'MARKER_ELEMS=') THEN
                  WRITE(*,*) 'Error! did not find marker elements.'
               ELSE
                  !WRITE(*,*) 'Found marker elements, with number of elements: ', NUMELEMS
               END IF

               DO J = 1, NUMELEMS
                  READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
               END DO
            END DO

         END IF
      END DO

      REWIND(in5)

      ALLOCATE(N_CELLS_WITH_NODE(U2D_GRID%NUM_NODES))
      ALLOCATE(IOF(U2D_GRID%NUM_NODES))

      N_CELLS_WITH_NODE = 0
      DO I = 1, U2D_GRID%NUM_CELLS
         DO V1 = 1, 3
            JN = U2D_GRID%CELL_NODES(V1,I)
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO
   
      IOF = -1
      IDX = 1
      DO JN = 1, U2D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) .NE. 0) THEN
            IOF(JN) = IDX
            IDX = IDX + N_CELLS_WITH_NODE(JN)
         END IF
      END DO
   
      ALLOCATE(CELL_WITH_NODE(IDX))
      
      N_CELLS_WITH_NODE = 0
      DO I = 1, U2D_GRID%NUM_CELLS
         DO V1 = 1, 3
            JN = U2D_GRID%CELL_NODES(V1,I)
            CELL_WITH_NODE(IOF(JN) + N_CELLS_WITH_NODE(JN)) = I
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO


      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE, NUM
         IF (ReasonEOF < 0) EXIT 
         !WRITE(*,*) 'Read line:', LINE, ' number ', NUM
         
         IF (LINE == 'NPOIN=') THEN
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
            END DO
         ELSE IF (LINE == 'NELEM=') THEN
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
            END DO
         ELSE IF (LINE == 'NMARK=') THEN

            ! Assign physical groups to cells/cell edges.

            ALLOCATE(GRID_BC(NUM)) ! Append the physical group to the list
            N_GRID_BC = NUM

            DO I = 1, NUM
               
               READ(in5,*, IOSTAT=ReasonEOF) LINE, GROUPNAME
      
               GRID_BC(I)%PHYSICAL_GROUP_NAME = GROUPNAME
         
               READ(in5,*, IOSTAT=ReasonEOF) LINE, NUMELEMS

               DO J = 1, NUMELEMS
                  READ(in5,'(A)', IOSTAT=ReasonEOF) LINE

                  READ(LINE,*) ELEM_TYPE

                  IF (ELEM_TYPE == 3) THEN ! element in physical group is a cell edge (segment).

                     READ(LINE,*) ELEM_TYPE, VLIST2
                     
                     VLIST2 = VLIST2 + 1

                     JN = VLIST2(1)
                     IF (N_CELLS_WITH_NODE(JN) > 0) THEN
                        DO IDX = 0, N_CELLS_WITH_NODE(JN) - 1
                           JC1 = CELL_WITH_NODE(IOF(JN) + IDX)
                           FOUND = 0
                           DO V1 = 1, 3
                              IF (ANY(VLIST2 == U2D_GRID%CELL_NODES(V1,JC1))) THEN
                                 FOUND = FOUND + 1
                                 WHICH1(FOUND) = V1
                              END IF
                           END DO
            
                           IF (FOUND == 2) THEN
                              IF (ANY(WHICH1 == 1) .AND. ANY(WHICH1 == 2)) THEN
                                 U2D_GRID%CELL_EDGES_PG(1, JC1) = I
                              ELSE IF (ANY(WHICH1 == 2) .AND. ANY(WHICH1 == 3)) THEN
                                 U2D_GRID%CELL_EDGES_PG(2, JC1) = I
                              ELSE IF (ANY(WHICH1 == 3) .AND. ANY(WHICH1 == 1)) THEN
                                 U2D_GRID%CELL_EDGES_PG(3, JC1) = I
                              END IF
                           END IF
                        END DO
                     END IF

                  ELSE IF (ELEM_TYPE == 5) THEN ! element in physical group is a cell (simplex).
                     READ(LINE,*) ELEM_TYPE, VLIST3
                     
                     VLIST3 = VLIST3 + 1

                     JN = VLIST3(1)
                     IF (N_CELLS_WITH_NODE(JN) > 0) THEN
                        DO IDX = 0, N_CELLS_WITH_NODE(JN) - 1
                           JC1 = CELL_WITH_NODE(IOF(JN) + IDX)
                           FOUND = 0
                           DO V1 = 1, 3
                              IF (ANY(VLIST3 == U2D_GRID%CELL_NODES(V1,JC1))) THEN
                                 FOUND = FOUND + 1
                              END IF
                           END DO
            
                           IF (FOUND == 3) THEN
                              U2D_GRID%CELL_PG(JC1) = I
                           END IF
                        END DO
                     END IF

                  ELSE
                     WRITE(*,*) 'Error! element type was not line or triangle.'
                  END IF

               END DO
            END DO

         END IF
      END DO

      ! Done reading
      CLOSE(in5)

      !WRITE(*,*) 'Read grid file. It contains ', U2D_GRID%NUM_NODES, &
      !           'points, and ', U2D_GRID%NUM_CELLS, 'cells.'

      ! Process the mesh: generate connectivity, normals and such...
      !XMIN, XMAX,...

      !DO I = 1, U2D_GRID%NUM_CELLS
      !   WRITE(*,*) U2D_GRID%CELL_NODES(:,I)
      !END DO


      WRITE(*,*) '==========================================='
      WRITE(*,*) 'Computing cell volumes.'
      WRITE(*,*) '==========================================='

      ! Compute cell volumes
      ALLOCATE(U2D_GRID%CELL_AREAS(U2D_GRID%NUM_CELLS))
      ALLOCATE(U2D_GRID%CELL_VOLUMES(U2D_GRID%NUM_CELLS))
      DO I = 1, U2D_GRID%NUM_CELLS
         A = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(1,I))
         B = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(2,I))
         C = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(3,I))

         U2D_GRID%CELL_AREAS(I) = 0.5*ABS(A(1)*(B(2)-C(2)) + B(1)*(C(2)-A(2)) + C(1)*(A(2)-B(2)))
         IF (DIMS == 2 .AND. .NOT. AXI) THEN
            U2D_GRID%CELL_VOLUMES(I) = U2D_GRID%CELL_AREAS(I) * 1.d0
            !WRITE(*,*) U2D_GRID%CELL_VOLUMES(I)
         END IF
         IF (DIMS == 2 .AND. AXI) THEN
            RAD = (A(2)+B(2)+C(2))/3.
            U2D_GRID%CELL_VOLUMES(I) = U2D_GRID%CELL_AREAS(I) * 2.*PI*RAD
         END IF
      END DO


      WRITE(*,*) '==========================================='
      WRITE(*,*) 'Computing cell centroids.'
      WRITE(*,*) '==========================================='

      ! Compute cell volumes
      ALLOCATE(U2D_GRID%CELL_CENTROIDS(3, U2D_GRID%NUM_CELLS))
      DO I = 1, U2D_GRID%NUM_CELLS
         A = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(1,I))
         B = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(2,I))
         C = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(3,I))

         U2D_GRID%CELL_CENTROIDS(:, I) = (A+B+C)/3.
      END DO


      WRITE(*,*) '==========================================='
      WRITE(*,*) 'Computing grid connectivity.'
      WRITE(*,*) '==========================================='


      ! Find cell connectivity
      ALLOCATE(TEMP_CELL_NEIGHBORS(3, U2D_GRID%NUM_CELLS))
      TEMP_CELL_NEIGHBORS = -1



      DO JN = 1, U2D_GRID%NUM_NODES
         !WRITE(*,*) 'Checking node ', JN, ' of ',  U2D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) > 1) THEN
            DO I = 0, N_CELLS_WITH_NODE(JN) - 1
               DO J = I, N_CELLS_WITH_NODE(JN) - 1
                  IF (I == J) CYCLE
                  JC1 = CELL_WITH_NODE(IOF(JN) + I)
                  JC2 = CELL_WITH_NODE(IOF(JN) + J)


                  FOUND = 0
                  DO V1 = 1, 3
                     DO V2 = 1, 3
                        IF (U2D_GRID%CELL_NODES(V1,JC1) == U2D_GRID%CELL_NODES(V2,JC2)) THEN
                           FOUND = FOUND + 1
                           IF (FOUND .GT. 2) CALL ERROR_ABORT('Error! Found duplicate cells in the mesh!')
                           WHICH1(FOUND) = V1
                           WHICH2(FOUND) = V2
                        END IF
                     END DO
                  END DO

                  IF (FOUND == 2) THEN
      
                     IF (ANY(WHICH1 == 1) .AND. ANY(WHICH1 == 2)) THEN
                        TEMP_CELL_NEIGHBORS(1, JC1) = JC2
                     ELSE IF (ANY(WHICH1 == 2) .AND. ANY(WHICH1 == 3)) THEN
                        TEMP_CELL_NEIGHBORS(2, JC1) = JC2
                     ELSE IF (ANY(WHICH1 == 3) .AND. ANY(WHICH1 == 1)) THEN
                        TEMP_CELL_NEIGHBORS(3, JC1) = JC2
                     END IF

                     IF (ANY(WHICH2 == 1) .AND. ANY(WHICH2 == 2)) THEN
                        TEMP_CELL_NEIGHBORS(1, JC2) = JC1
                     ELSE IF (ANY(WHICH2 == 2) .AND. ANY(WHICH2 == 3)) THEN
                        TEMP_CELL_NEIGHBORS(2, JC2) = JC1
                     ELSE IF (ANY(WHICH2 == 3) .AND. ANY(WHICH2 == 1)) THEN
                        TEMP_CELL_NEIGHBORS(3, JC2) = JC1
                     END IF
      
                  END IF


               END DO
            END DO
         END IF
      END DO

      U2D_GRID%CELL_NEIGHBORS = TEMP_CELL_NEIGHBORS



      !WRITE(*,*) 'Generated grid connectivity. '
      !DO I = 1, U2D_GRID%NUM_CELLS
      !   WRITE(*,*) 'Cell ', I, ' neighbors cells ', TEMP_CELL_NEIGHBORS(:, I)
      !END DO


      WRITE(*,*) '==========================================='
      WRITE(*,*) 'Computing face normals.'
      WRITE(*,*) '==========================================='


      ! Compute cell edge normals
      IND(1,:) = [1,2]
      IND(2,:) = [2,3]
      IND(3,:) = [3,1]
      ALLOCATE(U2D_GRID%EDGE_NORMAL(3, 3, U2D_GRID%NUM_CELLS))
      ALLOCATE(U2D_GRID%CELL_EDGES_LEN(3, U2D_GRID%NUM_CELLS))
      ALLOCATE(U2D_GRID%CELL_FACES_AREA(3, U2D_GRID%NUM_CELLS))

      DO I = 1, U2D_GRID%NUM_CELLS
         DO J = 1, 3
            X1 = U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(IND(J,1),I))
            X2 = U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(IND(J,2),I))
            Y1 = U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(IND(J,1),I))
            Y2 = U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(IND(J,2),I))
            LEN = SQRT((Y2-Y1)*(Y2-Y1) + (X2-X1)*(X2-X1))
            U2D_GRID%CELL_EDGES_LEN(J,I) = LEN
            U2D_GRID%EDGE_NORMAL(1,J,I) = (Y2-Y1)/LEN
            U2D_GRID%EDGE_NORMAL(2,J,I) = (X1-X2)/LEN
            U2D_GRID%EDGE_NORMAL(3,J,I) = 0.d0
            IF (AXI) THEN
               U2D_GRID%CELL_FACES_AREA(J,I) = PI*LEN*(Y1+Y2)
            ELSE
               U2D_GRID%CELL_FACES_AREA(J,I) = LEN
            END IF
         END DO
      END DO


      WRITE(*,*) '==========================================='
      WRITE(*,*) 'Checking ordering.'
      WRITE(*,*) '==========================================='


      DO I = 1, U2D_GRID%NUM_CELLS
         X1 = U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(2,I)) &
            - U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(1,I))
         X2 = U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(3,I)) &
            - U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(1,I))
         Y1 = U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(2,I)) &
            - U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(1,I))
         Y2 = U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(3,I)) &
            - U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(1,I))

         IF (X1*Y2-X2*Y1 < 0) CALL ERROR_ABORT('2D mesh triangles have negative z-normal.')

      END DO

      NCELLS = U2D_GRID%NUM_CELLS
      NNODES = U2D_GRID%NUM_NODES



      ALLOCATE(U2D_GRID%BASIS_COEFFS(3,3,NCELLS))

      DO I = 1, NCELLS
         V1 = U2D_GRID%CELL_NODES(1,I)
         V2 = U2D_GRID%CELL_NODES(2,I)
         V3 = U2D_GRID%CELL_NODES(3,I)

         X1 = U2D_GRID%NODE_COORDS(1, V1)
         X2 = U2D_GRID%NODE_COORDS(1, V2)
         X3 = U2D_GRID%NODE_COORDS(1, V3)
         Y1 = U2D_GRID%NODE_COORDS(2, V1)
         Y2 = U2D_GRID%NODE_COORDS(2, V2)
         Y3 = U2D_GRID%NODE_COORDS(2, V3)

         ! These are such that PSI_i = SUM_j [ x_j * BASIS_COEFFS(j,i,IC) ] + BASIS_COEFFS(3,i,IC)

         U2D_GRID%BASIS_COEFFS(1,1,I) =  Y2-Y3
         U2D_GRID%BASIS_COEFFS(2,1,I) = -(X2-X3)
         U2D_GRID%BASIS_COEFFS(3,1,I) =  X2*Y3 - X3*Y2

         U2D_GRID%BASIS_COEFFS(1,2,I) = -(Y1-Y3)
         U2D_GRID%BASIS_COEFFS(2,2,I) =  X1-X3
         U2D_GRID%BASIS_COEFFS(3,2,I) =  X3*Y1 - X1*Y3

         U2D_GRID%BASIS_COEFFS(1,3,I) = -(Y2-Y1)
         U2D_GRID%BASIS_COEFFS(2,3,I) =  X2-X1
         U2D_GRID%BASIS_COEFFS(3,3,I) =  X1*Y2 - X2*Y1

         U2D_GRID%BASIS_COEFFS(:,:,I) = 0.5*U2D_GRID%BASIS_COEFFS(:,:,I)/U2D_GRID%CELL_AREAS(I)

      END DO


      WRITE(*,*) '==========================================='
      WRITE(*,*) 'Computing least squares matrices.'
      WRITE(*,*) '==========================================='


      ALLOCATE(U2D_GRID%LSTSQ_COEFFS(2,3,NCELLS))

      DO I = 1, NCELLS
         V1 = U2D_GRID%CELL_NEIGHBORS(1,I)
         V2 = U2D_GRID%CELL_NEIGHBORS(2,I)
         V3 = U2D_GRID%CELL_NEIGHBORS(3,I)

         XP = U2D_GRID%CELL_CENTROIDS(1, I)
         YP = U2D_GRID%CELL_CENTROIDS(2, I)

         IF (V1 == -1) THEN
            W1 = 0.d0
            DX1 = 1.d6
            DY1 = 1.d6
         ELSE
            DX1 = U2D_GRID%CELL_CENTROIDS(1, V1) - XP
            DY1 = U2D_GRID%CELL_CENTROIDS(2, V1) - YP
            W1 = 1./SQRT(DX1*DX1 + DY1*DY1)
         END IF

         IF (V2 == -1) THEN
            W2 = 0.d0
            DX2 = 1.d6
            DY2 = 1.d6
         ELSE
            DX2 = U2D_GRID%CELL_CENTROIDS(1, V2) - XP
            DY2 = U2D_GRID%CELL_CENTROIDS(2, V2) - YP
            W2 = 1./SQRT(DX2*DX2 + DY2*DY2)
         END IF

         IF (V3 == -1) THEN
            W3 = 0.d0
            DX3 = 1.d6
            DY3 = 1.d6
         ELSE
            DX3 = U2D_GRID%CELL_CENTROIDS(1, V3) - XP
            DY3 = U2D_GRID%CELL_CENTROIDS(2, V3) - YP
            W3 = 1./SQRT(DX3*DX3 + DY3*DY3)
         END IF

         U2D_GRID%LSTSQ_COEFFS(1,1,I) = dx3*w3*(-dx1*dy1*w1**2 - dx2*dy2*w2**2 - dx3*dy3*w3**2)/(dx1**2*dy2**2*w1**2*w2**2 + & 
         dx1**2*dy3**2*w1**2*w3**2 - 2*dx1*dx2*dy1*dy2*w1**2*w2**2 - 2*dx1*dx3*dy1*dy3*w1**2*w3**2 + dx2**2*dy1**2*w1**2*w2**2 + &
         dx2**2*dy3**2*w2**2*w3**2 - 2*dx2*dx3*dy2*dy3*w2**2*w3**2 + dx3**2*dy1**2*w1**2*w3**2 + dx3**2*dy2**2*w2**2*w3**2) + &
         dy3*w3*(dx1**2*w1**2 + dx2**2*w2**2 + dx3**2*w3**2)/(dx1**2*dy2**2*w1**2*w2**2 + dx1**2*dy3**2*w1**2*w3**2 - &
         2*dx1*dx2*dy1*dy2*w1**2*w2**2 - 2*dx1*dx3*dy1*dy3*w1**2*w3**2 + dx2**2*dy1**2*w1**2*w2**2 + dx2**2*dy3**2*w2**2*w3**2 - &
         2*dx2*dx3*dy2*dy3*w2**2*w3**2 + dx3**2*dy1**2*w1**2*w3**2 + dx3**2*dy2**2*w2**2*w3**2)

         U2D_GRID%LSTSQ_COEFFS(1,2,I) = dx2*w2*(dy1**2*w1**2 + dy2**2*w2**2 + dy3**2*w3**2)/(dx1**2*dy2**2*w1**2*w2**2 + &
         dx1**2*dy3**2*w1**2*w3**2 - 2*dx1*dx2*dy1*dy2*w1**2*w2**2 - 2*dx1*dx3*dy1*dy3*w1**2*w3**2 + dx2**2*dy1**2*w1**2*w2**2 + &
         dx2**2*dy3**2*w2**2*w3**2 - 2*dx2*dx3*dy2*dy3*w2**2*w3**2 + dx3**2*dy1**2*w1**2*w3**2 + dx3**2*dy2**2*w2**2*w3**2) + &
         dy2*w2*(-dx1*dy1*w1**2 - dx2*dy2*w2**2 - dx3*dy3*w3**2)/(dx1**2*dy2**2*w1**2*w2**2 + dx1**2*dy3**2*w1**2*w3**2 - &
         2*dx1*dx2*dy1*dy2*w1**2*w2**2 - 2*dx1*dx3*dy1*dy3*w1**2*w3**2 + dx2**2*dy1**2*w1**2*w2**2 + dx2**2*dy3**2*w2**2*w3**2 - &
         2*dx2*dx3*dy2*dy3*w2**2*w3**2 + dx3**2*dy1**2*w1**2*w3**2 + dx3**2*dy2**2*w2**2*w3**2)

         U2D_GRID%LSTSQ_COEFFS(1,3,I) = dx3*w3*(dy1**2*w1**2 + dy2**2*w2**2 + dy3**2*w3**2)/(dx1**2*dy2**2*w1**2*w2**2 + &
         dx1**2*dy3**2*w1**2*w3**2 - 2*dx1*dx2*dy1*dy2*w1**2*w2**2 - 2*dx1*dx3*dy1*dy3*w1**2*w3**2 + dx2**2*dy1**2*w1**2*w2**2 + &
         dx2**2*dy3**2*w2**2*w3**2 - 2*dx2*dx3*dy2*dy3*w2**2*w3**2 + dx3**2*dy1**2*w1**2*w3**2 + dx3**2*dy2**2*w2**2*w3**2) + &
         dy3*w3*(-dx1*dy1*w1**2 - dx2*dy2*w2**2 - dx3*dy3*w3**2)/(dx1**2*dy2**2*w1**2*w2**2 + dx1**2*dy3**2*w1**2*w3**2 - &
         2*dx1*dx2*dy1*dy2*w1**2*w2**2 - 2*dx1*dx3*dy1*dy3*w1**2*w3**2 + dx2**2*dy1**2*w1**2*w2**2 + dx2**2*dy3**2*w2**2*w3**2 - &
         2*dx2*dx3*dy2*dy3*w2**2*w3**2 + dx3**2*dy1**2*w1**2*w3**2 + dx3**2*dy2**2*w2**2*w3**2)

         U2D_GRID%LSTSQ_COEFFS(2,1,I) = dx1*w1*(-dx1*dy1*w1**2 - dx2*dy2*w2**2 - dx3*dy3*w3**2)/(dx1**2*dy2**2*w1**2*w2**2 + &
         dx1**2*dy3**2*w1**2*w3**2 - 2*dx1*dx2*dy1*dy2*w1**2*w2**2 - 2*dx1*dx3*dy1*dy3*w1**2*w3**2 + dx2**2*dy1**2*w1**2*w2**2 + &
         dx2**2*dy3**2*w2**2*w3**2 - 2*dx2*dx3*dy2*dy3*w2**2*w3**2 + dx3**2*dy1**2*w1**2*w3**2 + dx3**2*dy2**2*w2**2*w3**2) + &
         dy1*w1*(dx1**2*w1**2 + dx2**2*w2**2 + dx3**2*w3**2)/(dx1**2*dy2**2*w1**2*w2**2 + dx1**2*dy3**2*w1**2*w3**2 - &
         2*dx1*dx2*dy1*dy2*w1**2*w2**2 - 2*dx1*dx3*dy1*dy3*w1**2*w3**2 + dx2**2*dy1**2*w1**2*w2**2 + dx2**2*dy3**2*w2**2*w3**2 - &
         2*dx2*dx3*dy2*dy3*w2**2*w3**2 + dx3**2*dy1**2*w1**2*w3**2 + dx3**2*dy2**2*w2**2*w3**2)

         U2D_GRID%LSTSQ_COEFFS(2,2,I) = dx2*w2*(-dx1*dy1*w1**2 - dx2*dy2*w2**2 - dx3*dy3*w3**2)/(dx1**2*dy2**2*w1**2*w2**2 + &
         dx1**2*dy3**2*w1**2*w3**2 - 2*dx1*dx2*dy1*dy2*w1**2*w2**2 - 2*dx1*dx3*dy1*dy3*w1**2*w3**2 + dx2**2*dy1**2*w1**2*w2**2 + &
         dx2**2*dy3**2*w2**2*w3**2 - 2*dx2*dx3*dy2*dy3*w2**2*w3**2 + dx3**2*dy1**2*w1**2*w3**2 + dx3**2*dy2**2*w2**2*w3**2) + &
         dy2*w2*(dx1**2*w1**2 + dx2**2*w2**2 + dx3**2*w3**2)/(dx1**2*dy2**2*w1**2*w2**2 + dx1**2*dy3**2*w1**2*w3**2 - &
         2*dx1*dx2*dy1*dy2*w1**2*w2**2 - 2*dx1*dx3*dy1*dy3*w1**2*w3**2 + dx2**2*dy1**2*w1**2*w2**2 + dx2**2*dy3**2*w2**2*w3**2 - &
         2*dx2*dx3*dy2*dy3*w2**2*w3**2 + dx3**2*dy1**2*w1**2*w3**2 + dx3**2*dy2**2*w2**2*w3**2)

         U2D_GRID%LSTSQ_COEFFS(2,2,I) = dx3*w3*(-dx1*dy1*w1**2 - dx2*dy2*w2**2 - dx3*dy3*w3**2)/(dx1**2*dy2**2*w1**2*w2**2 + &
         dx1**2*dy3**2*w1**2*w3**2 - 2*dx1*dx2*dy1*dy2*w1**2*w2**2 - 2*dx1*dx3*dy1*dy3*w1**2*w3**2 + dx2**2*dy1**2*w1**2*w2**2 + &
         dx2**2*dy3**2*w2**2*w3**2 - 2*dx2*dx3*dy2*dy3*w2**2*w3**2 + dx3**2*dy1**2*w1**2*w3**2 + dx3**2*dy2**2*w2**2*w3**2) + &
         dy3*w3*(dx1**2*w1**2 + dx2**2*w2**2 + dx3**2*w3**2)/(dx1**2*dy2**2*w1**2*w2**2 + dx1**2*dy3**2*w1**2*w3**2 - &
         2*dx1*dx2*dy1*dy2*w1**2*w2**2 - 2*dx1*dx3*dy1*dy3*w1**2*w3**2 + dx2**2*dy1**2*w1**2*w2**2 + dx2**2*dy3**2*w2**2*w3**2 - &
         2*dx2*dx3*dy2*dy3*w2**2*w3**2 + dx3**2*dy1**2*w1**2*w3**2 + dx3**2*dy2**2*w2**2*w3**2)


      END DO

      WRITE(*,*) '==========================================='
      WRITE(*,*) 'Creating boundary grid.'
      WRITE(*,*) '==========================================='


      ALLOCATE(U2D_GRID%CELL_EDGES_BOUNDARY_INDEX(3,NCELLS))
      U2D_GRID%CELL_EDGES_BOUNDARY_INDEX = -1
      ALLOCATE(NODE_ON_BOUNDARY(NNODES))
      NODE_ON_BOUNDARY = .FALSE.
      ALLOCATE(U2D_GRID%NODES_BOUNDARY_INDEX(NNODES))
      U2D_GRID%NODES_BOUNDARY_INDEX = -1
      NUM_BOUNDARY_NODES = 0
      NUM_BOUNDARY_ELEM = 0
      DO I = 1, NCELLS
         DO J = 1, 3
            ! If the edge belongs to any physical group, it should be part of the boundary grid
            ! Later, we may want to filter this further
            IF (U2D_GRID%CELL_EDGES_PG(J,I) .NE. -1) THEN
               NUM_BOUNDARY_ELEM = NUM_BOUNDARY_ELEM + 1
               V1 = U2D_GRID%CELL_NODES(J, I)
               IF (J == 3) THEN
                  V2 = U2D_GRID%CELL_NODES(1, I)
               ELSE
                  V2 = U2D_GRID%CELL_NODES(J+1, I)
               END IF
               IF (.NOT. NODE_ON_BOUNDARY(V1)) THEN
                  NUM_BOUNDARY_NODES = NUM_BOUNDARY_NODES + 1
                  U2D_GRID%NODES_BOUNDARY_INDEX(V1) = NUM_BOUNDARY_NODES
                  NODE_ON_BOUNDARY(V1) = .TRUE.
               END IF
               IF (.NOT. NODE_ON_BOUNDARY(V2)) THEN
                  NUM_BOUNDARY_NODES = NUM_BOUNDARY_NODES + 1
                  U2D_GRID%NODES_BOUNDARY_INDEX(V2) = NUM_BOUNDARY_NODES
                  NODE_ON_BOUNDARY(V2) = .TRUE.
               END IF
               
            END IF
         END DO
      END DO      

      U1D_GRID%NUM_CELLS = NUM_BOUNDARY_ELEM
      U1D_GRID%NUM_NODES = NUM_BOUNDARY_NODES
      ALLOCATE(U1D_GRID%CELL_NODES(2, NUM_BOUNDARY_ELEM))
      ALLOCATE(U1D_GRID%CELL_PG(NUM_BOUNDARY_ELEM))
      ALLOCATE(U1D_GRID%NODE_COORDS(3, NUM_BOUNDARY_NODES))

      DO I = 1, NNODES
         IF (NODE_ON_BOUNDARY(I)) THEN
            U1D_GRID%NODE_COORDS(:,U2D_GRID%NODES_BOUNDARY_INDEX(I)) = U2D_GRID%NODE_COORDS(:,I)
         END IF
      END DO

      NUM_BOUNDARY_ELEM = 0

      DO I = 1, NCELLS
         DO J = 1, 3
            IF (U2D_GRID%CELL_EDGES_PG(J,I) .NE. -1) THEN
               NUM_BOUNDARY_ELEM = NUM_BOUNDARY_ELEM + 1
               U1D_GRID%CELL_PG(NUM_BOUNDARY_ELEM) = U2D_GRID%CELL_EDGES_PG(J,I)
               U2D_GRID%CELL_EDGES_BOUNDARY_INDEX(J,I) = NUM_BOUNDARY_ELEM

               V1 = U2D_GRID%CELL_NODES(J, I)
               IF (J == 3) THEN
                  V2 = U2D_GRID%CELL_NODES(1, I)
               ELSE
                  V2 = U2D_GRID%CELL_NODES(J+1, I)
               END IF
               U1D_GRID%CELL_NODES(1, NUM_BOUNDARY_ELEM) = U2D_GRID%NODES_BOUNDARY_INDEX(V1)
               U1D_GRID%CELL_NODES(2, NUM_BOUNDARY_ELEM) = U2D_GRID%NODES_BOUNDARY_INDEX(V2)

            END IF
         END DO
      END DO
      
      DEALLOCATE(NODE_ON_BOUNDARY)

      NBOUNDCELLS = NUM_BOUNDARY_ELEM
      NBOUNDNODES = NUM_BOUNDARY_NODES

      ! Compute areas and lengths of boundary mesh
      ALLOCATE(U1D_GRID%SEGMENT_LENGTHS(U1D_GRID%NUM_CELLS))
      ALLOCATE(U1D_GRID%SEGMENT_AREAS(U1D_GRID%NUM_CELLS))
      DO I = 1, U1D_GRID%NUM_CELLS
         A = U1D_GRID%NODE_COORDS(:, U1D_GRID%CELL_NODES(1,I))
         B = U1D_GRID%NODE_COORDS(:, U1D_GRID%CELL_NODES(2,I))

         U1D_GRID%SEGMENT_LENGTHS(I) = SQRT((B(1) - A(1))**2 + (B(2) - A(2))**2)
         IF (.NOT. AXI) THEN
            U1D_GRID%SEGMENT_AREAS(I) = U1D_GRID%SEGMENT_LENGTHS(I) * 1.d0
         ELSE
            RAD = 0.5*(A(2)+B(2))
            U1D_GRID%SEGMENT_AREAS(I) = U1D_GRID%SEGMENT_LENGTHS(I) * 2.*PI*RAD
         END IF
      END DO


      WRITE(*,*) '============================================================='
      WRITE(*,*) 'Done reading grid file.'
      WRITE(*,*) 'It contains ', NNODES, ' nodes and ', NCELLS, ' cells.'
      WRITE(*,*) 'The boundary grid contains ', NBOUNDCELLS, ' lines and ', NBOUNDNODES, ' nodes.'
      WRITE(*,*) '============================================================='


   END SUBROUTINE READ_2D_UNSTRUCTURED_GRID_SU2


   SUBROUTINE ERROR_ABORT(string)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)  :: string

      PRINT*
      WRITE(*,*) '   ', string
      PRINT*
      STOP
      
   END SUBROUTINE ERROR_ABORT





   SUBROUTINE SPLIT_STR(STRING, DELIMITER, STRARRAY, N_STR)
      ! splitstring splits a string to an array of
      ! substrings based on a selected delimiter
      ! note any facing space/blank in substrings will be removed

      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: N_STR
      CHARACTER(LEN=80), ALLOCATABLE, INTENT(OUT) :: STRARRAY(:)
      CHARACTER(LEN=*), INTENT(IN) :: STRING
      INTEGER   :: n, i, j, idx
      CHARACTER(len=80) :: STRTMP = ''
      CHARACTER, INTENT(IN) :: DELIMITER

      n=LEN(STRING)
      ALLOCATE(STRARRAY(n))
      j = 1
      idx = 0
      DO i=1, n
         
         IF (STRING(i:i) /= DELIMITER) THEN
            STRTMP(j:j) = STRING(i:i)
            j = j + 1
            IF (i==n .OR. STRING(i+1:i+1) == DELIMITER) THEN
               j = 1
               idx = idx + 1
               STRARRAY(idx) = STRTMP
               STRTMP = ''
            END IF
         END IF
      END DO

      N_STR = idx

   END SUBROUTINE SPLIT_STR

end module
