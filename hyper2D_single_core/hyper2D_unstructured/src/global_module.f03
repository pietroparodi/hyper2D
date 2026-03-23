module global_module

   implicit none

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Constants !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   REAL(KIND=8) :: PI   = 3.1415926535897932d0           ! https://dlmf.nist.gov/3.12
   REAL(KIND=8) :: EPS0 = 8.8541878128d-12               ! https://physics.nist.gov/cgi-bin/cuu/Value?ep0
   REAL(KIND=8) :: MU0  = 1.25663706212d-6               ! https://physics.nist.gov/cgi-bin/cuu/Value?mu0
   REAL(KIND=8) :: KB   = 1.380649d-23                   ! https://physics.nist.gov/cgi-bin/cuu/Value?k
   REAL(KIND=8) :: QE   = 1.602176634d-19                ! https://physics.nist.gov/cgi-bin/cuu/Value?e
   REAL(KIND=8) :: NA   = 6.02214076e23                  ! https://physics.nist.gov/cgi-bin/cuu/Value?na
   REAL(KIND=8) :: ME   = 9.1093837139d-31               ! https://physics.nist.gov/cgi-bin/cuu/Value?me
   REAL(KIND=8) :: HP   = 6.62607015d-34                 ! https://physics.nist.gov/cgi-bin/cuu/Value?h

   CHARACTER(LEN=256) :: GRID_FILENAME
   LOGICAL :: BOOL_BINARY_OUTPUT = .TRUE.
   CHARACTER*256                           :: FLOWFIELD_SAVE_PATH = 'dumps/'
   ! Time integration

   real(kind=8) :: t_end      = 0.01d0 ! [s] total simulated time (from 0 to t_end)
   real(kind=8) :: CFL_target = 0.25
   real(kind=8) :: dtmax      = 1.d-2
   INTEGER      :: STATS_EVERY
   INTEGER      :: DUMP_GRID_EVERY
   INTEGER      :: DUMP_GRID_START

   
   integer, parameter :: Neq = 4 ! Number of equations (per-species)

   !!!! TEST TEST TEST !!! ! Reconstruction order 
   !!!! TEST TEST TEST !!! ! integer, parameter :: reconstr_order = 0 ! 0: no reconstruction -> first order in space
   !!!! TEST TEST TEST !!! integer, parameter :: reconstr_order = 1 ! 1: linear reconstruction -> second order in space

   ! Wall temperature
   real(kind=8), parameter :: Tw   = 300.0   ! [K]

   REAL(KIND=8) :: INITIAL_NRHO = 0.d0
   REAL(KIND=8) :: INITIAL_UX = 0.d0
   REAL(KIND=8) :: INITIAL_UY = 0.d0
   REAL(KIND=8) :: INITIAL_UZ = 0.d0
   REAL(KIND=8) :: INITIAL_TEMP = 0.d0
   INTEGER :: INITIAL_MIX_ID


   INTEGER :: FLUID_MIX_ID
   INTEGER :: N_SPECIES_FLUID = 0


   ! Utilities
   real(kind=8) :: invdt_adv, invdt_cond, invdt_diff


   REAL(KIND=8) :: BG_DENS
   REAL(KIND=8) :: BG_TEMP
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: BG_CELL_NRHO
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: BG_CELL_VX
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: BG_CELL_VY
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: BG_CELL_VZ
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: BG_CELL_TEMP
   INTEGER :: BG_MIX
   LOGICAL      :: BOOL_BG_FILE

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: VSS_GREFS ! Matrix of reference relative velocities for VSS
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: VSS_SIGMAS ! Matrix of reference cross sections for VSS
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: VSS_ALPHAS ! Matrix of reference scattering coeff. for VSS
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: VSS_OMEGAS ! Matrix of reference temperature exponent for VSS


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Multispecies !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   INTEGER            :: N_SPECIES = 0
   
   TYPE SPECIES_DATA_STRUCTURE
      CHARACTER*10 :: NAME
      REAL(KIND=8) :: MOLECULAR_MASS
      REAL(KIND=8) :: DIAM
      REAL(KIND=8) :: GAMMA
      REAL(KIND=8) :: CP
      REAL(KIND=8) :: KAPPA
      REAL(KIND=8) :: MU
      REAL(KIND=8) :: OMEGA
      REAL(KIND=8) :: TREF
      REAL(KIND=8) :: ALPHA
      REAL(KIND=8) :: SIGMA
   END TYPE SPECIES_DATA_STRUCTURE

   TYPE(SPECIES_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: SPECIES, TEMP_SPECIES


   INTEGER            :: N_MIXTURES = 0

   TYPE MIXTURE_COMPONENT
      INTEGER :: ID
      CHARACTER*64 :: NAME
      REAL(KIND=8) :: MOLFRAC
   END TYPE MIXTURE_COMPONENT

   TYPE MIXTURE
      CHARACTER*64 :: NAME
      INTEGER      :: N_COMPONENTS
      TYPE(MIXTURE_COMPONENT), DIMENSION(:), ALLOCATABLE :: COMPONENTS
   END TYPE MIXTURE

   TYPE(MIXTURE), DIMENSION(:), ALLOCATABLE :: MIXTURES
 


   CHARACTER(LEN=256) :: WALL_REACTIONS_FILENAME
   INTEGER :: N_WALL_REACTIONS = 0
   
   TYPE WALL_REACTIONS_DATA_STRUCTURE
      INTEGER :: R_SP_ID
      INTEGER :: P1_SP_ID
      INTEGER :: P2_SP_ID
      REAL(KIND=8) :: PROB
      INTEGER :: N_PROD
   END TYPE WALL_REACTIONS_DATA_STRUCTURE

   TYPE(WALL_REACTIONS_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: WALL_REACTIONS, TEMP_WALL_REACTIONS



   CHARACTER(LEN=256) :: REACTIONS_FILENAME
   INTEGER :: N_REACTIONS = 0

   ENUM, BIND(C)
      ENUMERATOR FIXED_RATE, TCE, LXCAT
   END ENUM
   
   TYPE REACTIONS_DATA_STRUCTURE
      INTEGER(KIND(FIXED_RATE)) :: TYPE
      INTEGER :: R1_SP_ID
      INTEGER :: R2_SP_ID
      INTEGER :: P1_SP_ID
      INTEGER :: P2_SP_ID
      INTEGER :: P3_SP_ID
      INTEGER :: P4_SP_ID
      REAL(KIND=8) :: A, N, EA
      REAL(KIND=8) :: C1, C2, C3
      INTEGER :: N_PROD
      LOGICAL :: IS_CEX
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TABLE_ENERGY
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TABLE_CS
      REAL(KIND=8) :: MAX_SIGMA
      REAL(KIND=8) :: CONSTANT_CS
      INTEGER :: COUNTS
   END TYPE REACTIONS_DATA_STRUCTURE

   TYPE(REACTIONS_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: REACTIONS, TEMP_REACTIONS


   ENUM, BIND(C)
      ENUMERATOR AUSM, HLL
   END ENUM

   INTEGER(KIND(AUSM)) :: FLUX_FUNCTION = AUSM
   

end module
