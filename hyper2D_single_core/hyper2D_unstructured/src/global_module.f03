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

   LOGICAL                                   :: BOOL_DUMP_MOMENTS = .FALSE.
   CHARACTER(LEN=256) :: GRID_FILENAME
   LOGICAL :: BOOL_BINARY_OUTPUT = .TRUE.
   CHARACTER*256                           :: FLOWFIELD_SAVE_PATH = 'dumps/'
   ! Time integration

   real(kind=8), parameter :: t_end      = 50.0d-0 ! [s] total simulated time (from 0 to t_end)
   real(kind=8), parameter :: CFL_target = 0.1
   real(kind=8), parameter :: dtmax      = 1.d-2


   
   integer, parameter :: Neq = 4 ! Number of equations (per-species)

   !!!! TEST TEST TEST !!! ! Reconstruction order 
   !!!! TEST TEST TEST !!! ! integer, parameter :: reconstr_order = 0 ! 0: no reconstruction -> first order in space
   !!!! TEST TEST TEST !!! integer, parameter :: reconstr_order = 1 ! 1: linear reconstruction -> second order in space

   ! Free stream
   real(kind=8), parameter :: rho0 = 254.0*1.66e-27 * 4.0/(1.380649d-23*300.0)   ! [kg/m3]
   real(kind=8), parameter :: ux0  = 0.0    ! [m/s]
   real(kind=8), parameter :: uy0  = 0.0  ! [m/s]
   real(kind=8), parameter :: T0   = 300.0   ! [K]
   real(kind=8), parameter :: Tw   = 300.0   ! [K]

   ! Utilities
   real(kind=8) :: ws_over_sqrtA_maxabs




   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Multispecies !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   INTEGER            :: NSPECIES = 0
   
   TYPE SPECIES_DATA_STRUCTURE
      CHARACTER*10 :: NAME
      REAL(KIND=8) :: MOLECULAR_MASS
      REAL(KIND=8) :: DIAM
      REAL(KIND=8) :: GAMMA
      REAL(KIND=8) :: CP
      REAL(KIND=8) :: KAPPA
      REAL(KIND=8) :: MU
   END TYPE SPECIES_DATA_STRUCTURE

   TYPE(SPECIES_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: SPECIES



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Wall reactions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   CHARACTER(LEN=256) :: WALL_REACTIONS_FILENAME
   INTEGER :: N_WALL_REACTIONS = 0
   
   TYPE WALL_REACTIONS_DATA_STRUCTURE
      INTEGER :: R_SP_ID
      INTEGER :: P1_SP_ID
      INTEGER :: P2_SP_ID
      REAL(KIND=8) :: PROB
      INTEGER :: N_PROD
   END TYPE WALL_REACTIONS_DATA_STRUCTURE

   TYPE(WALL_REACTIONS_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: WALL_REACTIONS


end module
