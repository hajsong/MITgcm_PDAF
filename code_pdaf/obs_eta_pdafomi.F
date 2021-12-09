C$Id: obs_TYPE_pdafomi_TEMPLATE.F90 579 2020-11-23 07:32:00Z lnerger $
#include "CPP_OPTIONS.h"
CBOP
C> PDAF-OMI template observation module 
C 
C  This module handles operations for one data type (called 'module-type' below).
C  
C  The subroutines in this module are for the particular handling of
C  a single observation type.
C  The routines are called by the different call-back routines of PDAF.
C  Most of the routines are generic so that in practice only 2 routines
C  need to be adapted for a particular data type. These are the routines
C  for the initialization of the observation information (init_dim_obs)
C  and for the observation operator (obs_op).
C 
C  The module and the routines are named according to the observation type.
C  This allows to distinguish the observation type and the routines in this
C  module from other observation types.
C 
C  The module uses two derived data type (obs_f and obs_l), which contain
C  all information about the full and local observations. Only variables
C  of the type obs_f need to be initialized in this module. The variables
C  in the type obs_l are initilized by the generic routines from PDAFomi.
C 
C 
C  **Using this template:**
C    To be able to distinguish the observation type and the routines in this module,
C    we recommend to rename the module according to the observation module-type.
C    Further,we recommend to replace 'TYPE' in the routine names according to the
C    type of the observation so that they can be identified when calling them from 
C    the call-back routines.
C 
C 
C  These 2 routines need to be adapted for the particular observation type:
C  * init_dim_obs_TYPE \n
C            Count number of process-local and full observations; 
C            initialize vector of observations and their inverse variances;
C            initialize coordinate array and index array for indices of
C            observed elements of the state vector.
C  * obs_op_TYPE \n
C            observation operator to get full observation vector of this type. Here
C            one has to choose a proper observation operator or implement one.
C 
C  In addition, there are two optional routine, which are required if filters 
C  with localization are used:
C  * init_dim_obs_l_TYPE \n
C            Only required if domain-localized filters (e.g. LESTKF, LETKF) are used:
C            Count number of local observations of module-type according to
C            their coordinates (distance from local analysis domain). Initialize
C            module-internal distances and index arrays.
C  * localize_covar_TYPE \n
C            Only required if the localized EnKF is used:
C            Apply covariance localization in the LEnKF.
C 
C  __Revision history:__
C  * 2019-06 - Lars Nerger - Initial code
C  * Later revisions - see repository log
C 
      MODULE obs_eta_pdafomi

      USE mod_parallel_pdaf,
     &    ONLY: mype_filter      ! Rank of filter process
      USE PDAFomi,
     &    ONLY: obs_f, obs_l     ! Declaration of observation data types
      IMPLICIT NONE
      SAVE

#include "SIZE.h"
#include "EEPARAMS.h"
#include "DYNVARS.h"
C#include "PARAMS.h"
C#include "GRID.h"
     
C  Variables which are inputs to the module (usually set in init_pdaf)
      LOGICAL :: assim_eta         ! Whether to assimilate this data type
      REAL*8  :: rms_obs_eta       ! Observation error standard deviation
                                   ! (for constant errors)
      CHARACTER(len=32) :: eta_obsfile  ! observation file name
C  One can declare further variables, e.g. for file names which can
C  be use-included in init_pdaf() and initialized there.


C ***********************************************************************
C *** The following two data types are used in PDAFomi                ***
C *** They are declared in PDAFomi and only listed here for reference ***
C ***********************************************************************

C Data type to define the full observations by internally shared variables of the module
C   TYPE obs_f
C           Mandatory variables to be set in INIT_DIM_OBS
C      INTEGER :: doassim                   C Whether to assimilate this observation type
C      INTEGER :: disttype                  C Type of distance computation to use for localization
C                                           C (0) Cartesian, (1) Cartesian periodic
C                                           C (2) simplified geographic, (3) geographic haversine function
C      INTEGER :: ncoord                    C Number of coordinates use for distance computation
C      INTEGER, ALLOCATABLE :: id_obs_p(:,:) C Indices of observed field in state vector (process-local)
C           
C           Optional variables - they can be set in INIT_DIM_OBS
C      REAL, ALLOCATABLE :: icoeff_p(:,:)   C Interpolation coefficients for obs. operator
C      REAL, ALLOCATABLE :: domainsize(:)   C Size of domain for periodicity (<=0 for no periodicity)
C
C           Variables with predefined values - they can be changed in INIT_DIM_OBS
C      INTEGER :: obs_err_type=0            C Type of observation error: (0) Gauss, (1) Laplace
C      INTEGER :: use_global_obs=1          C Whether to use (1) global full obs. 
C                                           C or (0) obs. restricted to those relevant for a process domain
C
C           The following variables are set in the routine PDAFomi_gather_obs
C      INTEGER :: dim_obs_p                 C number of PE-local observations
C      INTEGER :: dim_obs_f                 C number of full observations
C      INTEGER :: dim_obs_g                 C global number of observations
C      INTEGER :: off_obs_f                 C Offset of this observation in overall full obs. vector
C      INTEGER :: off_obs_g                 C Offset of this observation in overall global obs. vector
C      INTEGER :: obsid                     C Index of observation over all assimilated observations
C      REAL, ALLOCATABLE :: obs_f(:)        C Full observed field
C      REAL, ALLOCATABLE :: ocoord_f(:,:)   C Coordinates of full observation vector
C      REAL, ALLOCATABLE :: ivar_obs_f(:)   C Inverse variance of full observations
C      INTEGER, ALLOCATABLE :: id_obs_f_lim(:) C Indices of domain-relevant full obs. in global vector of obs.
C                                           C (only if full obs. are restricted to process domain))
C   END TYPE obs_f

C Data type to define the local observations by internally shared variables of the module
C   TYPE obs_l
C      INTEGER :: dim_obs_l                 C number of local observations
C      INTEGER :: off_obs_l                 C Offset of this observation in overall local obs. vector
C      INTEGER, ALLOCATABLE :: id_obs_l(:)  C Indices of local observations in full obs. vector 
C      REAL, ALLOCATABLE :: distance_l(:)   C Distances of local observations
C      REAL, ALLOCATABLE :: ivar_obs_l(:)   C Inverse variance of local observations
C      INTEGER :: locweight                 C Specify localization function
C      REAL :: lradius                      C localization radius
C      REAL :: sradius                      C support radius for localization function
C   END TYPE obs_l
C ***********************************************************************

C Declare instances of observation data types used here
C We use generic names here, but one could renamed the variables
      TYPE(obs_f), TARGET, PUBLIC :: thisobs      ! full observation
      TYPE(obs_l), TARGET, PUBLIC :: thisobs_l    ! local observation

C$OMP THREADPRIVATE(thisobs_l)


C-------------------------------------------------------------------------------

      CONTAINS

C> Initialize information on the module-type observation
CC
CC The routine is called by each filter process.
CC at the beginning of the analysis step before 
CC the loop through all local analysis domains.
CC 
CC It has to count the number of observations of the
CC observation type handled in this module according
CC to the current time step for all observations 
CC required for the analyses in the loop over all local 
CC analysis domains on the PE-local state domain.
CC
CC The following four variables have to be initialized in this routine
CC * thisobs\%doassim     - Whether to assimilate this type of observations
CC * thisobs\%disttype    - type of distance computation for localization with this observaton
CC * thisobs\%ncoord      - number of coordinates used for distance computation
CC * thisobs\%id_obs_p    - index of module-type observation in PE-local state vector
CC
CC Optional is the use of
CC * thisobs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
CC * thisobs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
CC * thisobs\%obs_err_type - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
CC * thisobs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
CC                          (default: 1=use global full observations)
CC
CC Further variables are set when the routine PDAFomi_gather_obs is called.
CC
CC **Adapting the template**
CC In this routine the variables listed above have to be initialized. One
CC can include modules from the model with 'use', e.g. for mesh information.
CC Alternatively one could include these as subroutine arguments
CC
      SUBROUTINE init_dim_obs_eta(step, dim_obs)

      USE PDAFomi,
     &    ONLY: PDAFomi_gather_obs
      USE mod_assimilation,
     &    ONLY: filtertype, local_range,
     &    myIter, myThid
           

      IMPLICIT NONE

C *** Arguments ***
      INTEGER, INTENT(in)    :: step       !< Current time step
      INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

C *** Local variables ***
      INTEGER :: i, j, k, bi, bj                   ! Counters
      INTEGER :: dim_obs_p                 ! Number of process-local observations
      INTEGER :: cnt_p
      REAL*8, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
      REAL*8, ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
      REAL*8, ALLOCATABLE :: ocoord_p(:,:)   ! PE-local observation coordinates 
      REAL, ALLOCATABLE :: obs_field(:,:) 
      CHARACTER*(11) suffix

C Template reminder - delete when implementing functionality
      WRITE (*,*) 'TEMPLATE init_TYPE_pdafomi_TEMPLATE.F90: 
     &             Initialize observations'

C *********************************************
C *** Initialize full observation dimension ***
C *********************************************

      IF (mype_filter==0)
     & WRITE (*,'(8x,a)') 'Assimilate observations - OBS_eta'

C Store whether to assimilate this observation type (used in routines below)
      IF (assim_eta==.TRUE.)
     & thisobs%doassim = 1

C Specify type of distance computation
      thisobs%disttype = 0   ! 0=Cartesian

C Number of coordinates used for distance computation
C The distance compution starts from the first row
      thisobs%ncoord = 2
    
      ALLOCATE(thisobs%domainsize(2))
      thisobs%domainsize(1) = REAL(sNx)
      thisobs%domainsize(2) = REAL(sNy)

C      thisobs%use_global_obs = 1

C **********************************
C *** Read PE-local observations ***
C **********************************

C read observation values and their coordinates
C also read observation error information if available
      ALLOCATE(obs_field(sNy,sNx))
  
      IF ( eta_obsfile .NE. ' ' ) THEN
C       WRITE(suffix,'(A,I10.10)') '.', myIter
C       DO bj = myByLo(myThid), myByHi(myThid)
C        DO bi = myBxLo(myThid), myBxHi(myThid)
C         DO j=1-OLy,sNy+OLy
C          DO i=1-OLx,sNx+OLx
C           etaN(i,j,bi,bj) = 0. _d 0
C          ENDDO
C         ENDDO
C        ENDDO
C       ENDDO
C  Read an initial state
C       CALL READ_FLD_XY_RL( eta_obsfile, suffix, etaN, 0, myThid )
C  fill the overlap (+ BARRIER)
C       _EXCH_XY_RL(etaN, myThid)
       obs_field(:,:) = 0.5
      ELSE  
       obs_field(:,:) = 0.5
      ENDIF

C ***********************************************************
C *** Count available observations for the process domain ***
C *** and initialize index and coordinate arrays.         ***
C ***********************************************************

C *** Count valid observations that lie within the process sub-domain ***

      dim_obs_p = sNy*sNx


C *** Initialize vector of observations on the process sub-domain ***


C *** Initialize coordinate array of observations on the process sub-domain ***

    
      ALLOCATE(obs_p(dim_obs_p)) 
      ALLOCATE(ocoord_p(2, dim_obs_p))
      ALLOCATE(ivar_obs_p(dim_obs_p))
      ALLOCATE(thisobs%id_obs_p(1,dim_obs_p))
  
      cnt_p = 0
      Do j = 1,sNx
       Do i = 1,sNy
        cnt_p = cnt_p+1
        thisobs%id_obs_p(1,cnt_p) = cnt_p
        obs_p(cnt_p) = obs_field(i,j)
        ocoord_p(1,cnt_p) = REAL(j)
        ocoord_p(2,cnt_p) = REAL(i)
       END DO
      END DO
C *** Initialize process local index array                         ***
C *** This array holds the information which elements of the state ***
C *** vector are used in the observation operator.                 ***
C *** It has a many rows as required for the observation operator, ***
C *** i.e. 1 if observations are at grid points; >1 if             ***
C *** interpolation is required                                    ***

C The initialization is done locally for each process sub-domain and later
C used in the observation operator. 
C Examples:
C 1. If the observations are model fields located at grid points, one should
C   initialize the index array thisobs%id_obs_p with one row so that it contains 
C   the indices of the observed field values in the process-local state vector
C   (state_p). Then one can use the observation operator OBS_OP_GRIDPOINT 
C   provided by the module PDAFomi.
C 2. If the observations are the average of model fields located at grid points,
C   one should initialize the index array thisobs%id_obs_p with as many rows as 
C   values to be averaged. Each column of the arrays then contains the indices of
C   the elements of the process-local state vector that have to be averaged. With
C   this index array one can then use the observation operator OBS_OP_GRIDAVG
C   provided by the module PDAFomi.
C 3. If model values need to be interpolated to the observation location
C   one should initialize the index array thisobs%id_obs_p with as many rows as 
C   values are required in the interpolationto be averaged. Each column of the 
C   array then contains the indices of elements of the process-local state vector 
C   that are used in the interpolation.
C Below, you need to replace NROWS by the number of required rows

C    ALLOCATE(thisobs%id_obs_p(1, dim_obs_p))

C    thisobs%id_obs_p = ...


C **********************************************************************
C *** Initialize interpolation coefficients for observation operator ***
C **********************************************************************

C This initialization is only required if an observation operator
C with interpolation is used. The coefficients should be determined
C here instead of the observation operator, because the operator is
C called for each ensemble member while init_dim_obs is only called
C once.
C Allocate array of interpolation coefficients. As thisobs%id_obs_p, the number
C of rows corresponds to the number of grid points using the the interpolation

C    ALLOCATE(thisobs%icoeff_p( NROWS , dim_obs_p))

C Ensure that the order of the coefficients is consistent with the
C indexing in thisobs%id_obs_p. Further ensure that the order is consistent
C with the assumptions used in the observation operator.

C    thisobs%icoeff_p = ...


C ****************************************************************
C *** Define observation errors for process-local observations ***
C ****************************************************************

C    ALLOCATE(ivar_obs_p(dim_obs_p))
    
      ivar_obs_p(:) = 1.0 / (rms_obs_eta*rms_obs_eta)


C ****************************************
C *** Gather global observation arrays ***
C ****************************************

C NOTE FOR DIM_OBS_P=0
C For the call to PDAFomi_gather_obs_f, obs_p, ivar_obs_p, ocoord_p,
C and thisobs%id_obs_p need to be allocated. Thus, if dim_obs_p=0 can 
C happen in your application you should explicitly handle this case.
C You can introduce an IF block in the initializations above:
C  IF dim_obs_p>0 THEN
C     regular allocation and initialization of obs_p, ivar_obs_p, ocoord_p
C  ELSE
C     allocate obs_p, ivar_obs_p, ocoord_p, thisobs%id_obs_p with size=1
C  ENDIF


C This routine is generic for the case that only the observations, 
C inverse variances and observation coordinates are gathered

      CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p,
     &     ocoord_p, thisobs%ncoord, local_range, dim_obs)


C *********************************************************
C *** For twin experiment: Read synthetic observations  ***
C *********************************************************

C   IF (twin_experiment .AND. filtertype/=11) THEN
C      CALL read_syn_obs(file_syntobs_TYPE, dim_obs, thisobs%obs_f, 0, 1-mype_filter)
C   END IF


C ********************
C *** Finishing up ***
C ********************

C Deallocate all local arrays
      DEALLOCATE(obs_field)
      DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

C Arrays in THISOBS have to be deallocated after the analysis step
C by a call to deallocate_obs() in prepoststep_pdaf.

      END SUBROUTINE init_dim_obs_eta



C-------------------------------------------------------------------------------
C> Implementation of observation operator 
CC
CC This routine applies the full observation operator
CC for the type of observations handled in this module.
CC
CC One can choose a proper observation operator from
CC PDAFOMI_OBS_OP or add one to that module or 
CC implement another observation operator here.
CC
CC The routine is called by all filter processes.
CC
      SUBROUTINE obs_op_eta(dim_p, dim_obs, state_p, ostate)

      USE PDAFomi,
     & ONLY: PDAFomi_obs_op_gridpoint

      IMPLICIT NONE

C *** Arguments ***
      INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
      INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
      REAL*8, INTENT(in)    :: state_p(dim_p)        !< PE-local model state
      REAL*8, INTENT(inout) :: ostate(dim_obs)       !< Full observed state


C Template reminder - delete when implementing functionality
      WRITE (*,*) 'TEMPLATE init_TYPE_pdafomi_TEMPLATE.F90: 
     &             Apply observation operator'

C ******************************************************
C *** Apply observation operator H on a state vector ***
C ******************************************************

      IF (thisobs%doassim==1) THEN 
C+++  Choose suitable observation operator from the
C+++  module PDAFomi_obs_op or implement your own

C=== hajoon
C   print *,'ostate_0:', ostate
C       print *, 'before PDAFomi_obs_op_gridpoint'
C       print *, 'state_p(299) = ', state_p(299)
C=== hajoon
C Example: Observation operator for observed grid point values
       CALL PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)
       WRITE (*,*) 'obs_op_eta is running'
C=== hajoon
C   print *,'ostate_1:', ostate
C       print *, 'after PDAFomi_obs_op_gridpoint'
C       print *, 'state_p(299) = ', state_p(299)
C=== hajoon
      END IF

      END SUBROUTINE obs_op_eta



C-------------------------------------------------------------------------------
C> Initialize local information on the module-type observation
CC
CC The routine is called during the loop over all local
CC analysis domains. It has to initialize the information
CC about local observations of the module type. It returns
CC number of local observations of the module type for the
CC current local analysis domain in DIM_OBS_L and the full
CC and local offsets of the observation in the overall
CC observation vector.
CC
CC This routine calls the routine PDAFomi_init_dim_obs_l
CC for each observation type. The call allows to specify a
CC different localization radius and localization functions
CC for each observation type and  local analysis domain.
CC
      SUBROUTINE init_dim_obs_l_eta(domain_p, step, dim_obs, dim_obs_l)

C Include PDAFomi function
      USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l

C Include localization radius and local coordinates
C one can also set observation-specific values for the localization.
      USE mod_assimilation,
     & ONLY: coords_l, local_range, locweight, srange

      IMPLICIT NONE

C *** Arguments ***
      INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
      INTEGER, INTENT(in)  :: step         !< Current time step
      INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
      INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector


C Template reminder - delete when implementing functionality
      WRITE (*,*) 'TEMPLATE init_TYPE_pdafomi_TEMPLATE.F90: 
     &             Initialize local observations'

C **********************************************
C *** Initialize local observation dimension ***
C **********************************************

C Here one has to specify the coordinates of the local analysis domain
C (coords_l) and the localization variables, which can be different for
C each observation type and can be made dependent on the index DOMAIN_P.
C coords_l should be set in the call-back routine init_dim_l.

      CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l,
     &   locweight, local_range, srange, dim_obs_l)

      END SUBROUTINE init_dim_obs_l_eta



C-------------------------------------------------------------------------------
C> Perform covariance localization for local EnKF on the module-type observation
CC
CC The routine is called in the analysis step of the localized
CC EnKF. It has to apply localization to the two matrices
CC HP and HPH of the analysis step for the module-type
CC observation.
CC
CC This routine calls the routine PDAFomi_localize_covar
CC for each observation type. The call allows to specify a
CC different localization radius and localization functions
CC for each observation type.
CC
      SUBROUTINE localize_covar_eta(dim_p, dim_obs, HP_p, HPH, coords_p)

C Include PDAFomi function
      USE PDAFomi, ONLY: PDAFomi_localize_covar

C Include localization radius and local coordinates
      USE mod_assimilation,
     & ONLY: local_range, locweight, srange

      IMPLICIT NONE

C *** Arguments ***
      INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
      INTEGER, INTENT(in) :: dim_obs               !< Dimension of observation vector
      REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
      REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
      REAL, INTENT(in)    :: coords_p(:,:)         !< Coordinates of state vector elements


C Template reminder - delete when implementing functionality
      WRITE (*,*) 'TEMPLATE init_TYPE_pdafomi_TEMPLATE.F90: 
     &             Apply covariance localization'

C *************************************
C *** Apply covariance localization ***
C *************************************

C Here one has to specify the three localization variables
C which can be different for each observation type.

      CALL PDAFomi_localize_covar(thisobs, dim_p, locweight,
     &     local_range, srange, coords_p, HP_p, HPH)

      END SUBROUTINE localize_covar_eta

      END MODULE obs_eta_pdafomi
