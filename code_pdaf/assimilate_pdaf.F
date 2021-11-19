C     $Id: assimilate_pdaf.F 7 2018-01-31 18:11:38Z lnerger $
#include "CPP_OPTIONS.h"
C
C     !ROUTINE: assimilate_pdaf - Routine to control perform analysis step
C
C     !INTERFACE:
      SUBROUTINE assimilate_pdaf( myTimeLoc, myIterLoc, myThidLoc )

C     !DESCRIPTION:
C     This routine is called during the model integrations at each time 
C     step. It check whether the forecast phase is completed. If so, 
C     PDAF_put_state_X is called to perform the analysis step.
C
C     !REVISION HISTORY:
C     2013-08 - Lars Nerger - Initial code for NEMO
C     Later revisions - see svn log
C
C     !USES:
      USE mod_parallel_model,
     &     ONLY: mype_world, abort_parallel
      USE mod_assimilation,
     &     ONLY: filtertype, myIter, myThid, myTime
C-hajoon
      ! Interface definitions to PDAF core routines
      USE pdaf_interfaces_module,
     &     ONLY: PDAFomi_assimilate_local, PDAFomi_assimilate_global,
     &     PDAFomi_assimilate_lenkf, PDAFomi_generate_obs,
     &     PDAF_get_localfilter
C-hajoon
      
      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: myIterLoc
      INTEGER, INTENT(in) :: myThidLoc
      _RL , INTENT(in)    :: myTimeLoc

C     !CALLING SEQUENCE:
C     Called by: step
C     CAlls: PDAF_assimilate_X

C     Local variables
      INTEGER status_pdaf    ! PDAF status flag
      INTEGER localfilter
C     External subroutines
      EXTERNAL :: collect_state_pdaf, ! Routine to collect a state vector from model fields
C     &     init_dim_obs_pdaf,   ! Initialize Dimension Of Observation Vector
C     &     obs_op_pdaf,         ! Implementation of the Observation operator
     &     init_obs_pdaf,       ! Routine to provide vector of measurements
     &     prepoststep_ens_pdaf,  ! User supplied pre/poststep routine
     &     prodRinvA_pdaf,      ! Provide product R^-1 A for some matrix A
     &     init_obsvar_pdaf,    ! Initialize mean observation error variance
     &     next_observation_pdaf, ! Provide time step of next observation
     &     distribute_state_pdaf  ! Routine to distribute a state vector to model fields
      EXTERNAL :: init_n_domains_pdaf, ! Provide number of local analysis domains
     &     init_dim_l_pdaf,     ! Initialize state dimension for local ana. domain
C     &     init_dim_obs_l_pdaf, ! Initialize dim. of obs. vector for local ana. domain
     &     g2l_state_pdaf,      ! Get state on local ana. domain from global state
     &     l2g_state_pdaf,      ! Init global state from state on local analysis domain
     &     g2l_obs_pdaf,        ! Restrict a global obs. vector to local analysis domain
     &     init_obs_l_pdaf,     ! Provide vector of measurements for local ana. domain
     &     prodRinvA_l_pdaf,    ! Provide product R^-1 A for some local matrix A
     &     init_obsvar_l_pdaf,  ! Initialize local mean observation error variance
     &     init_obs_f_pdaf,     ! Provide full vector of measurements for PE-local domain
C     &     obs_op_f_pdaf,       ! Obs. operator for full obs. vector for PE-local domain
     &     init_dim_obs_f_pdaf  ! Get dimension of full obs. vector for PE-local domain
C-hajoon
C    ! Interface to PDAF-OMI for local and global filters
      EXTERNAL :: init_dim_obs_pdafomi, ! Get dimension of full obs. vector for PE-local domain
     &     obs_op_pdafomi,          ! Obs. operator for full obs. vector for PE-local domain
     &     init_dim_obs_l_pdafomi,  ! Get dimension of obs. vector for local analysis domain
     &     localize_covar_pdafomi   ! Apply localization to covariance matrix in LEnKF
C-hajoon
C      EXTERNAL :: get_obs_f_pdaf
      
C     *********************************
C     *** Call assimilation routine ***
C     *********************************
      CALL PDAF_get_localfilter(localfilter)

      myIter = myIterLoc
      myThid = myThidLoc
      myTime = myTimeLoc

      IF (localfilter==1) THEN
       CALL PDAFomi_assimilate_local(collect_state_pdaf,
     &      distribute_state_pdaf, init_dim_obs_pdafomi, 
     &      obs_op_pdafomi, prepoststep_ens_pdaf,
     &      init_n_domains_pdaf, init_dim_l_pdaf,
     &      init_dim_obs_l_pdafomi, g2l_state_pdaf,
     &      l2g_state_pdaf, next_observation_pdaf,status_pdaf)

      ELSE
         IF (filtertype==8) THEN
            CALL PDAFomi_assimilate_lenkf(collect_state_pdaf, 
     &      distribute_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi,
     &      prepoststep_ens_pdaf, localize_covar_pdafomi, 
     &      next_observation_pdaf, status_pdaf)
C         ELSE IF (filtertype == 11) THEN
C            CALL PDAFomi_generate_obs(collect_state_pdaf, 
C     &      distribute_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi,
C     &      get_obs_f_pdaf, prepoststep_ens_pdaf, 
C     &      next_observation_pdaf, status_pdaf)
         ELSE  
            CALL PDAFomi_assimilate_global(collect_state_pdaf, 
     &      distribute_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi,
     &      prepoststep_ens_pdaf, 
     &      next_observation_pdaf, status_pdaf)
         END IF
      END IF

      
C     Check for errors during execution of PDAF
      IF (status_pdaf /= 0) THEN
         WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') 
     &        'ERROR ', status_pdaf, 
     &        ' in PDAF_put_state - stopping! (PE ', mype_world,')'
         CALL  abort_parallel()
      END IF

      END 
