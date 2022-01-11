C     $Id: init_ens_pdaf.F 864 2021-11-19 12:30:17Z lnerger $
C
C     !ROUTINE: init_ens_pdaf --- Initialize ensemble
C
C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C
C     The routine is called when the filter is
C     initialized in PDAF\_filter\_init.  It has
C     to initialize an ensemble of dim\_ens states.
C     Typically, the ensemble will be directly read from files.
C
C     The routine is called by all filter processes and
C     initializes the ensemble for the PE-local domain.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     2022-01 - Hajoon Song - Add an ability to read initial ensemble
C                             from files
C     Later revisions - see svn log
C
C     !INTERFACE:
      SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p,
     &     Uinv, ens_p, flag)

C     !USES:
      USE mod_parallel_pdaf,
     &     ONLY: mype_filter
      USE mod_assimilation,
     &     ONLY: myIter, myThid, read_ens

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: filtertype                  ! Type of filter to initialize
      INTEGER, INTENT(in) :: dim_p                       ! PE-local state dimension
      INTEGER, INTENT(in) :: dim_ens                     ! Size of ensemble
      REAL*8, INTENT(inout) :: state_p(dim_p)            ! PE-local model state
C     It is not necessary to initialize the array 'state_p' for SEIK. 
C     It is available here only for convenience and can be used freely.
      REAL*8, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
      REAL*8, INTENT(out)   :: ens_p(dim_p, dim_ens)     ! PE-local state ensemble
      INTEGER, INTENT(inout) :: flag                     ! PDAF status flag

C     !CALLING SEQUENCE:
C     Called by: PDAF_filter_init    (as U_ens_init)


#include "CPP_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DYNVARS.h"
#include "PARAMS.h"
#include "GRID.h"

C     *** local variables ***
      INTEGER :: i, j, k, bi, bj, member  ! Counters
      INTEGER :: localWarnings
      REAL*8 :: Tfreezing
      CHARACTER(len=4) :: ensstr          ! String for ensemble member
      CHARACTER*(MAX_LEN_MBUF) msgBuf

C     i,j,k  :: Loop counters
C     bi,bj  :: tile indices

C     **********************
C     *** INITIALIZATION ***
C     **********************

C     *** Generate full ensemble on filter-PE 0 ***
      IF (mype_filter==0) THEN
         WRITE (*, '(/a, 9x, a)') 'PDAFuser: ',
     &        'Initialize state ensemble'
         WRITE (*, '(a, 9x, a)') 'PDAFuser: ',
     &        '--- read ensemble from files'
         WRITE (*, '(a, 9x, a, i5)') 'PDAFuser: ',
     &        '--- Ensemble size:  ', dim_ens
      END IF


C     ********************************
C     *** Read ensemble from files ***
C     ********************************
      IF (read_ens) THEN
       IF (mype_filter==0) THEN
         WRITE (*, '(a, 9x, a)') 'PDAFuser: ',
     &        '--- read ensemble from files'
       END IF
C     Set model variables to initial/restart conditions
       DO member = 1, dim_ens
        IF ( startTime .EQ. baseTime .AND.  nIter0 .EQ. 0
     &     .AND. pickupSuff .EQ. ' ' ) THEN
         WRITE(ensstr,'(A,I3.3)') '_',member
         CALL INI_PSURF_ENS( myThid, ensstr )
         CALL INI_THETA_ENS( myThid, ensstr )
         CALL INI_SALT_ENS( myThid, ensstr )
         CALL INI_VEL_ENS( myThid, ensstr )

C      Or Read pickup files
        ELSEIF ( .NOT.useOffLine .OR. nonlinFreeSurf.GT.0 ) THEN
         CALL READ_PICKUP( nIter0, myThid )
        ENDIF

        CALL collect_state_pdaf(dim_p, state_p)
        ens_p(:, member)  = state_p(:)
       ENDDO
C     ********************************
C     *** Simply create ensemble   ***
C     ********************************
      ELSE
       IF (mype_filter==0) THEN
         WRITE (*, '(a, 9x, a)') 'PDAFuser: ',
     &        '--- create ensemble naively'
       END IF
       CALL collect_state_pdaf(dim_p, state_p)

       DO member = 1, dim_ens
          ens_p(:, member)  = state_p(:) + float(member-1)*1.d-2
       ENDDO

      ENDIF

C     ****************
C     *** clean up ***
C     ****************

      END SUBROUTINE init_ens_pdaf
