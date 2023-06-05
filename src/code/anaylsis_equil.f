      SUBROUTINE SUM(HB,HH)

#include "def.h"

      USE analysis_GLOBAL
	use inputreadin

      IMPLICIT NONE

      INTEGER HB,HH
        
      IF (COLL .LE. 2*QUARTER) THEN
	 HB_SUM(2) = HB_SUM(2) + HB
	 HH_SUM(2) = HH_SUM(2) + HH
	 QT_COUNT(2) = QT_COUNT(2) + 1
      ELSEIF (COLL .LE. 3*QUARTER) THEN
         HB_SUM(3) = HB_SUM(3) + HB
         HH_SUM(3) = HH_SUM(3) + HH
         QT_COUNT(3) = QT_COUNT(3) + 1
      ELSE
         HB_SUM(4) = HB_SUM(4) + HB
         HH_SUM(4) = HH_SUM(4) + HH
         QT_COUNT(4) = QT_COUNT(4) + 1
      ENDIF

      RETURN
      
      END

!*********************************************************************************

      SUBROUTINE AVERAGE(PER_HB,PER_HH)

#include "def.h"

      USE GLOBAL

      IMPLICIT NONE

      INTEGER K
      REAL PER_HB,PER_HH
      REAL AV_HB(4),AV_HH(4),SQSUM_HB,SUMSQ_HB,SQSUM_HH,SUMSQ_HH
      character*64 filename
      
      SQSUM_HB = 0.0
      SUMSQ_HB = 0.0
      SQSUM_HH = 0.0
      SUMSQ_HH = 0.0

      DO K = 2, 4
	 IF (QT_COUNT(K) .NE. 0) THEN
            AV_HB(K) = REAL(HB_SUM(K))/QT_COUNT(K)
            AV_HH(K) = REAL(HH_SUM(K))/QT_COUNT(K)
	 ELSE
	    AV_HB(K) = 0.0
	    AV_HH(K) = 0.0
	 END IF
         SQSUM_HB = SQSUM_HB + AV_HB(K)**2
         SQSUM_HH = SQSUM_HH + AV_HH(K)**2
	 SUMSQ_HB = SUMSQ_HB + AV_HB(K)
	 SUMSQ_HH = SUMSQ_HH + AV_HH(K)
      END DO      
      SUMSQ_HB = SUMSQ_HB**2
      SUMSQ_HH = SUMSQ_HH**2
      
      IF (SUMSQ_HB .NE. 0.0) THEN
         PER_HB = SQRT(ABS(3.0*SQSUM_HB - SUMSQ_HB)/(2.0/3.0*SUMSQ_HB))*100.0
      ELSE
	 PER_HB = 0.0
      ENDIF
      IF (SUMSQ_HH .NE. 0.0) THEN
         PER_HH = SQRT(ABS(3.0*SQSUM_HH - SUMSQ_HH)/(2.0/3.0*SUMSQ_HH))*100.0
      ELSE
	 PER_HH = 0.0
      ENDIF
	
      WRITE(6,*) 'AV_HB', AV_HB(2), AV_HB(3), AV_HB(4), PER_HB
      WRITE(6,*) 'AV_HH', AV_HH(2), AV_HH(3), AV_HH(4), PER_HH

      IF (AV_HB(4) .LE. 10.0) PER_HB = 0.0
      IF (AV_HH(4) .LE. 10.0) PER_HH = 0.0

      RETURN

      END

!*********************************************************************************

      SUBROUTINE DOUBLING

#include "def.h"

      USE GLOBAL

      IMPLICIT NONE

      HB_SUM(2) = HB_SUM(3) + HB_SUM(4)
      HH_SUM(2) = HH_SUM(3) + HH_SUM(4)
      QT_COUNT(2) = QT_COUNT(3) + QT_COUNT(4)

      HB_SUM(3) = 0
      HB_SUM(4) = 0

      HH_SUM(3) = 0
      HH_SUM(4) = 0

      QT_COUNT(3) = 0
      QT_COUNT(4) = 0

      RETURN

      END

!*********************************************************************************
