c     UC Copyright Notice
c
c     This software is Copyright (c) 2014-2015 The Regents of the 
c     University of California. All Rights Reserved.
c
c     Permission to copy and modify this software and its documentation
c     for educational, research and non-profit purposes, without fee, 
c     and without a written agreement is hereby granted, provided that
c     the above copyright notice, this paragraph and the following three
c     paragraphs appear in all copies.
c
c     Permission to make commercial use of this software may be obtained
c     by contacting:
c
c     Technology Transfer Office
c     9500 Gilman Drive, Mail Code 0910
c     University of California
c     La Jolla, CA 92093-0910
c     (858) 534-5815
c     invent@ucsd.edu
c
c     This software program and documentation are copyrighted by The
c     Regents of the University of California. The software program and
c     documentation are supplied "as is", without any accompanying
c     services from The Regents. The Regents does not warrant that the
c     operation of the program will be uninterrupted or error-free. The
c     end-user understands that the program was developed for research
c     purposes and is advised not to rely exclusively on the program for
c     any reason.
c
c     IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY 
c     PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL 
c     DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS 
c     SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF 
c     CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
c     THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY 
c     WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
c     OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE 
c     SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE 
c     UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE 
c     MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

#include "cvFlowsolverOptions.h"

#if (VER_CLOSEDLOOP == 1)

!> This subroutine initiate the General boundary condition
 
      SUBROUTINE initSvZeroD

      USE GeneralBC
      INCLUDE "global.h"
      INCLUDE "mpif.h"
      INCLUDE "common_blocks/workfc.h"
      INCLUDE "common_blocks/conpar.h"
      INCLUDE "common_blocks/nomodule.h"

      REAL*8 temp(nshg)
      INTEGER loopA
      LOGICAL ierr

      ALLOCATE(nsrflistCoupled(0:MAXSURF), ECoupled(0:MAXSURF),
     2   QCoupled(0:MAXSURF), QnCoupled(0:MAXSURF), ACoupled(0:MAXSURF),
     2   PCoupled(0:MAXSURF), PnCoupled(0:MAXSURF), PGenDer(0:MAXSURF))

      CHARACTER(len=50), ALLOCATABLE :: svzd_blk_names(:)

      ALLOCATE(svzd_blk_names(0:numCoupledSrfs))

      IF (ipvsq .GE. 2) THEN
         PGenDerFlag = .TRUE.
      ELSE
         PGenDerFlag = .FALSE.
      END IF
      PGenDer = 0D0     

      nsrflistCoupled = nsrflistDirichlet
      nsrflistCoupled(numDirichletSrfs+1:numCoupledSrfs) = 
     2   nsrflistNeumann(1:numNeumannSrfs)
      
c     Eval Area of Coupled Surfaces
      temp = 1D0
      CALL integrScalar(ACoupled, temp, nsrflistCoupled, numCoupledSrfs)
      
      IF (myrank .EQ. master) THEN

c       CHECK AREA
        DO loopA = 1,numCoupledSrfs
          IF(ACoupled(loopA).LE.0.0D0)THEN
            PRINT *,'' 
            PRINT *,'WARNING: Zero Area for Coupled Surface ',loopA
            PRINT *,'' 
          ENDIF
        ENDDO

        !TODO
        ! Define library path in modules probably
        ! Read svzerod input file from solver.inp
        !CALL lpn_interface_add_model_(const char* lpn_library_name, const
        !char* lpn_json_file, int* block_id)
        ! Allocate solution vector
        ! Need to figure out how to map the order of *Coupled(:) to
        ! svZeroD blocks
        ! Save block names in correct order (for parameter updates) and
        ! position of relevant variables in state vector (for return)
        !TODO

!        IF (iGenFromFile .EQ. 1) THEN
!           INQUIRE(FILE='GenBC',EXIST=ierr)
!           IF (.NOT.ierr) THEN
!              INQUIRE(FILE='../GenBC',EXIST=ierr)
!              IF (ierr) THEN
!                 CALL system('cp ../GenBC ./')
!              ELSE
!                 PRINT *, 'ERROR: No executable GenBC has been found'
!                 PRINT *, 'Please, make sure that you had this', 
!    2               ' executable file inside runnig directory'
!                 STOP
!              END IF
!           END IF
!        END IF

!        INQUIRE(FILE='InitialData',EXIST=ierr)
!        IF (ierr) THEN
!           PRINT *, ''
!           PRINT *, 'NOTE: Initializing General BC form previous simulation'
!           PRINT *, ''
!        ELSE
!           INQUIRE(FILE='../InitialData',EXIST=ierr)
!           IF (ierr) THEN
!              CALL system('cp ../InitialData ./')
!              PRINT *, ''
!              PRINT *,'NOTE: General BC has been initialized from provided',
!    2            ' file'
!              PRINT *, ''
!           ELSE
!              PRINT *, ''
!              PRINT *, 'NOTE: Self GenBC initializiation'
!              PRINT *, ''
!           END IF
!        END IF
      END IF

      RETURN
      END SUBROUTINE initSvZeroD


!> Main subroutine for gathering required information and 
!! communicating with 0D to acquire necessary information

      SUBROUTINE calcSvZeroD (y, yold)
       
      USE GeneralBC
      INCLUDE "global.h"
      INCLUDE "mpif.h"
      INCLUDE "common_blocks/conpar.h"
      INCLUDE "common_blocks/workfc.h"
      INCLUDE "common_blocks/inpdat.h"
      INCLUDE "common_blocks/matdat.h"
      INCLUDE "common_blocks/nomodule.h"

      INTEGER MPIstat(MPI_STATUS_SIZE), i, ia

      REAL*8 y(nshg,ndof), yold(nshg,ndof), QInitial, en(nshg,ndof), 
     2       enOfEle,rho

      REAL*8 testnorm
      INTEGER ierr

c     Get Density
      rho = datmat(1,1,1)

      DO ia=1, nshg
         enOfEle = yold(ia,4)
         DO i=1, nsd
           enOfEle = enOfEle + 5D-1*rho*yold(ia,i)*yold(ia,i)
         END DO
         DO i=1, nsd
            en(ia,i) = enOfEle*yold(ia,i)
         END DO
      END DO
      
      CALL GetFlowQ (QCoupled,  y,    nsrflistCoupled, numCoupledSrfs) 
      CALL GetFlowQ (QnCoupled, yold, nsrflistCoupled, numCoupledSrfs)
      CALL GetFlowQ (ECoupled,  en,   nsrflistCoupled, numCoupledSrfs)

      CALL integrScalar (PCoupled, y(:,4),     nsrflistCoupled, 
     2   numCoupledSrfs)

      CALL integrScalar (PnCoupled, yold(:,4), nsrflistCoupled, 
     2   numCoupledSrfs)
      
      DO i=1, numCoupledSrfs
         PnCoupled(i) = PnCoupled(i)/ACoupled(i)
         PCoupled(i) = PCoupled(i)/ACoupled(i)
      END DO

      IF (myrank .EQ. master) THEN
         IF (GenFlag .EQ. 'L') THEN
            i = numCoupledSrfs
            CALL printSvZeroD (i, nsrflistCoupled(1:i), QCoupled(1:i), 
     2         PCoupled(1:i), ECoupled(1:i))
         END IF

         ! TODO
         ! Update block paramaters by sending PnCoupled, PCoupled,
         ! QnCoupled, QCoupled to svZeroD
         ! Run svZeroD
         ! Extract PnCoupled, PCoupled, QnCoupled, QCoupled from state
         ! vector or from new function that returns this (in
         ! lpn_interface)
         ! TODO

!        OPEN (1,FILE='GenBC.int', STATUS='UNKNOWN', FORM='UNFORMATTED')
!        WRITE (1) GenFlag
!        WRITE (1) Delt(1)
!        WRITE (1) numDirichletSrfs
!        WRITE (1) numNeumannSrfs
!        DO i=1, numCoupledSrfs
!           IF (i .LE. numDirichletSrfs) THEN
!              WRITE (1) PnCoupled(i), PCoupled(i)
!           ELSE
!              WRITE (1) QnCoupled(i), QCoupled(i)
!           END IF
!        END DO
!        CLOSE (1)

!        
!        EXTERNAL CALL TO GENBC
!        IF (iGenFromFile .EQ. 1) THEN            
!           CALL system('./GenBC')
!        ELSE
!           CALL system('GenBC')
!        END IF

!        OPEN (1,FILE='GenBC.int',STATUS='OLD', FORM='UNFORMATTED')
!        DO i=1, numCoupledSrfs
!           IF (i .LE. numDirichletSrfs) THEN
!              READ (1) QCoupled(i)
!           ELSE
!              READ (1) PCoupled(i)
!           END IF
!        END DO
!        CLOSE(1)
!     END IF
      
      i = MAXSURF + 1
      CALL MPI_BCAST(QCoupled, i, MPI_DOUBLE_PRECISION, master,
     2   MPI_COMM_WORLD, ierr)
      
      CALL MPI_BCAST(PCoupled, i, MPI_DOUBLE_PRECISION, master,
     2   MPI_COMM_WORLD, ierr)
      
      RETURN
      END SUBROUTINE calcsvZeroD


!> Writting results to the disk

      SUBROUTINE printSvZeroD (nSrfs, surfID, Q, P, E)
       
      IMPLICIT NONE

      LOGICAL ierr

      INTEGER nSrfs, nParam, i
      PARAMETER (nParam = 3)
      INTEGER surfID(nSrfs)

      REAL*8 Q(nSrfs), P(nSrfs), E(nSrfs), R(nParam,nSrfs)
 
      CHARACTER*32 fileNames(nParam)
      CHARACTER(len=32) :: myFMT1,myFMT2

      IF (nSrfs .EQ. 0) RETURN

      fileNames = (/'QGeneral','PGeneral','EGeneral'/)

      R(1,:) = Q
      R(2,:) = P
      R(3,:) = E

c     SET FORMATS
      WRITE(myFMT1, '("(",I0,"(E13.5))")') nSrfs
      WRITE(myFMT2, '("(",I0,"(I13))")') nSrfs

      DO i=1, nParam
         INQUIRE(FILE=TRIM(fileNames(i)), EXIST=ierr)
         IF (ierr) THEN
            OPEN (1, FILE=TRIM(fileNames(i)), STATUS='OLD', 
     2         ACCESS='APPEND')

            WRITE (1,fmt=myFMT1) R(i,:)
            CLOSE (1)
         ELSE
            OPEN (1, FILE=TRIM(fileNames(i)), STATUS='NEW')
            WRITE (1,fmt=myFMT2) surfID
            WRITE (1,fmt=myFMT1) R(i,:)
            CLOSE (1)
         END IF
      END DO

      RETURN
      END SUBROUTINE printSvZeroD

#endif
