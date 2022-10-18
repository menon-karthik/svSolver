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

      USE svZeroD
      INCLUDE "global.h"
      INCLUDE "mpif.h"
      INCLUDE "common_blocks/workfc.h"
      INCLUDE "common_blocks/conpar.h"
      INCLUDE "common_blocks/nomodule.h"
      INCLUDE "common_blocks/inpdat.h"

      REAL*8 temp(nshg)
      LOGICAL ierr
      
      CHARACTER(len=160) svzerod_library, svzerod_file
      INTEGER :: ids(0:1)
      CHARACTER(len=50), ALLOCATABLE :: svzd_blk_names_unsrtd(:)
      INTEGER i,j,found

      ALLOCATE(nsrflistCoupled(0:MAXSURF), ECoupled(0:MAXSURF),
     2   QCoupled(0:MAXSURF), QnCoupled(0:MAXSURF), ACoupled(0:MAXSURF),
     2   PCoupled(0:MAXSURF), PnCoupled(0:MAXSURF), PDer(0:MAXSURF))

      IF (ipvsq .GE. 2) THEN
         PDerFlag = .TRUE.
      ELSE
         PDerFlag = .FALSE.
      END IF
      PDer = 0D0     

      nsrflistCoupled = nsrflistDirichlet
      nsrflistCoupled(numDirichletSrfs+1:numCoupledSrfs) = 
     2   nsrflistNeumann(1:numNeumannSrfs)
      
c     Eval Area of Coupled Surfaces
      temp = 1D0
      CALL integrScalar(ACoupled, temp, nsrflistCoupled, numCoupledSrfs)
      
      IF (myrank .EQ. master) THEN

c       CHECK AREA
        DO i = 1,numCoupledSrfs
          IF (ACoupled(i) .LE. 0.0D0) THEN
            PRINT *,'' 
            PRINT *,'WARNING: Zero Area for Coupled Surface ',i
            PRINT *,'' 
          END IF
        END DO

        ! Open the interface file
        INQUIRE(FILE='svZeroD_interface.dat',EXIST=ierr)
        IF (.NOT.ierr) THEN
          PRINT *, 'ERROR: svZeroD_interface.dat not found'
          STOP
        END IF
        OPEN (113,FILE='svZeroD_interface.dat',STATUS='OLD')
        
        ! Read the svZeroD library location
        READ (113, *) ! Header
        READ (113, '(A)') svzerod_library
        READ (113, *) ! Blank line
        
        ! Read svZeroD input config file name
        READ (113, *) ! Header
        READ (113, '(A)') svzerod_file
        READ (113, *) ! Blank line
        
        ! Read svZeroD blocks names and surface IDs
        ALLOCATE(svzd_blk_names(numCoupledSrfs))
        ALLOCATE(svzd_blk_names_unsrtd(numCoupledSrfs))
        ALLOCATE(svzd_blk_ids(numCoupledSrfs))
        ALLOCATE(svzd_blk_name_len(numCoupledSrfs))
        READ (113, *) ! Header
        DO i = 1, numCoupledSrfs
          READ (113, *) svzd_blk_names_unsrtd(i), svzd_blk_ids(i)
        END DO
        CLOSE(113)

        ! Arrange svzd_blk_names in the same order as surface IDs in nsrflistCoupled
        DO i = 1, numCoupledSrfs
          found = 0
          DO j = 1, numCoupledSrfs
            IF (svzd_blk_ids(j) .EQ. nsrflistCoupled(i)) THEN
              found = 1
              svzd_blk_names(i) = svzd_blk_names_unsrtd(j)
              svzd_blk_name_len(i) = LEN(TRIM(svzd_blk_names(i)))
              EXIT
            END IF
          ENDDO
          IF (found .EQ. 0) THEN
            WRITE(*,*) "ERROR: Did not find block name for
     &                surface ID: ", nsrflistCoupled(i)
            STOP
          END IF
        END DO

        CALL lpn_interface_add_model(TRIM(svzerod_library),
     &                           LEN(TRIM(svzerod_library)),
     &           TRIM(svzerod_file),LEN(TRIM(svzerod_file)),
     &               model_id,num_output_steps,system_size)

        CALL lpn_interface_set_external_step_size(model_id,Delt(1))

        ! Allocate variables for svZeroD solution
        ALLOCATE(lpn_times(0:num_output_steps))
        ALLOCATE(lpn_solutions(0:num_output_steps*system_size))
        ALLOCATE(lpn_state_y(0:system_size))
        ALLOCATE(last_state_y(0:system_size))
        ALLOCATE(last_state_ydot(0:system_size))
        svZeroDTime = 0.0D0

        ! Save IDs of relevant variables in the solution vector
        ALLOCATE(sol_IDs(2*numCoupledSrfs))
        DO i = 1, numCoupledSrfs
          CALL lpn_interface_get_variable_ids(model_id,
     &                         TRIM(svzd_blk_names(i)),
     &                        svzd_blk_name_len(i),ids)
          sol_IDs(2*(i-1)+1) = ids(0)
          sol_IDs(2*(i-1)+2) = ids(1)
        END DO

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
       
      USE svZeroD
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
      INTEGER error_code

      REAL*8 :: params(2), times(2)

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
         IF (svzerodFlag .EQ. 'L') THEN !Last iteration
            i = numCoupledSrfs
            CALL printSvZeroD (i, nsrflistCoupled(1:i), QCoupled(1:i), 
     2         PCoupled(1:i), ECoupled(1:i))
         END IF

         IF (svzerodFlag .NE. 'I') THEN
           
            IF (svZeroDTime > 0.0D0) THEN
               ! Set initial condition from previous state
               CALL lpn_interface_update_state(model_id, last_state_y,
     &                                                last_state_ydot)
            END IF

            times = (/svZeroDTime, svZeroDTime+Delt(1)/)

            ! Update pressure and flow in the zeroD model
            DO i=1, numCoupledSrfs
               IF (i .LE. numDirichletSrfs) THEN
                  params = (/PnCoupled(i), PCoupled(i)/)
               ELSE
                  params = (/QnCoupled(i), QCoupled(i)/)
               END IF
               CALL lpn_interface_update_block_params(model_id,
     &            TRIM(svzd_blk_names(i)),svzd_blk_name_len(i),
     &                                          times,params,2)
            END DO

            ! Run zeroD simulation
            CALL lpn_interface_run_simulation(model_id, svZeroDTime, 
     &                         lpn_times, lpn_solutions, error_code)

            ! Extract pressure and flow from zeroD solution
            lpn_state_y = lpn_solutions((num_output_steps-1)*system_size
     &                                   :num_output_steps*system_size)
            DO i=1, numCoupledSrfs
               IF (i .LE. numDirichletSrfs) THEN
                  QCoupled(i) = lpn_state_y(sol_IDs(2*(i-1)+1))
               ELSE
                  PCoupled(i) = lpn_state_y(sol_IDs(2*(i-1)+2))
               END IF
            END DO

            IF (svzerodFlag .EQ. 'L') THEN !Last iteration
               ! Save state and update time only after last iteration
               CALL lpn_interface_return_ydot_(model_id,
     &                                  last_state_ydot)
               last_state_y = lpn_state_y
               ! Keep track of current time
               svZeroDTime = svZeroDTime + Delt(1)
            END IF
         ENDIF


         ! TODO
         ! Figure out initial condition
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

      END IF
      
      i = MAXSURF + 1
      CALL MPI_BCAST(QCoupled, i, MPI_DOUBLE_PRECISION, master,
     2   MPI_COMM_WORLD, ierr)
      
      CALL MPI_BCAST(PCoupled, i, MPI_DOUBLE_PRECISION, master,
     2   MPI_COMM_WORLD, ierr)
      
      RETURN
      END SUBROUTINE calcsvZeroD

!> Main subroutine for calculating the surface pressure 
!! derivative with respect to flow rate

      SUBROUTINE calcsvZeroDBCDerivative
       
      USE svZeroD
      INCLUDE "global.h"
      INCLUDE "mpif.h"
      INCLUDE "common_blocks/workfc.h"
      INCLUDE "common_blocks/inpdat.h"
      INCLUDE "common_blocks/nomodule.h"

      INTEGER MPIstat(MPI_STATUS_SIZE)
      REAL*8 diff, PBase(numCoupledSrfs), garbage
      REAL*8, PARAMETER :: absTol = 1D-8, relTol = 1D-5
      INTEGER ierr,i,j
      INTEGER error_code

      REAL*8 :: params(2), times(2)

      IF (numNeumannSrfs .EQ. 0) RETURN
      
      PDerFlag = .FALSE.
      !i0 = numDirichletSrfs

      diff = SUM(ABS(QCoupled(numDirichletSrfs+1:numCoupledSrfs)))
     2   /REAL(numNeumannSrfs,8)

      IF (diff*relTol .LT. absTol) THEN
         diff = absTol
      ELSE
         diff = diff*relTol
      END IF

      IF (myrank .EQ. master) THEN
         IF (svZeroDTime > 0.0D0) THEN
            ! Set initial condition from previous state
            CALL lpn_interface_update_state(model_id, last_state_y,
     &                                             last_state_ydot)
         END IF

         times = (/svZeroDTime, svZeroDTime+Delt(1)/)
         DO j = numDirichletSrfs, numCoupledSrfs
            ! Update pressure and flow in the zeroD model
            DO i=1, numCoupledSrfs
               IF (i .LE. numDirichletSrfs) THEN
                  params = (/PnCoupled(i), PCoupled(i)/)
               ELSE
                  IF (i .NE. j) THEN
                     params = (/QnCoupled(i), QCoupled(i)/)
                  ELSE
                     params = (/QnCoupled(i), QCoupled(i) + diff/)
                  END IF
               END IF
               CALL lpn_interface_update_block_params(model_id,
     &            TRIM(svzd_blk_names(i)),svzd_blk_name_len(i),
     &                                          times,params,2)
            END DO

            ! Run zeroD simulation
            CALL lpn_interface_run_simulation(model_id, svZeroDTime, 
     &                         lpn_times, lpn_solutions, error_code)

            ! Extract pressure and flow from zeroD solution
            lpn_state_y = lpn_solutions((num_output_steps-1)*system_size
     &                                   :num_output_steps*system_size)
            DO i=1, numCoupledSrfs
               IF (i > numDirichletSrfs) THEN
                  IF (i .NE. j) THEN
                     IF (j .EQ. numDirichletSrfs) THEN
                        PBase(i) = lpn_state_y(sol_IDs(2*(i-1)+2))
                     END IF
                  ELSE
                     PDer(i) = lpn_state_y(sol_IDs(2*(i-1)+2))
                     PDer(i) = (PDer(i) - PBase(i))/diff
                  ENDIF
               END IF
            END DO
         END DO

!        IF (svzerodFlag .EQ. 'L') THEN !Last iteration
!           ! Save state and update time only after last iteration
!           CALL lpn_interface_return_ydot_(model_id,
!    &                               last_state_ydot)
!           last_state_y = lpn_state_y
!           ! Keep track of current time
!           svZeroDTime = svZeroDTime + Delt(1)
!        END IF

!        DO j=i0, numCoupledSrfs
!           OPEN (1, FILE='GenBC.int', STATUS='UNKNOWN', 
!    2         FORM='UNFORMATTED')
!           WRITE (1) 'D'
!           WRITE (1) Delt(1)
!           WRITE (1) numDirichletSrfs
!           WRITE (1) numNeumannSrfs
!           DO i=1, numCoupledSrfs
!              IF (i .LE. i0) THEN
!                 WRITE (1) PnCoupled(i), PCoupled(i)
!              ELSE
!                 IF (i .NE. j) THEN
!                    WRITE (1) QnCoupled(i), QCoupled(i)
!                 ELSE
!                    WRITE (1) QnCoupled(i), QCoupled(i) + diff
!                 END IF
!              END IF
!           END DO
!           CLOSE (1)
!           
!           IF (iGenFromFile .EQ. 1) THEN
!              CALL system('./GenBC')
!           ELSE
!              CALL system('GenBC')
!           END IF

!           OPEN (1,FILE='GenBC.int',STATUS='OLD', FORM='UNFORMATTED')
!           DO i=1, numCoupledSrfs
!              IF (i .LE. i0) THEN
!                 READ (1) garbage
!              ELSE
!                 IF (i .NE. j) THEN
!                    IF (j .EQ. i0) THEN
!                       READ (1) PBase(i)
!                    ELSE
!                       READ (1) garbage
!                    END IF
!                 ELSE
!                    READ (1) PDer(i)
!                    PDer(i) = (PDer(i) - PBase(i))/diff
!                 END IF
!              END IF
!           END DO
!           CLOSE(1)
!        END DO
      END IF

      i = MAXSURF + 1
      CALL MPI_BCAST(PDer, i, MPI_DOUBLE_PRECISION, master,
     2   MPI_COMM_WORLD, ierr)

      RETURN
      END SUBROUTINE calcsvZeroDBCDerivative



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
            OPEN (113, FILE=TRIM(fileNames(i)), STATUS='OLD', 
     2         ACCESS='APPEND')

            WRITE (113,fmt=myFMT1) R(i,:)
            CLOSE (113)
         ELSE
            OPEN (113, FILE=TRIM(fileNames(i)), STATUS='NEW')
            WRITE (113,fmt=myFMT2) surfID
            WRITE (113,fmt=myFMT1) R(i,:)
            CLOSE (113)
         END IF
      END DO

      RETURN
      END SUBROUTINE printSvZeroD

#endif
