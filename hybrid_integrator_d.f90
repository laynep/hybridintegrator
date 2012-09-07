!*****************************************************************
!Layne Price, University of Auckland, May 15, 2012.

!This is a program that will probe the space of initial conditions for 2 field hybrid inflation in order to
!determine which initial conditions will yield more than 60 e-folds of inflation.  This
!will be done by solving the dynamical equations for a number of initial conditions and
!then plotting the results.  The equations will be solved by calling a differential
!equation solver called FCVODE.  FCVODE will solve the system from time t=T to t=TOUT
!and spits out the values of Y, so it's necessary to loop this over all required values
!of TOUT. 

!Y(1)=N, Y(2)=phi, Y(3)=psi, Y(4)=phi_dot, Y(5)=psi_dot.
!*****************************************************************


PROGRAM hybrid_integrator_d
USE d_hybrid_initialconditions
use hybrid_subroutines
use mpi
IMPLICIT NONE

	INTEGER:: i,j, points, success, counter, iccounter, sucunit,&
		& failunit, failcount, localcount
	INTEGER :: errorcount, iend
	INTEGER :: badfieldcounter, badfieldlocal, successlocal, faillocal,&
		& errorlocal, ierr, rc
	INTEGER :: numtasks, rank
	INTEGER :: IC
	DOUBLE PRECISION :: check, V, ratio, dt
	DOUBLE PRECISION :: V_0
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: sample_table
	INTEGER :: samp_len, samp_wid
	logical :: leave, allfailcheck, printing

	!*****************************
	!FCVODE PARAMS
	DOUBLE PRECISION :: Y(5), RPAR(5)
	INTEGER :: IPAR(5), METH, ITMETH
	REAL ::  ROUT(6)
	INTEGER(KIND=8) :: NEQ, NGLOBAL
	INTEGER :: IER, IATOL, IOUT(21), ITASK
	DOUBLE PRECISION :: T0, Y0(5), T, TOUT, RTOL, ATOL(5)
	!*****************************

	NAMELIST /ics/ points, IC
	NAMELIST /sample/ samp_len, samp_wid

	!Read numb of data points per numb of processes & IC type from file.
	OPEN(unit=10000, file="parameters_hybrid.txt", status="old", delim = "apostrophe")
	READ(unit=10000, nml=ics)
	CLOSE(unit=10000)
!	points = 100000

!	IC = 1		!IC with zero vel slice.
!	IC = 2		!IC on eq energy slice.
!	IC = 3		!IC as a slicing of EQEN slice.
!	IC = 4		!IC from Metropolis Algorithm.

	!Do we want to print to stdout?  Excludes some error messages...
	printing = .true.

	!Global counters.
	counter = 0
	failcount = 0
	badfieldcounter=0
	iccounter = 0
	errorcount = 0
	!Local counters.
	errorlocal = 0
	badfieldlocal = 0
	successlocal = 0
	faillocal = 0
	localcount = 0

	!Parallelizes.
	call MPI_INIT(ierr)
		IF(ierr .NE. MPI_SUCCESS) THEN
			if(printing) print*,"Error parallelizing."
			call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
			stop
		END IF
	!Obtains info on processors.
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
	if(printing) print*,'Number of tasks=',numtasks,' My rank=',rank
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	!Opens success and fail files.
	call open_hybridfiles(rank,numtasks,sucunit,failunit)

	!Set potential parameters.
	call parameters_hybrid()

	!Print stats.
	if(rank==0) call hybrid_initstats(ic, printing)

	!Set seed for each thread.
	call init_random_seed1(rank)

	!Set some params for FCVODE integrator.
	call set_paramsFCVODE(rpar, neq, nglobal, numtasks, iatol, atol, rtol, &
		& meth, itmeth, t0, t, itask, tout)

	!Set IC.
	IF (IC == 1) THEN
		call D_IC_ZEROV(Y0)
	ELSE IF (IC == 2) THEN
		call D_IC_EQEN(Y0,iccounter)
	ELSE IF (IC==3) THEN
		call EQEN_SLICING(Y0)
	ELSE IF (IC==4) THEN
		!Set IC at random on EQEN slice.
		call D_IC_EQEN(Y0,iccounter)
		!Load the data which we will later sample via nearest neighbor interpolation.
		OPEN(unit=10000, file="parameters_hybrid.txt", status="old",&
		& delim = "apostrophe")
		READ(unit=10000, nml=sample)
		CLOSE(unit=10000)
		
		ALLOCATE(sample_table(samp_len,samp_wid))
		sample_table=0D0

		!Data *must* be in form Y(2),...,Y(5), where Y(1)=0D0 assumed.
		OPEN(unit=10001, file="datasample.bin", status="old", form="UNFORMATTED")
		DO i=1,SIZE(sample_table,1)		
			READ(unit=10001,END=10) (sample_table(i,j),j=2,5)
		END DO
10		CLOSE(unit=10001)

		!Burn in.
		DO i=1, 10000
			!Get new IC from sample_table.
			call IC_METR(Y0,sample_table,iccounter)
		END DO
	END IF
	Y=Y0

	!Initialize FCVODE integrator.
	call FCVMALLOC(T0, Y0, METH, ITMETH, IATOL, RTOL, ATOL,&
		&IOUT, ROUT, IPAR, RPAR, IER)
	call FCVSETIIN("MAX_NSTEPS", 5000000, IER)
	call FCVDENSE(NEQ, IER)
	call FCVDENSESETJAC (1, IER)




	!Loop over ICs until achieve numb of desired points.
do1: 	DO WHILE (successlocal<points) 

		localcount = localcount + 1
		!Get new point if on second or greater run.
		IF (localcount>1) THEN
			call new_point(y0,iccounter,sample_table, ic)
			
			Y=Y0
			T0=0D0
			TOUT = 1D1
			T=T0
			ITASK = 1
			!Reinitialize integrator.
			CALL FCVREINIT(T0, Y0, IATOL, RTOL, ATOL, IER)			
		END IF
		
		success = 0
		iccounter = iccounter + 1
	
		!Perform the integration.
		iend=3000000
	do3:	DO i=1,iend

			!*********************************
			!Perform the integration.
			CALL FCVODE(TOUT,T,Y,ITASK,IER)
			dt=1D4
			TOUT = TOUT + dt
			!*********************************

			!Check succ or fail condition: N>65 success=1, and N<65 failure=0.
			call succ_or_fail(Y, success, successlocal, faillocal,&
				& sucunit, failunit, ic, leave, printing)
			!If condition met, leave=.true.
			if (leave) exit do3
			
			!Tells if not enough time for fields to evolve.
			IF (i==iend) THEN
				errorlocal = errorlocal + 1
				if (printing) print*, "Error: field didn't reach minima."
			END IF
			IF(printing .and. MOD(i,iend/10)==0) print*,"i is getting big...",i
		END DO do3

		!Check if the integrator isn't finding any succ points.
		call all_fail_check(successlocal, faillocal, allfailcheck, printing)
		if (allfailcheck) exit do1


	END DO do1

	!Halts processors here.
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	!Gives slave data to master.
	call MPI_REDUCE(badfieldlocal,badfieldcounter,1,MPI_INTEGER,&
		&MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(successlocal,counter,1,MPI_INTEGER,&
		&MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(faillocal,failcount,1,MPI_INTEGER,&
		&MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(errorlocal,errorcount,1,MPI_INTEGER,MPI_SUM,&
		&0,MPI_COMM_WORLD,ierr)

	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	!Print from master.
	if(rank==0) call hybrid_finalstats(ic, counter, failcount, &
			&badfieldcounter, errorcount, printing)

	!Clean up.
	call FCVFREE

	!End parallel.
	call MPI_FINALIZE(ierr)

	!Why this is necessary I have no idea.  It works fine without it on my computer, but MPI gives an error if this isn't here on the cluster.
	stop

END PROGRAM hybrid_integrator_d




!**************************************************************************************


SUBROUTINE init_random_seed1(rank)
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
		INTEGER, INTENT(IN) :: rank
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37*rank* (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
END SUBROUTINE


!***********************************************************************
!RHS
!***********************************************************************
SUBROUTINE FCVFUN(T, Y, YDOT, IPAR, RPAR, IER)
IMPLICIT NONE
	
	!******************
	DOUBLE PRECISION :: Y(*), YDOT(*), T
	INTEGER :: IPAR(*), IER
	DOUBLE PRECISION :: RPAR(*)
	!******************

	DOUBLE PRECISION :: V, D_phi_V, D_psi_V, DD_phi_V, DD_psi_V, D_phi_D_psi_V, Hub


	!Potential and its derivatives.
	V = (RPAR(1)**4D0)*((1D0 - ((Y(3)*Y(3))/(RPAR(2)*RPAR(2))))**2D0 &
		&+ ((Y(2)*Y(2))/(RPAR(3)*RPAR(3)))&
		& + ((Y(2)*Y(2)*Y(3)*Y(3))/ (RPAR(4)**4D0)))

	D_phi_V = 2D0*(RPAR(1)**4D0)*(((Y(2))/(RPAR(3)*RPAR(3)))+ &
		&((Y(2)*Y(3)*Y(3))/(RPAR(4)**4D0)))

	D_psi_V = 2D0*(RPAR(1)**4D0)*(((-2D0*Y(3))/(RPAR(2)*RPAR(2)))*(1D0 &
		&-((Y(3)*Y(3))/(RPAR(2)*RPAR(2))))+&
		& ((Y(2)*Y(2)*Y(3))/(RPAR(4)**4D0)))

	DD_phi_V = 2D0*(RPAR(1)**4D0)*((1D0/(RPAR(3)*RPAR(3)))+&
		&((Y(3)*Y(3))/(RPAR(4)**4D0)))

	DD_psi_V = 2D0*(RPAR(1)**4D0)*(((-2D0)/(RPAR(2)*RPAR(2)))+ &
		&((4D0*Y(3)*Y(3))/(RPAR(2)**4D0))+&
		& ((Y(2)*Y(2))/(RPAR(4)**4D0)))

	D_phi_D_psi_V = (4D0*(RPAR(1)**4D0)*Y(2)*Y(3))/(RPAR(4)**4D0)

	Hub = SQRT(RPAR(5)*(.5D0*((Y(4)*Y(4)) + (Y(5)*Y(5))) + V))

	
 	!Equations of motion.
	YDOT(1) = Hub
	YDOT(2) = Y(4)
	YDOT(3) = Y(5)
	YDOT(4) = -3D0*Hub*Y(4)-D_phi_V
	YDOT(5) = -3D0*Hub*Y(5)-D_psi_V

	!Success
	IER = 0

END SUBROUTINE FCVFUN

!***************************************************************************
SUBROUTINE FCVDJAC (NEQ, T, Y, FY, DJAC, H, IPAR, RPAR,&
			&WK1, WK2, WK3, IER)
IMPLICIT NONE

	!**************************
	DOUBLE PRECISION :: Y(*), FY(*), DJAC(NEQ,*), T, H
	INTEGER :: IPAR(*), NEQ, IER
	DOUBLE PRECISION :: RPAR(*)
	DOUBLE PRECISION :: WK1(*), WK2(*), WK3(*)
	!**************************

	DOUBLE PRECISION :: V, D_phi_V, D_psi_V, DD_phi_V, DD_psi_V, D_phi_D_psi_V, Hub

	!Potential, V, and its derivatives.
	V = (RPAR(1)**4D0)*((1D0 - ((Y(3)*Y(3))/(RPAR(2)*RPAR(2))))**2D0 +&
		& ((Y(2)*Y(2))/(RPAR(3)*RPAR(3)))&
		& + ((Y(2)*Y(2)*Y(3)*Y(3))/ (RPAR(4)**4D0)))

	D_phi_V = 2D0*(RPAR(1)**4D0)*(((Y(2))/(RPAR(3)*RPAR(3)))+&
		& ((Y(2)*Y(3)*Y(3))/(RPAR(4)**4D0)))

	D_psi_V = 2D0*(RPAR(1)**4D0)*(((-2D0*Y(3))/(RPAR(2)*RPAR(2)))*(1D0 &
		&-((Y(3)*Y(3))/(RPAR(2)*RPAR(2))))+ &
		&((Y(2)*Y(2)*Y(3))/(RPAR(4)**4D0)))

	DD_phi_V = 2D0*(RPAR(1)**4D0)*((1D0/(RPAR(3)*RPAR(3)))+&
		&((Y(3)*Y(3))/(RPAR(4)**4D0)))

	DD_psi_V = 2D0*(RPAR(1)**4D0)*(((-2D0)/(RPAR(2)*RPAR(2)))+ &
		&((4D0*Y(3)*Y(3))/(RPAR(2)**4D0))+&
		& ((Y(2)*Y(2))/(RPAR(4)**4D0)))

	D_phi_D_psi_V = (4D0*(RPAR(1)**4D0)*Y(2)*Y(3))/(RPAR(4)**4D0)

	Hub = SQRT(RPAR(5)*(.5D0*((Y(4)*Y(4)) + (Y(5)*Y(5))) + V))

	!Partial derivatives of the RHS vector in system of equations.
	DJAC(1,1)= 0D0
	DJAC(1,2)= .5D0*(1D0/Hub)*RPAR(5)*D_phi_V
	DJAC(1,3)= .5D0*(1D0/Hub)*RPAR(5)*D_psi_V
	DJAC(1,4)= .5D0*(1D0/Hub)*RPAR(5)*Y(4)
	DJAC(1,5)= .5D0*(1D0/Hub)*RPAR(5)*Y(5)
	DJAC(2,1)= 0D0
	DJAC(2,2)= 0D0
	DJAC(2,3)= 0D0
	DJAC(2,4)= 1D0
	DJAC(2,5)= 0D0
	DJAC(3,1)= 0D0
	DJAC(3,2)= 0D0
	DJAC(3,3)= 0D0
	DJAC(3,4)= 0D0
	DJAC(3,5)= 1D0
	DJAC(4,1)= 0D0
	DJAC(4,2)= -1.5D0*(1D0/Hub)*RPAR(5)*Y(4)*D_phi_V - DD_phi_V
	DJAC(4,3)= -1.5D0*(1D0/Hub)*RPAR(5)*Y(4)*D_psi_V - D_phi_D_psi_V
	DJAC(4,4)= -1.5D0*(1D0/Hub)*RPAR(5)*Y(4)*Y(4) - 3D0*Hub
	DJAC(4,5)= -1.5D0*(1D0/Hub)*RPAR(5)*Y(4)*Y(5)
	DJAC(5,1)= 0D0
	DJAC(5,2)= -1.5D0*(1D0/Hub)*RPAR(5)*Y(5)*D_phi_V - D_phi_D_psi_V
	DJAC(5,3)= -1.5D0*(1D0/Hub)*RPAR(5)*Y(5)*D_psi_V - DD_psi_V
	DJAC(5,4)= -1.5D0*(1D0/Hub)*RPAR(5)*Y(4)*Y(5)
	DJAC(5,5)= -1.5D0*(1D0/Hub)*RPAR(5)*Y(5)*Y(5) - 3D0*Hub

	!Success
	IER = 0

END SUBROUTINE FCVDJAC








