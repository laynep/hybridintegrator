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


program hybrid_integrator_d
use d_hybrid_initialconditions
use hybrid_subroutines
use rng
use mpi
use linked_list
implicit none

	integer:: i,j, points, success, counter, iccounter, sucunit,&
		& failunit, failcount, localcount
	integer :: errorcount, iend
	integer :: badfieldcounter, badfieldlocal, successlocal, faillocal,&
		& errorlocal, ierr, rc
	integer :: numtasks, rank
	integer :: ic, trajnumb
	double precision :: check, v, ratio, dt
	double precision :: v_0
	double precision, dimension(:,:), allocatable :: sample_table
	logical :: leave, allfailcheck, printing, traj
	type(llnode), pointer :: ytraj_head, ytraj_tail

	!*****************************
	!FCVODE PARAMS
	double precision :: y(5), rpar(5)
	integer :: ipar(5), meth, itmeth
	real ::  rout(6)
	integer(kind=8) :: neq, nglobal
	integer :: ier, iatol, iout(21), itask
	double precision :: t0, y0(5), t, tout, rtol, atol(5)
	!*****************************

	namelist /ics/ points, IC, dt, printing, traj

	!Read numb of (data points per numb of processes) & IC type from file.
	!Do we want to print to stdout?  Do we want to record the trajectory?
	open(unit=10000, file="parameters_hybrid.txt", status="old", delim = "apostrophe")
	read(unit=10000, nml=ics)
	close(unit=10000)

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
		if(ierr .ne. MPI_SUCCESS) then
			if(printing) print*,"Error parallelizing."
			call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
			stop
		end if
	!Obtains info on processors.
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
	if(printing) print*,'Number of tasks=',numtasks,' My rank=',rank
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	!Opens success and fail files. Optionally traj files.  Defaults to unform binary.
	call open_hybridfiles(rank,numtasks,sucunit,failunit)
	if (traj) call open_trajfiles(rank, trajnumb)

	!Set potential parameters.
	call parameters_hybrid()

	!Print stats.
	if(rank==0) call hybrid_initstats(ic, printing)

	!Set seed for each thread from module rng.
	call init_random_seed(rank)

	!Set some params for FCVODE integrator.
	call set_paramsFCVODE(rpar, neq, nglobal, numtasks, iatol, atol, rtol, &
		& meth, itmeth, t0, t, itask, tout)

	!If recording trajs, then initialize linked list.
	if (traj) call ll_init(ytraj_head,ytraj_tail)

	!Set IC.
	if (IC == 1) then
		call D_IC_ZEROV(Y0)
	else if (IC == 2) then
		call D_IC_EQEN(Y0,iccounter)
	else if (IC==3) then
		call EQEN_SLICING(Y0)
	else if (IC==4) then
		call IC_METR_INIT(Y0, iccounter, sample_table, 10000)
	end if
	Y=Y0

	!Initialize FCVODE integrator.
	call FCVMALLOC(T0, Y0, METH, ITMETH, IATOL, RTOL, ATOL,&
		&IOUT, ROUT, IPAR, RPAR, IER)
	call FCVSETIIN("MAX_NSTEPS", 5000000, IER)
	call FCVDENSE(NEQ, IER)
	call FCVDENSESETJAC (1, IER)


	!Loop over ICs until achieve numb of desired points.
do1: 	do while (successlocal<points) 
		localcount = localcount + 1
		!Get new point if on second or greater run.
		if (localcount>1) then
			call new_point(y0,iccounter,sample_table, ic)
			!Reinit time and Y
			Y=Y0
			T0=0D0
			TOUT = 1D1
			T=T0
			!Reinitialize integrator.
			ITASK = 1
			call FCVREINIT(T0, Y0, IATOL, RTOL, ATOL, IER)			
		end if
		
		!Counters
		success = 0
		iccounter = iccounter + 1
	
		!Perform the integration.
		iend=3000000
	do3:	do i=1,iend
			!Take field values if recording trajectory. Previously initialized
			if (traj) call rec_traj(Y, ytraj_head, ytraj_tail)

			!*********************************
			!Perform the integration. dt set in namelist ics.
			call FCVODE(TOUT,T,Y,ITASK,IER)
			TOUT = TOUT + dt
			!*********************************

			!Check succ or fail condition: N>65 success=1, and N<65 failure=0.
			call succ_or_fail(Y, success, successlocal, faillocal,&
				& sucunit, failunit, ic, leave, printing)
			!If condition met, leave=.true.
			if (leave) exit do3
			
			!Tells if not enough time for fields to evolve.
			if (i==iend) then
				errorlocal = errorlocal + 1
				if (printing) print*, "Error: field didn't reach minima."
			end if
			if(printing .and. MOD(i,iend/10)==0) print*,"i is getting big...",i
		end do do3
		
		!Print the traj & delete O(2n).
		if (traj) call print_del_traj(ytraj_head, ytraj_tail, trajnumb)

		!Check if the integrator isn't finding any succ points.
		call all_fail_check(successlocal, faillocal, allfailcheck, printing)
		if (allfailcheck) exit do1

	end do do1

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

	!Clean up integrator.
	call FCVFREE

	!End parallel.
	call MPI_FINALIZE(ierr)

	!Why this is necessary I have no idea.  It works fine without it on my computer, but MPI gives an error if this isn't here on the cluster.
	stop

end program hybrid_integrator_d




!***********************************************************************
!RHS of equation to integrate.
!***********************************************************************
subroutine FCVFUN(T, Y, YDOT, IPAR, RPAR, IER)
implicit none
	
	!******************
	double precision, intent(in) :: Y(*), t
	double precision, intent(inout) :: ydot(*)
	integer, intent(in) :: IPAR(*)
	integer, intent(out) :: IER
	double precision, intent(in) :: RPAR(*)
	!******************

	double precision :: V, D_phi_V, D_psi_V, DD_phi_V, DD_psi_V, D_phi_D_psi_V, Hub


	!Potential and its derivatives.
	V = (RPAR(1)*RPAR(1)*RPAR(1)*RPAR(1))*((1D0 - &
		&((Y(3)*Y(3))/(RPAR(2)*RPAR(2))))**2D0 &
		&+ ((Y(2)*Y(2))/(RPAR(3)*RPAR(3)))&
		& + ((Y(2)*Y(2)*Y(3)*Y(3))/ (RPAR(4)**4D0)))

	D_phi_V = 2D0*(RPAR(1)*RPAR(1)*RPAR(1)*RPAR(1))*(((Y(2))/(RPAR(3)*RPAR(3)))+ &
		&((Y(2)*Y(3)*Y(3))/(RPAR(4)**4D0)))

	D_psi_V = 2D0*(RPAR(1)*RPAR(1)*RPAR(1)*RPAR(1))*(((-2D0*Y(3))/&
		&(RPAR(2)*RPAR(2)))*(1D0 &
		&-((Y(3)*Y(3))/(RPAR(2)*RPAR(2))))+&
		& ((Y(2)*Y(2)*Y(3))/(RPAR(4)**4D0)))

	DD_phi_V = 2D0*(RPAR(1)*RPAR(1)*RPAR(1)*RPAR(1))*((1D0/(RPAR(3)*RPAR(3)))+&
		&((Y(3)*Y(3))/(RPAR(4)**4D0)))

	DD_psi_V = 2D0*(RPAR(1)*RPAR(1)*RPAR(1)*RPAR(1))*(((-2D0)/(RPAR(2)*RPAR(2)))+ &
		&((4D0*Y(3)*Y(3))/(RPAR(2)**4D0))+&
		& ((Y(2)*Y(2))/(RPAR(4)**4D0)))

	D_phi_D_psi_V = (4D0*(RPAR(1)*RPAR(1)*RPAR(1)*RPAR(1))*Y(2)*Y(3))/&
		&(RPAR(4)*RPAR(4)*RPAR(4)*RPAR(4))

	Hub = sqrt(RPAR(5)*(.5D0*((Y(4)*Y(4)) + (Y(5)*Y(5))) + V))

	
 	!Equations of motion.
	YDOT(1) = Hub
	YDOT(2) = Y(4)
	YDOT(3) = Y(5)
	YDOT(4) = -3D0*Hub*Y(4)-D_phi_V
	YDOT(5) = -3D0*Hub*Y(5)-D_psi_V

	!Success
	IER = 0

end subroutine FCVFUN

!***************************************************************************
!Jacobian of RHS of equation to integrate.
subroutine FCVDJAC (NEQ, T, Y, FY, DJAC, H, IPAR, RPAR,&
			&WK1, WK2, WK3, IER)
implicit none

	!**************************
	double precision, intent(in) :: Y(*), FY(*), T, H
	integer, intent(in) :: IPAR(*), NEQ
	integer, intent(out) :: IER
	double precision, intent(inout) :: DJAC(NEQ,*)
	double precision, intent(in) :: RPAR(*)
	double precision, intent(in) :: WK1(*), WK2(*), WK3(*)
	!**************************

	double precision :: V, D_phi_V, D_psi_V, DD_phi_V, DD_psi_V, D_phi_D_psi_V, Hub

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

end subroutine FCVDJAC








