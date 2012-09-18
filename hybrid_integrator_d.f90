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
use types, only : dp
implicit none

	!Main variables.
	real(dp), dimension(5) :: y, y0		!The fields.
	real(dp), dimension(5) :: yref		!Ref field for zooming proc.
	real(dp) :: t0, t, tout, dt		!Timers.
	!Counters & unit numbers
	integer:: i, success, counter, iccounter, sucunit,&
		& failunit, failcount, localcount
	integer :: errorcount, iend
	integer :: badfieldcounter, badfieldlocal, successlocal, faillocal,&
		& errorlocal, ierr, rc
	!Program variables.
	integer :: ic, trajnumb, points
	real(dp) :: check, v, ratio, toler
	logical :: leave, allfailcheck, printing, traj
	logical :: integr_ch
	!Variables to load IC from file for direct read or interpolation (ic=4,5)
	real(dp), dimension(:,:), allocatable :: sample_table, ic_table
	!List to record trajectory.
	type(linkedlist) :: ytraj
	!Parallel variables.
	integer :: numtasks, rank
	!FCVODE params
	real(dp) :: rpar(5)
	integer :: ipar(5), meth, itmeth
	real ::  rout(6)
	integer(kind=8) :: neq, nglobal
	integer :: ier, iatol, iout(21), itask
	real(dp) :: rtol, atol(5)

	namelist /ics/ points, IC, dt, printing, traj

	!***********************************************************************
	!***********************************************************************

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
		& meth, itmeth, t0, t, itask, tout, dt)

	!Set first IC.
	if (IC == 1) then
		!Zero vel slice.
		call D_IC_ZEROV(Y0)
	else if (IC == 2) then
		!Eq energy slice.
		call D_IC_EQEN(Y0,iccounter)
	else if (IC==3) then
		!Subslice of eq en slice.
		call EQEN_SLICING(Y0)
	else if (IC==4) then
		!Metropolis sample a dataset.
		call IC_METR_INIT(Y0, iccounter, sample_table, 10000)
	else if (IC==5) then
		!Read IC from a file.
		call ic_file_init(y0, rank,numtasks,ic_table)
	else if (IC==6) then
		!Zoom in on one point on eq en surface.
		call ic_zoom_init(y0, yref, iccounter, toler)
  else if (ic==7) then
    !Get a fixed initial condition.
    call fixed_ic(y0)
	end if
	Y=Y0

	!Initialize FCVODE integrator.
	call FCVMALLOC(T0, Y0, METH, ITMETH, IATOL, RTOL, ATOL,&
		&IOUT, ROUT, IPAR, RPAR, IER)
	call FCVSETIIN("MAX_NSTEPS", 5000000, IER)
	call FCVDENSE(NEQ, IER)
	call FCVDENSESETJAC (1, IER)

	!Loop over ICs until achieve numb of desired points.
	integr_ch=.true.	!Exit condition.
do1: 	do while (integr_ch)

		!Count numb of times each thread goes through loop.
		localcount = localcount + 1
		!Get new point if on second or greater run.
		if (localcount>1) then
			call new_point(y0,iccounter,sample_table, ic, &
				&ic_table, yref, toler)

			!Reinit time and Y
			Y=Y0
			T0=0_dp
			TOUT = dt
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
			if (traj) call rec_traj(Y, ytraj)


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

		!Print the traj & delete -- O(2n).
		if (success==1 .and. traj) then
      call print_del_traj(ytraj, trajnumb)
    !If not successful, then just delete list.
    else if (success==0 .and. traj) then
      call ll_del_all(ytraj)
    end if

		!Check if the integrator isn't finding any succ points.
		call all_fail_check(successlocal, faillocal, allfailcheck, printing, check)
		if (allfailcheck) exit do1
		!Determine loop exit condition.
		if (ic<5 .or. ic==6) then
			integr_ch=(successlocal<points)
		else if (ic==5) then
			integr_ch=(localcount<size(ic_table,1))
		end if

	end do do1

	!Halts processors here.
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!Gives child data to master.
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
	if(rank==0) then
		call hybrid_finalstats(ic, counter, failcount, &
		&badfieldcounter, errorcount, printing)
	end if

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
  use types, only : dp
  implicit none
	
	!******************
	real(dp), intent(in) :: Y(*), t
	real(dp), intent(out) :: ydot(*)
	integer, intent(in) :: IPAR(*)
	integer, intent(out) :: IER
	real(dp), intent(in) :: RPAR(*)
	!******************

	real(dp) :: V, D_phi_V, D_psi_V, DD_phi_V, DD_psi_V, D_phi_D_psi_V, Hub


	!Potential and its derivatives.
	v = (rpar(1)*rpar(1)*rpar(1)*rpar(1))*((1_dp - &
		&((y(3)*y(3))/(rpar(2)*rpar(2))))**2 &
		&+ ((y(2)*y(2))/(rpar(3)*rpar(3)))&
		& + ((y(2)*y(2)*y(3)*y(3))/ (rpar(4)**4)))

	d_phi_v = 2_dp*(rpar(1)*rpar(1)*rpar(1)*rpar(1))*(((y(2))/(rpar(3)*rpar(3)))+ &
		&((y(2)*y(3)*y(3))/(rpar(4)**4)))

	d_psi_v = 2_dp*(rpar(1)*rpar(1)*rpar(1)*rpar(1))*(((-2_dp*y(3))/&
		&(rpar(2)*rpar(2)))*(1_dp &
		&-((y(3)*y(3))/(rpar(2)*rpar(2))))+&
		& ((y(2)*y(2)*y(3))/(rpar(4)**4)))

	dd_phi_v = 2_dp*(rpar(1)*rpar(1)*rpar(1)*rpar(1))*((1_dp/(rpar(3)*rpar(3)))+&
		&((y(3)*y(3))/(rpar(4)**4)))

	dd_psi_v = 2_dp*(rpar(1)*rpar(1)*rpar(1)*rpar(1))*(((-2_dp)/(rpar(2)*rpar(2)))+ &
		&((4_dp*y(3)*y(3))/(rpar(2)**4))+&
		& ((y(2)*y(2))/(rpar(4)**4)))

	d_phi_d_psi_v = (4_dp*(rpar(1)*rpar(1)*rpar(1)*rpar(1))*y(2)*y(3))/&
		&(rpar(4)*rpar(4)*rpar(4)*rpar(4))

	hub = sqrt(rpar(5)*(.5_dp*((y(4)*y(4)) + (y(5)*y(5))) + v))
	if (hub<0) then
		print*,"Hubble parameter < 0"
		ier=1
		return
	end if
	
 	!Equations of motion.
	ydot(1) = hub
	ydot(2) = y(4)
	ydot(3) = y(5)
	ydot(4) = -3_dp*hub*y(4)-d_phi_v
	ydot(5) = -3_dp*hub*y(5)-d_psi_v

	!Success
	IER = 0

end subroutine FCVFUN

!***************************************************************************
!Jacobian of RHS of equation to integrate.
subroutine FCVDJAC (NEQ, T, Y, FY, DJAC, H, IPAR, RPAR,&
			&WK1, WK2, WK3, IER)
  use types, only : dp
  implicit none

	!**************************
	real(dp), intent(in) :: Y(*), FY(*), T, H
	integer, intent(in) :: IPAR(*), NEQ
	integer, intent(out) :: IER
	real(dp), intent(inout) :: DJAC(NEQ,*)
	real(dp), intent(in) :: RPAR(*)
	real(dp), intent(in) :: WK1(*), WK2(*), WK3(*)
	!**************************

	real(dp) :: V, D_phi_V, D_psi_V, DD_phi_V, DD_psi_V, D_phi_D_psi_V, Hub

	!Potential, V, and its derivatives.
	v = (rpar(1)**4)*((1_dp - ((y(3)*y(3))/(rpar(2)*rpar(2))))**2 +&
		& ((y(2)*y(2))/(rpar(3)*rpar(3)))&
		& + ((y(2)*y(2)*y(3)*y(3))/ (rpar(4)**4)))

	d_phi_v = 2_dp*(rpar(1)**4)*(((y(2))/(rpar(3)*rpar(3)))+&
		& ((y(2)*y(3)*y(3))/(rpar(4)**4)))

	d_psi_v = 2_dp*(rpar(1)**4)*(((-2_dp*y(3))/(rpar(2)*rpar(2)))*(1_dp &
		&-((y(3)*y(3))/(rpar(2)*rpar(2))))+ &
		&((y(2)*y(2)*y(3))/(rpar(4)**4)))

	dd_phi_v = 2_dp*(rpar(1)**4)*((1_dp/(rpar(3)*rpar(3)))+&
		&((y(3)*y(3))/(rpar(4)**4)))

	dd_psi_v = 2_dp*(rpar(1)**4)*(((-2_dp)/(rpar(2)*rpar(2)))+ &
		&((4_dp*y(3)*y(3))/(rpar(2)**4))+&
		& ((y(2)*y(2))/(rpar(4)**4)))

	d_phi_d_psi_v = (4_dp*(rpar(1)**4)*y(2)*y(3))/(rpar(4)**4)

	hub = sqrt(rpar(5)*(.5_dp*((y(4)*y(4)) + (y(5)*y(5))) + v))

	!Partial derivatives of the RHS vector in system of equations.
	djac(1,1)= 0_dp
	djac(1,2)= .5_dp*(1_dp/hub)*rpar(5)*d_phi_v
	djac(1,3)= .5_dp*(1_dp/hub)*rpar(5)*d_psi_v
	djac(1,4)= .5_dp*(1_dp/hub)*rpar(5)*y(4)
	djac(1,5)= .5_dp*(1_dp/hub)*rpar(5)*y(5)
	djac(2,1)= 0_dp
	djac(2,2)= 0_dp
	djac(2,3)= 0_dp
	djac(2,4)= 1_dp
	djac(2,5)= 0_dp
	djac(3,1)= 0_dp
	djac(3,2)= 0_dp
	djac(3,3)= 0_dp
	djac(3,4)= 0_dp
	djac(3,5)= 1_dp
	djac(4,1)= 0_dp
	djac(4,2)= -1.5_dp*(1_dp/hub)*rpar(5)*y(4)*d_phi_v - dd_phi_v
	djac(4,3)= -1.5_dp*(1_dp/hub)*rpar(5)*y(4)*d_psi_v - d_phi_d_psi_v
	djac(4,4)= -1.5_dp*(1_dp/hub)*rpar(5)*y(4)*y(4) - 3_dp*hub
	djac(4,5)= -1.5_dp*(1_dp/hub)*rpar(5)*y(4)*y(5)
	djac(5,1)= 0_dp
	djac(5,2)= -1.5_dp*(1_dp/hub)*rpar(5)*y(5)*d_phi_v - d_phi_d_psi_v
	djac(5,3)= -1.5_dp*(1_dp/hub)*rpar(5)*y(5)*d_psi_v - dd_psi_v
	djac(5,4)= -1.5_dp*(1_dp/hub)*rpar(5)*y(4)*y(5)
	djac(5,5)= -1.5_dp*(1_dp/hub)*rpar(5)*y(5)*y(5) - 3_dp*hub

	!Success
	IER = 0

end subroutine FCVDJAC








