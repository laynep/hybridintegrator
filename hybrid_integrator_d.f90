!*******************************************************************************
!Layne Price, University of Auckland, May 15, 2012.
!*******************************************************************************

!SUMMARY:
!A program that does integration for two field hybrid inflation.  Uses the
!integrator FCVODE from the LLNL SUNDIALS package, which requires the RHS of the
!ODE to be expressed in the external subroutine FCVFUN and the Jacobian in the
!external subroutine FCVDJAC.  These are included below the main program.
!FCVODE functions are stored in a library.

!OUTLINE:
!The program architecture is as follows: We parallelize using OpenMPI; load ICs
!according to the method specified in the namelist; integrate the initial
!condition; sort the initial condition into a "success" array if it reaches N>65
!and into a "fail" array if it reaches the minimum of the potential without
!inflating.  A new IC is chosen and the integration is repeated until we get
!enough points.  Counters on the number of points found is collected on the
!master thread and stats are printed.

!OPTIONS:
!Options for the program are contained in the namelists in the parameter file,
!parameters_hybrid.txt.  Program options are contained in the namelist &ics: how
!many successful initial conditions do we
!want to find (points); what method should we use to sample the IC space (IC);
!what should the time-step be (dt); do we want to print to stdout (printing); do we
!want to record the trajectories (traj).  If we are getting ICs via Metropolis
!sampling of the IC space, the namelist &sample can be used to specify the
!dimensions of the array we should be sampling and the name of the file it is
!stored in.  Similarly, the namelist &filetoread provides similar information if
!we are reading ICs explicitly from a given file, without sampling.  The &zoom
!namelist gives one point which we are to zoom in on, providing a high
!resolution sample near the given point with tolerance specified.  Note that
!this point should be put in the namelist as (e-fold,psi,phi,psi_dot,phi_dot).
!The namelist &parameters specify the particle physics parameters for hybrid
!inflation.

!DEPENDENCIES:
!The program subroutines are contained in the module hybrid_subroutines and the
!routines that calculate the initial conditions are contained in the
!d_hybrid_initialconditions module.  Random number generation is done with the
!module rng.  MPI functions are called with the module mpi.  When storing
!trajectories we don't a priori know how many steps it will take to reach an end
!state, so we store these as a linked list, with type and methods defined in the
!module linked_list.  The types module defines the amount of working precision;
!and the newunit function from the features module gives us the ability to open
!a new file without worrying about whether the unit has already been used.  The
!libraries from the SUNDIALS package are necessary to do the integration.  The
!modules use sorting and location routines that are collected in the sorters
!module.

!NOTE:  The Y vector corresponds to:
!Y(1)=N, Y(2)=phi, Y(3)=psi, Y(4)=phi_dot, Y(5)=psi_dot.
!*******************************************************************************


program hybrid_integrator_d
  use d_hybrid_initialconditions
  use hybrid_subroutines
  use rng, only : init_random_seed
  use mpi
  use linked_list
  use types, only : dp
  use features, only : newunit
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
	integer :: ic, trajnumb, points, u
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
	open(unit=newunit(u), file="parameters_hybrid.txt", status="old", delim = "apostrophe")
	read(unit=u, nml=ics)
	close(unit=u)

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

	!Opens success and fail files. Optionally traj files.  Defaults to unformatted binary.
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
  if (.not. traj) then
	  call FCVSETIIN("MAX_NSTEPS", 5000000, IER)
  else
    !If recording trajectories, we want the integrator to return at every step
    !taken so that we can use the adaptive step size control in FCVODE to get
    !all important data. Set with itask=2.
    call FCVSETIIN("MAX_NSTEPS", 10, ier)
    itask=2
  end if
	call FCVDENSE(NEQ, IER)
	call FCVDENSESETJAC (1, IER)

	!Loop over ICs until achieve numb of desired points.
	integr_ch=.true.	!Exit condition.
icloop: 	do while (integr_ch)

		!Count numb of times each thread goes through loop.
		localcount = localcount + 1
		!Get new point if on second or greater run.
		if (localcount>1) then
			call new_point(y0,iccounter,sample_table, ic, &
				&ic_table, yref, toler)

			!Reinit time and Y
			Y=Y0
			T0=0e0_dp
      tout =t0+ dt
			T=T0
			!Reinitialize integrator.
      if (.not. traj) then
        itask = 1
      else
        itask=2
      end if
			call FCVREINIT(T0, Y0, IATOL, RTOL, ATOL, IER)			
		end if
		!Counters
		success = 0
		iccounter = iccounter + 1

    !###################################################################
		!Perform the integration.
    iend=3000000
intloop:	do i=1,iend

			!Take field values if recording trajectory. Previously initialized
			if (traj .and. mod(i,10)==0) then
        call rec_traj(Y, ytraj)
      end if

			!*********************************
			!Perform the integration. dt set in namelist ics.
			call FCVODE(TOUT,T,Y,ITASK,IER)
      tout = tout + dt
      !*********************************
			
      !Check succ or fail condition: N>65 success=1, and N<65 failure=0.
			call succ_or_fail(Y, success, successlocal, faillocal,&
				& sucunit, failunit, ic, leave, printing)
			!If condition met, leave=.true.
			if (leave) exit intloop
			
			!Tells if not enough time for fields to evolve.
			if (i==iend) then
				errorlocal = errorlocal + 1
				if (printing) print*, "Error: field didn't reach minima."
			end if
			if(printing .and. MOD(i,iend/10)==0 .and. .not. traj) print*,"i is getting big...",i
		end do intloop
    !###################################################################


		!Print the traj & delete -- O(2n).
		if (success==1 .and. traj) then
      call print_del_traj(ytraj, trajnumb)
    !If not successful, then just delete list.
    else if (success==0 .and. traj) then
      call ll_del_all(ytraj)
    end if

		!Check if the integrator isn't finding any succ points.
		call all_fail_check(successlocal, faillocal, allfailcheck, printing, check)
		if (allfailcheck) exit icloop
		!Determine loop exit condition.
		if (ic<5 .or. ic==6) then
			integr_ch=(successlocal<points)
		else if (ic==5) then
			integr_ch=(localcount<size(ic_table,1))
		end if

	end do icloop

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

	!Why this is necessary I have no idea.  It works fine without it on my computer,
  !but MPI gives an error if this isn't here on the cluster.
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
	v = (rpar(1)*rpar(1)*rpar(1)*rpar(1))*((1e0_dp - &
		&((y(3)*y(3))/(rpar(2)*rpar(2))))**2 &
		&+ ((y(2)*y(2))/(rpar(3)*rpar(3)))&
		& + ((y(2)*y(2)*y(3)*y(3))/ (rpar(4)**4)))

	d_phi_v = 2e0_dp*(rpar(1)*rpar(1)*rpar(1)*rpar(1))*(((y(2))/(rpar(3)*rpar(3)))+ &
		&((y(2)*y(3)*y(3))/(rpar(4)**4)))

	d_psi_v = 2e0_dp*(rpar(1)*rpar(1)*rpar(1)*rpar(1))*(((-2e0_dp*y(3))/&
		&(rpar(2)*rpar(2)))*(1e0_dp &
		&-((y(3)*y(3))/(rpar(2)*rpar(2))))+&
		& ((y(2)*y(2)*y(3))/(rpar(4)**4)))

	dd_phi_v = 2e0_dp*(rpar(1)*rpar(1)*rpar(1)*rpar(1))*((1e0_dp/(rpar(3)*rpar(3)))+&
		&((y(3)*y(3))/(rpar(4)**4)))

	dd_psi_v = 2e0_dp*(rpar(1)*rpar(1)*rpar(1)*rpar(1))*(((-2e0_dp)/(rpar(2)*rpar(2)))+ &
		&((4e0_dp*y(3)*y(3))/(rpar(2)**4))+&
		& ((y(2)*y(2))/(rpar(4)**4)))

	d_phi_d_psi_v = (4e0_dp*(rpar(1)*rpar(1)*rpar(1)*rpar(1))*y(2)*y(3))/&
		&(rpar(4)*rpar(4)*rpar(4)*rpar(4))

	hub = sqrt(rpar(5)*(.5e0_dp*((y(4)*y(4)) + (y(5)*y(5))) + v))
	if (hub<0) then
		print*,"Hubble parameter < 0"
		ier=1
		return
	end if
	
 	!Equations of motion.
	ydot(1) = hub
	ydot(2) = y(4)
	ydot(3) = y(5)
	ydot(4) = -3e0_dp*hub*y(4)-d_phi_v
	ydot(5) = -3e0_dp*hub*y(5)-d_psi_v

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
	v = (rpar(1)**4)*((1e0_dp - ((y(3)*y(3))/(rpar(2)*rpar(2))))**2 +&
		& ((y(2)*y(2))/(rpar(3)*rpar(3)))&
		& + ((y(2)*y(2)*y(3)*y(3))/ (rpar(4)**4)))

	d_phi_v = 2e0_dp*(rpar(1)**4)*(((y(2))/(rpar(3)*rpar(3)))+&
		& ((y(2)*y(3)*y(3))/(rpar(4)**4)))

	d_psi_v = 2e0_dp*(rpar(1)**4)*(((-2e0_dp*y(3))/(rpar(2)*rpar(2)))*(1e0_dp &
		&-((y(3)*y(3))/(rpar(2)*rpar(2))))+ &
		&((y(2)*y(2)*y(3))/(rpar(4)**4)))

	dd_phi_v = 2e0_dp*(rpar(1)**4)*((1e0_dp/(rpar(3)*rpar(3)))+&
		&((y(3)*y(3))/(rpar(4)**4)))

	dd_psi_v = 2e0_dp*(rpar(1)**4)*(((-2e0_dp)/(rpar(2)*rpar(2)))+ &
		&((4e0_dp*y(3)*y(3))/(rpar(2)**4))+&
		& ((y(2)*y(2))/(rpar(4)**4)))

	d_phi_d_psi_v = (4e0_dp*(rpar(1)**4)*y(2)*y(3))/(rpar(4)**4)

	hub = sqrt(rpar(5)*(.5e0_dp*((y(4)*y(4)) + (y(5)*y(5))) + v))

	!Partial derivatives of the RHS vector in system of equations.
	djac(1,1)= 0e0_dp
	djac(1,2)= .5e0_dp*(1e0_dp/hub)*rpar(5)*d_phi_v
	djac(1,3)= .5e0_dp*(1e0_dp/hub)*rpar(5)*d_psi_v
	djac(1,4)= .5e0_dp*(1e0_dp/hub)*rpar(5)*y(4)
	djac(1,5)= .5e0_dp*(1e0_dp/hub)*rpar(5)*y(5)
	djac(2,1)= 0e0_dp
	djac(2,2)= 0e0_dp
	djac(2,3)= 0e0_dp
	djac(2,4)= 1e0_dp
	djac(2,5)= 0e0_dp
	djac(3,1)= 0e0_dp
	djac(3,2)= 0e0_dp
	djac(3,3)= 0e0_dp
	djac(3,4)= 0e0_dp
	djac(3,5)= 1e0_dp
	djac(4,1)= 0e0_dp
	djac(4,2)= -1.5e0_dp*(1e0_dp/hub)*rpar(5)*y(4)*d_phi_v - dd_phi_v
	djac(4,3)= -1.5e0_dp*(1e0_dp/hub)*rpar(5)*y(4)*d_psi_v - d_phi_d_psi_v
	djac(4,4)= -1.5e0_dp*(1e0_dp/hub)*rpar(5)*y(4)*y(4) - 3e0_dp*hub
	djac(4,5)= -1.5e0_dp*(1e0_dp/hub)*rpar(5)*y(4)*y(5)
	djac(5,1)= 0e0_dp
	djac(5,2)= -1.5e0_dp*(1e0_dp/hub)*rpar(5)*y(5)*d_phi_v - d_phi_d_psi_v
	djac(5,3)= -1.5e0_dp*(1e0_dp/hub)*rpar(5)*y(5)*d_psi_v - dd_psi_v
	djac(5,4)= -1.5e0_dp*(1e0_dp/hub)*rpar(5)*y(4)*y(5)
	djac(5,5)= -1.5e0_dp*(1e0_dp/hub)*rpar(5)*y(5)*y(5) - 3e0_dp*hub

	!Success
	IER = 0

end subroutine FCVDJAC








