!*******************************************************************************
!Layne Price, University of Auckland, May 15, 2012.
!*******************************************************************************

!SUMMARY:
!A program that does integration for two field hybrid inflation.  Uses the
!integrator FCVODE from the LLNL SUNDIALS package, which requires the RHS of the
!ODE to be expressed in the external subroutine FCVFUN and the Jacobian in the
!external subroutine FCVDJAC.  These are included below the main program.
!SUNDIALS functions are stored in a library.

!OUTLINE:
!The program architecture is as follows: We parallelize using OpenMPI; load ICs
!according to the method specified in the namelist; integrate the initial
!condition; sort the initial condition into a "success" array if it reaches N>65
!and into a "fail" array if it reaches the minimum of the potential without
!inflating.  A new IC is chosen and the integration is repeated until we get
!enough points.  Counters on the number of points found are collected on the
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
		& errorlocal, ierr, rc, burnin
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
	call mpi_init(ierr)
		if(ierr .ne. mpi_success) then
			print*,"Error parallelizing."
			call mpi_abort(mpi_comm_world, rc, ierr)
			stop
		end if
	!Obtains info on processors.
	call mpi_comm_rank(mpi_comm_world, rank, ierr)
	call mpi_comm_size(mpi_comm_world, numtasks, ierr)
	if(printing) print*,'Number of tasks=',numtasks,' My rank=',rank
	call mpi_barrier(mpi_comm_world,ierr)

	!Opens success and fail files. Optionally traj files.  Defaults to unformatted binary.
	call open_hybridfiles(rank,numtasks,sucunit,failunit)
	if (traj) call open_trajfiles(rank, trajnumb)

	!Set potential's parameters.
	call parameters_hybrid()

	!Print stats.
	if(rank==0) call hybrid_initstats(ic, printing)

	!Set seed for each thread.
	call init_random_seed(rank)

	!Set some params for FCVODE integrator.
	call set_paramsFCVODE(rpar, neq, nglobal, numtasks, iatol, atol, rtol, &
		& meth, itmeth, t0, t, itask, tout, dt)

	!Set first IC.
  burnin=10000
  call ic_init(ic, y0, iccounter, sample_table, burnin, rank,&
  &numtasks, ic_table, yref, toler, printing)
	Y=Y0

	!Initialize FCVODE integrator.
	call fcvmalloc(t0, y0, meth, itmeth, iatol, rtol, atol,&
		&iout, rout, ipar, rpar, ier)
  if (.not. traj) then
	  call fcvsetiin("MAX_NSTEPS", 5000000, ier)
  else
    !If recording trajectories, we want the integrator to return at every step
    !taken so that we can use the adaptive step size control in FCVODE to get
    !important data. Set with itask=2.
    itask=2
  end if
	call fcvdense(neq, ier)
	call fcvdensesetjac (1, ier)

	!Loop over ICs until achieve numb of desired points.
	integr_ch=.true.	!Exit condition.
  icloop: 	do while (integr_ch)

		!Count numb of times each thread gets a new IC.
		localcount = localcount + 1
		!Get new point if on second or greater run.
		if (localcount>1) then
			call new_point(y0,iccounter,sample_table, ic, &
				&ic_table, yref, toler)

			!Reinit time and Y
			y=y0
			t0=0e0_dp
      tout =t0+ dt
			t=t0
			!Reinitialize integrator.
      if (.not. traj) then
        itask = 1
      else
        itask=2
      end if
			call fcvreinit(t0, y0, iatol, rtol, atol, ier)
		end if
		!Counters
		success = 0
		iccounter = iccounter + 1

    !###################################################################
		!Perform the integration.
    iend=3000000
    intloop:	do i=1,iend

			!Take field values if recording trajectory.
			if (traj .and. mod(i,10)==0) call rec_traj(Y, ytraj)

			!*********************************
			!Perform the integration. dt set in namelist ics.
			call fcvode(tout,t,y,itask,ier)
print*,t
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
			if(printing .and. mod(i,iend/10)==0 .and. .not. traj) print*,"i is getting big...",i
		end do intloop
    !###################################################################


		!Print the traj & delete -- O(2n).
		if (success==1 .and. traj) then
      call print_del_traj(ytraj, trajnumb)
    !If not successful, then just delete list.
    else if (success==0 .and. traj) then
      call ll_del_all(ytraj)
    end if

		!Check if not finding any succ points.
		call all_fail_check(successlocal, faillocal, allfailcheck, printing, check)
		if (allfailcheck) exit icloop

		!Determine loop exit condition.
		if (ic<5 .or. ic==6) then
			integr_ch=(successlocal<points)
		else if (ic==5) then
			integr_ch=(localcount<size(ic_table,1))
    else if (ic==7) then
      exit
		end if

	end do icloop

	!Halts processors here.
	call mpi_barrier(mpi_comm_world,ierr)
	!Gives child data to master.
	call mpi_reduce(badfieldlocal,badfieldcounter,1,mpi_integer,&
		&mpi_sum,0,mpi_comm_world,ierr)
	call mpi_reduce(successlocal,counter,1,mpi_integer,&
		&mpi_sum,0,mpi_comm_world,ierr)
	call mpi_reduce(faillocal,failcount,1,mpi_integer,&
		&mpi_sum,0,mpi_comm_world,ierr)
	call mpi_reduce(errorlocal,errorcount,1,mpi_integer,mpi_sum,&
		&0,mpi_comm_world,ierr)
	call mpi_barrier(mpi_comm_world,ierr)

	!Print from master.
	if(rank==0) then
		call hybrid_finalstats(ic, counter, failcount, &
		&badfieldcounter, errorcount, printing)
	end if

	!Clean up integrator.
	call fcvfree

	!End parallel.
	call mpi_finalize(ierr)

	!Why this is necessary I have no idea.  It works fine without it on my computer,
  !but MPI gives an error if this isn't here on the cluster.
	stop

end program hybrid_integrator_d




!***********************************************************************
!RHS of equation to integrate.
!***********************************************************************
subroutine fcvfun(t, y, ydot, ipar, rpar, ier)
  use types, only : dp
  implicit none
	
	!******************
	real(dp), intent(in) :: y(*), t
	real(dp), intent(out) :: ydot(*)
	integer, intent(in) :: ipar(*)
	integer, intent(out) :: ier
	real(dp), intent(in) :: rpar(*)
	!******************

	real(dp) :: v, d_phi, d_psi, hub

	!Potential and its derivatives.
	v = (rpar(1)*rpar(1)*rpar(1)*rpar(1))*((1e0_dp - &
		&((y(3)*y(3))/(rpar(2)*rpar(2))))**2 &
		&+ ((y(2)*y(2))/(rpar(3)*rpar(3)))&
		& + ((y(2)*y(2)*y(3)*y(3))/ (rpar(4)**4)))

	d_phi = 2e0_dp*(rpar(1)*rpar(1)*rpar(1)*rpar(1))*(((y(2))/(rpar(3)*rpar(3)))+ &
		&((y(2)*y(3)*y(3))/(rpar(4)**4)))

	d_psi = 2e0_dp*(rpar(1)*rpar(1)*rpar(1)*rpar(1))*(((-2e0_dp*y(3))/&
		&(rpar(2)*rpar(2)))*(1e0_dp &
		&-((y(3)*y(3))/(rpar(2)*rpar(2))))+&
		& ((y(2)*y(2)*y(3))/(rpar(4)**4)))
	
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
	ydot(4) = -3e0_dp*hub*y(4)-d_phi
	ydot(5) = -3e0_dp*hub*y(5)-d_psi

	!Success
	ier = 0

end subroutine fcvfun

!***************************************************************************
!Jacobian of RHS of equation to integrate.
subroutine fcvdjac (neq, t, y, fy, djac, h, ipar, rpar,&
			&wk1, wk2, wk3, ier)
  use types, only : dp
  implicit none

	!**************************
	real(dp), intent(in) :: y(*), fy(*), t, h
	integer, intent(in) :: ipar(*), neq
	integer, intent(out) :: ier
	real(dp), intent(inout) :: djac(neq,*)
	real(dp), intent(in) :: rpar(*)
	real(dp), intent(in) :: wk1(*), wk2(*), wk3(*)
	!**************************

	real(dp) :: v, d_phi, d_psi, dd_phi, dd_psi, d_phi_d_psi, Hub

	!Potential, V, and its derivatives.
	v = (rpar(1)**4)*((1e0_dp - ((y(3)*y(3))/(rpar(2)*rpar(2))))**2 +&
		& ((y(2)*y(2))/(rpar(3)*rpar(3)))&
		& + ((y(2)*y(2)*y(3)*y(3))/ (rpar(4)**4)))

	d_phi = 2e0_dp*(rpar(1)**4)*(((y(2))/(rpar(3)*rpar(3)))+&
		& ((y(2)*y(3)*y(3))/(rpar(4)**4)))

	d_psi = 2e0_dp*(rpar(1)**4)*(((-2e0_dp*y(3))/(rpar(2)*rpar(2)))*(1e0_dp &
		&-((y(3)*y(3))/(rpar(2)*rpar(2))))+ &
		&((y(2)*y(2)*y(3))/(rpar(4)**4)))

	dd_phi = 2e0_dp*(rpar(1)**4)*((1e0_dp/(rpar(3)*rpar(3)))+&
		&((y(3)*y(3))/(rpar(4)**4)))

	dd_psi = 2e0_dp*(rpar(1)**4)*(((-2e0_dp)/(rpar(2)*rpar(2)))+ &
		&((4e0_dp*y(3)*y(3))/(rpar(2)**4))+&
		& ((y(2)*y(2))/(rpar(4)**4)))

	d_phi_d_psi = (4e0_dp*(rpar(1)**4)*y(2)*y(3))/(rpar(4)**4)

	hub = sqrt(rpar(5)*(.5e0_dp*((y(4)*y(4)) + (y(5)*y(5))) + v))
  if (hub<0) then
		print*,"Hubble parameter < 0"
		ier=1
		return
	end if

	!Partial derivatives of the RHS vector in system of equations.
	djac(1,1)= 0e0_dp
	djac(1,2)= .5e0_dp*(1e0_dp/hub)*rpar(5)*d_phi
	djac(1,3)= .5e0_dp*(1e0_dp/hub)*rpar(5)*d_psi
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
	djac(4,2)= -1.5e0_dp*(1e0_dp/hub)*rpar(5)*y(4)*d_phi - dd_phi
	djac(4,3)= -1.5e0_dp*(1e0_dp/hub)*rpar(5)*y(4)*d_psi - d_phi_d_psi
	djac(4,4)= -1.5e0_dp*(1e0_dp/hub)*rpar(5)*y(4)*y(4) - 3e0_dp*hub
	djac(4,5)= -1.5e0_dp*(1e0_dp/hub)*rpar(5)*y(4)*y(5)
	djac(5,1)= 0e0_dp
	djac(5,2)= -1.5e0_dp*(1e0_dp/hub)*rpar(5)*y(5)*d_phi - d_phi_d_psi
	djac(5,3)= -1.5e0_dp*(1e0_dp/hub)*rpar(5)*y(5)*d_psi - dd_psi
	djac(5,4)= -1.5e0_dp*(1e0_dp/hub)*rpar(5)*y(4)*y(5)
	djac(5,5)= -1.5e0_dp*(1e0_dp/hub)*rpar(5)*y(5)*y(5) - 3e0_dp*hub

	!Success
	ier = 0

end subroutine fcvdjac
