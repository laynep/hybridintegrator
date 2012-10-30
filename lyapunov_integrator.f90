!*******************************************************************************
!Layne Price, University of Auckland, October 8, 2012.
!*******************************************************************************

!SUMMARY:
!A program that does integration for two field hybrid inflation.  Uses the
!integrator LESNLS from http://people.math.gatech.edu/~dieci/software-les.html
!by Luca Dieci to perform nonlinear integration and to calculate the Lyapunov
!exponents.

!CAVEAT:
!A good calculation of the Lyapunov exponents requires integration over a large
!range of time.  However, because the system that we are integrating is
!dissipative, we know that at (very) late times all orbits will be at one of the
!fixed points corresponding to the minima of the potential.  Furthermore,
!friction will inevitably force the late time Lyapunov exponents --> 0.
!Consequently, any calculation of the Lyapunov exponents is merely suggestive
!and is truly representative only for systems with very little to no friction
!and in the late time limit.  Furthermore, the Lyapunov exponents are not
!invariant with respect to spacetime diffeomorphism, so the value we calculate
!here will be meaningful only for comoving observers.  However, the overall sign
!of the LE will remain the same, implying that if we find a LE>0, then any
!observer should be able to find a positive LE.

!OUTLINE:
!The program architecture is as follows: We parallelize using OpenMPI; load ICs
!according to the method specified in the namelist; integrate the initial
!condition and calculate the Lyapunov exponent; sort the initial condition into
!a "success" array if it reaches N>65 and into a "fail" array if it reaches the
!minimum of the potential without inflating. The Lyapunov exponents can be
!either recorded in a linked list at each integration step or ????.
!A new IC is then chosen and the integration is repeated until we get
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
!module.  All routines needed for the calculation of Lyapunov exponents are
!contained in the lyapunov_subroutines module.

!NOTE:  The Y vector corresponds to:
!Y(1)=N, Y(2)=phi, Y(3)=psi, Y(4)=phi_dot, Y(5)=psi_dot.
!*******************************************************************************


program lyapunov_integrator
  use d_hybrid_initialconditions
  use hybrid_subroutines
  use lyapunov_subroutines
  use rng, only : init_random_seed
  use mpi
  use linked_list
  use types, only : dp
  use features, only : newunit
  implicit none

	!Main variables.
	real(dp), dimension(5) :: y, y0		!The fields.
	real(dp), dimension(5) :: yref		!Ref field for zooming proc.
	real(dp) :: t0, tout, dt		!Timers.
	!Counters & unit numbers
	integer:: i, j, success, counter, iccounter, sucunit,&
		& failunit, failcount, localcount, burnin
	integer :: errorcount, iend
	integer :: badfieldcounter, badfieldlocal, successlocal, faillocal,&
		& errorlocal, ierr, rc
	!Program variables.
	integer :: ic, trajnumb, points, u, numb
	real(dp) :: check, v, ratio, toler
	logical :: leave, allfailcheck, printing, traj, lyap
	logical :: integr_ch
  character(len=12) :: fname
	!Variables to load IC from file for direct read or interpolation (ic=4,5)
	real(dp), dimension(:,:), allocatable :: sample_table, ic_table
	!Lists to record trajectory and lyapunov exponents.
	type(linkedlist) :: ytraj, le_list
	!Parallel variables.
	integer :: numtasks, rank
  !LESNLS parameters
  external :: getdf, getf
  integer, parameter :: le_m=5, n=3 !Dimn of problem=le_m, and numb of LEs=n.
  integer, parameter :: ifdim=le_m*le_m+11*le_m*n+13*le_m+8*n+63
  real(dp) :: tolt,tolq,toll(n)
  real(dp) :: lyap_exp(n), dt_lyap
  integer :: ipar(13), iflag, inarr(5)
  real(dp) :: fwork(ifdim), x0(le_m,n), rearr(5)

	namelist /ics/ points, IC, dt, printing, traj

	!***********************************************************************
	!***********************************************************************

	!Read numb of (data points per numb of processes) & IC type from file.
	!Do we want to print to stdout?  Do we want to record the trajectory?
	open(unit=newunit(u), file="parameters_hybrid.txt", status="old", delim = "apostrophe")
	read(unit=u, nml=ics)
	close(unit=u)

  !Record Lyapunov Exponents at each step?
  lyap=.true.
  !Test with or without friction?
  friction=.true.

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

	!Set some params for LESNLS integrator.
  call lesnls_parameters(ipar,rearr,tolt,tolq,toll,ifdim,x0,lyap)
  call reinit_lesnls(t0,tout,dt,ipar)

	!Set first IC.
  burnin=10000
  call ic_init(ic, y0, iccounter, sample_table, burnin, rank,&
  &numtasks, ic_table, yref, toler, printing)
	Y=Y0

  !Initialize integrator.
  call init(le_m,n,ipar,t0,tout,fwork,iflag)
    if (iflag .ne. 0) then 
      print*,"ERROR in INIT. IFLAG = ", iflag
      stop
    end if

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
      call reinit_lesnls(t0,tout,dt,ipar)
		end if
		!Counters
		success = 0
		iccounter = iccounter + 1

    !###################################################################
		!Perform the integration.
    if (friction) then
      iend=30000000
    else
      iend=20000
    end if

    intloop:	do i=1,iend

			!Take field values if recording trajectory. Previously initialized
			if (traj .and. mod(i,10)==0) call rec_traj(Y, ytraj)

			!*********************************
			!Perform the integration. dt set in namelist ics.
      call lesnls(getf,getdf,le_m,n,lyap_exp,t0,tout,dt_lyap,y,x0,&
        &tolt,tolq,toll,ipar,fwork,iflag,inarr,rearr)
        if (iflag .ne. 0 ) then
          print*, "ERROR in lesnls. IFLAG is ", iflag
          stop
        end if
      tout = tout + dt
      !*********************************

      !Record Lyapunov exponents and time as linked list.
      if (lyap) then
        if (i<600) then
          call rec_LE(lyap_exp,t0,y,le_list)
        else if (mod(i,50)==0) then
          call rec_LE(lyap_exp,t0,y,le_list)
        end if
      end if

      !Check succ or fail condition: N>65 success=1, and N<65 failure=0.
			call succ_or_fail(Y, success, successlocal, faillocal,&
				& sucunit, failunit, ic, leave, printing)
			!If condition met, leave=.true.
			if (leave) exit intloop
			
			!Tells if not enough time for fields to evolve.
			!if (i==iend) then
			!	errorlocal = errorlocal + 1
			!	if (printing) print*, "Error: field didn't reach minima."
			!end if
			!if(printing .and. MOD(i,iend/10)==0 .and. .not. traj) print*,"i is getting big...",i
		end do intloop
    !###################################################################


		!Print the traj & delete -- O(2n).
		if (success==1 .and. traj) then
      call print_del_traj(ytraj, trajnumb)
    !If not successful, then just delete list.
    else if (success==0 .and. traj) then
      call ll_del_all(ytraj)
    end if

    !Open file and write LEs.
    numb=newunit(numb)
    write(fname,'(a,i4.4,a)') "lexp",numb,".bin"
    open(unit=numb,file=fname,status="new",form="unformatted")
    if (lyap) then
      !Print the Lyap exps.
      call print_del_traj(le_list,numb)
    else
      write(unit=numb) Y(1),psi_0,phi_0,psi_dot_0,phi_dot_0,&
        &(lyap_exp(j),j=1,size(lyap_exp))
    end if

		!Check if the integrator isn't finding any succ points.
		call all_fail_check(successlocal, faillocal, allfailcheck, printing, check)
		if (allfailcheck) exit icloop
		!Determine loop exit condition.
    if (.not. friction) then
      integr_ch=(localcount<points)
    else if (ic<5 .or. ic==6) then
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
    if (printing) call le_stats(ipar, iflag, t0)
    if (printing) print*, "Columns in LE files: ", size(lyap_exp)+6
	end if

	!End parallel.
	call MPI_FINALIZE(ierr)

	!Why this is necessary I have no idea.  It works fine without it on my computer,
  !but MPI gives an error if this isn't here on the cluster.
	stop

end program lyapunov_integrator


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !subroutine to specify the rhs of linear equation of y'=a*y.

  subroutine getdf(le_m,y,a,inarr,rearr)
    use lyapunov_subroutines, only : friction
  	use d_hybrid_initialconditions
  	implicit none

  	integer, intent(inout) :: le_m
  	real(dp), dimension(le_m), intent(inout) :: y
  	real(dp), dimension(le_m,le_m), intent(inout) :: a
  	integer, dimension(5), intent(inout) :: inarr
  	real(dp), dimension(5), intent(inout) :: rearr
  	real(dp) :: v, d_phi, d_psi, dd_phi, dd_psi, d_phi_d_psi, hub

  	!potential, v, and its derivatives.
  	v = (rearr(1)**4e0_dp)*((1e0_dp - ((y(3)*y(3))/(rearr(2)*rearr(2))))**2e0_dp &
  	&+ ((y(2)*y(2))/(rearr(3)*rearr(3))) + ((y(2)*y(2)*y(3)*y(3))/ (rearr(4)**4e0_dp)))

  	d_phi = 2e0_dp*(rearr(1)**4e0_dp)*(((y(2))/(rearr(3)*rearr(3)))+ &
  	&((y(2)*y(3)*y(3))/(rearr(4)**4e0_dp)))

  	d_psi = 2e0_dp*(rearr(1)**4e0_dp)*(((-2e0_dp*y(3))/(rearr(2)*rearr(2)))*(1e0_dp&
  	& -((y(3)*y(3))/(rearr(2)*rearr(2))))+ ((y(2)*y(2)*y(3))/(rearr(4)**4e0_dp)))

  	dd_phi = 2e0_dp*(rearr(1)**4e0_dp)*((1e0_dp/(rearr(3)*rearr(3)))+&
      &((y(3)*y(3))/(rearr(4)**4e0_dp)))

  	dd_psi = 2e0_dp*(rearr(1)**4e0_dp)*(((-2e0_dp)/(rearr(2)*rearr(2)))+ &
  	&((4e0_dp*y(3)*y(3))/(rearr(2)**4e0_dp))+ ((y(2)*y(2))/(rearr(4)**4e0_dp)))

  	d_phi_d_psi = (4e0_dp*(rearr(1)**4e0_dp)*y(2)*y(3))/(rearr(4)**4e0_dp)

  	hub = sqrt(rearr(5)*(.5e0_dp*((y(4)*y(4)) + (y(5)*y(5))) + v))
	
  	!partial derivatives of the rhs vector in system of equations.
    if (friction) then
  	  a(1,1)= hub
  	  a(1,2)= .5e0_dp*(1e0_dp/hub)*rearr(5)*d_phi
  	  a(1,3)= .5e0_dp*(1e0_dp/hub)*rearr(5)*d_psi
  	  a(1,4)= .5e0_dp*(1e0_dp/hub)*rearr(5)*y(4)
  	  a(1,5)= .5e0_dp*(1e0_dp/hub)*rearr(5)*y(5)
  	  a(2,1)= 0e0_dp
  	  a(2,2)= 0e0_dp
  	  a(2,3)= 0e0_dp
  	  a(2,4)= 1e0_dp
  	  a(2,5)= 0e0_dp
  	  a(3,1)= 0e0_dp
  	  a(3,2)= 0e0_dp	
  	  a(3,3)= 0e0_dp
  	  a(3,4)= 0e0_dp
  	  a(3,5)= 1e0_dp
  	  a(4,1)= 0e0_dp
  	  a(4,2)= -1.5e0_dp*(1e0_dp/hub)*rearr(5)*y(4)*d_phi - dd_phi
  	  a(4,3)= -1.5e0_dp*(1e0_dp/hub)*rearr(5)*y(4)*d_psi - d_phi_d_psi
  	  a(4,4)= -1.5e0_dp*(1e0_dp/hub)*rearr(5)*y(4)*y(4) - 3e0_dp*hub
  	  a(4,5)= -1.5e0_dp*(1e0_dp/hub)*rearr(5)*y(4)*y(5)
  	  a(5,1)= 0e0_dp
  	  a(5,2)= -1.5e0_dp*(1e0_dp/hub)*rearr(5)*y(5)*d_phi - d_phi_d_psi
  	  a(5,3)= -1.5e0_dp*(1e0_dp/hub)*rearr(5)*y(5)*d_psi - dd_psi
  	  a(5,4)= -1.5e0_dp*(1e0_dp/hub)*rearr(5)*y(4)*y(5)
  	  a(5,5)= -1.5e0_dp*(1e0_dp/hub)*rearr(5)*y(5)*y(5) - 3e0_dp*hub
    else
      a=0e0_dp
      a(4,2)=-1e0_dp*dd_phi
      a(4,3)=-1e0_dp*d_phi_d_psi
      a(5,2)=-1e0_dp*d_phi_d_psi
      a(5,3)=-1e0_dp*dd_psi
      a(2,4)=1e0_dp
      a(3,5)=1e0_dp
    end if


  end subroutine getdf

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !subroutine to specify the rhs of y'=f[y].
  !parameters are defined in d_hybrid_initialconditions

  subroutine getf(le_m,y,ydot,inarr,rearr)
    use lyapunov_subroutines, only : friction
  	use d_hybrid_initialconditions
  	implicit none

  	integer, intent(in) :: le_m
  	real(dp), dimension(5), intent(inout) :: y
    real(dp), dimension(5), intent(inout) :: ydot
  	integer, dimension(5), intent(in) :: inarr
  	real(dp), dimension(5), intent(in) :: rearr
  	real(dp) :: v, d_phi, d_psi, dd_phi, dd_psi, d_phi_d_psi, hub

  	!Potential, v, and its derivatives.
  	v = (rearr(1)**4e0_dp)*((1e0_dp - ((y(3)*y(3))/(rearr(2)*rearr(2))))**2e0_dp + &
  		&((y(2)*y(2))/(rearr(3)*rearr(3))) + ((y(2)*y(2)*y(3)*y(3))/ (rearr(4)**4e0_dp)))

  	d_phi = 2e0_dp*(rearr(1)**4e0_dp)*(((y(2))/(rearr(3)*rearr(3)))+ &
      &((y(2)*y(3)*y(3))/(rearr(4)**4e0_dp)))

  	d_psi = 2e0_dp*(rearr(1)**4e0_dp)*(((-2e0_dp*y(3))/(rearr(2)*rearr(2)))*(1e0_dp -&
      &((y(3)*y(3))/(rearr(2)*rearr(2))))&
  		&+ ((y(2)*y(2)*y(3))/(rearr(4)**4e0_dp)))

  	hub = sqrt(rearr(5)*(.5e0_dp*((y(4)*y(4)) + (y(5)*y(5))) + v))

   	!Equations of motion with and without friction.
    if (.not. friction) hub = 0e0_dp
  	ydot(1) = hub
  	ydot(2) = y(4)
  	ydot(3) = y(5)
  	ydot(4) = -3e0_dp*hub*y(4)-d_phi
  	ydot(5) = -3e0_dp*hub*y(5)-d_psi	

  end subroutine getf
