!*******************************************************************************
!Layne Price, University of Auckland, May 15, 2012.
!*******************************************************************************

!Module that contains the subroutines necessary to run the main program
!hybrid_integrator_d.  However, this module does not contain the procedures
!necessary to obtain the ICs.  These are further specified in
!d_hybrid_initialconditions.f90.

module hybrid_subroutines
  use types, only : dp
  use d_hybrid_initialconditions
  implicit none

contains

!Open the files to write the successful and failure points to for each thread.
subroutine open_hybridfiles(rank,numtasks,sucunit,failunit, datatype, formt)
	implicit none

	integer, intent(in) :: rank, numtasks
	integer, intent(out) :: sucunit, failunit
	character(len=*), optional, intent(inout) :: datatype, formt
	character(len=12) :: failname, sucname

	!Get the units to write to.
	sucunit = rank + 20
	failunit = rank + 20 + numtasks

	!Create the file name.
	if (present(datatype)) then
		write(sucname,'(a,i4.4,a)')'succ',sucunit,datatype
		write(failname,'(a,i4.4,a)')'fail',failunit,datatype
	else
		write(sucname,'(a,i4.4,a)')'succ',sucunit,'.bin'
		write(failname,'(a,i4.4,a)')'fail',failunit,'.bin'
	end if

	!Open.
	if (present(formt)) then
		open(unit=sucunit,status='new',file=sucname, form=formt)
		open(unit=failunit,status='new',file=failname, form=formt)
	else
		open(unit=sucunit,status='new',file=sucname, form='unformatted')
		open(unit=failunit,status='new',file=failname, form='unformatted')
	end if

end subroutine open_hybridfiles

!Opens files for storing the trajectories.
subroutine open_trajfiles(rank, trajnumb, datatype, formt)
	implicit none

	integer, intent(in) :: rank
	integer, intent(out) :: trajnumb
	character(len=12) :: trajname
	character(len=*), optional, intent(inout) :: datatype, formt
	
	!Write to file.
	trajnumb=3141+rank

	!Create the file name.
	if (present(datatype)) then
		write(trajname,'(a,i4.4,a)')'traj',trajnumb,datatype
	else
		write(trajname,'(a,i4.4,a)')'traj',trajnumb,'.bin'
	end if

	!Open.
	if (present(formt)) then
		open(unit=trajnumb,status='new',file=trajname, form=formt)
	else
		open(unit=trajnumb,status='new',file=trajname, form='unformatted')
	end if

end subroutine open_trajfiles

!Print the initial stats to file info.200.
subroutine hybrid_initstats(ic,printing,infounit)
	implicit none
	
	integer, intent(in) :: ic
	integer, optional, intent(in) :: infounit
	integer :: u
	logical, intent(in) :: printing

	if(present(infounit)) then
		u=infounit
	else
		u=200
	end if

	open(unit=u,status='new',file='info.200')
	if(printing) print*, "Hybrid Inflation -- Homogeneous"
	write(unit=u,fmt=*) "Hybrid Inflation -- Homogeneous"
	if (ic == 1) then
		if(printing) print*,"ZERO VELOCITY SLICE"
		write(unit=u,fmt=*) "ZERO VELOCITY SLICE"
	else if (ic == 2) then
		if(printing) print*,"EQUAL ENERGY SLICE"
		write(unit=u,fmt=*) "EQUAL ENERGY SLICE"
	else if (ic == 3) then
		if(printing) print*,"2D SLICING OF EQEN SLICE"
        write(unit=u,fmt=*) "2D SLICING OF EQEN SLICE"
	else if (ic == 4) then
		if(printing) print*,"IC FROM METROPOLIS SAMPLING"
        write(unit=u,fmt=*) "IC FROM METROPOLIS SAMPLING"
	else if (ic==5) then
		if(printing) print*,"IC FROM FILE"
        write(unit=u,fmt=*) "IC FROM FILE"
    else if (ic==6) then
        if (printing) print*, "IC FROM ZOOMING IN ON POINT"
        write(unit=u,fmt=*) "IC FROM ZOOMING IN ON POINT"
	end if

	if(printing) print*,"Planck mass is ",m_planck
	write(unit=u,fmt=*) "Planck mass is ",m_planck
	if(printing) print*,"E=",energy_scale
	write(unit=u,fmt=*) "E=",energy_scale
	if(printing) print*,"Lambda is ",lambda
	write(unit=u,fmt=*) "Lambda is ",lambda
       	if(printing) print*,"M is ",m
	write(unit=u,fmt=*) "M is ",m
       	if(printing) print*,"mu is ",mu
	write(unit=u,fmt=*) "mu is ",mu
       	if(printing) print*,"nu is ",nu
	write(unit=u,fmt=*) "nu is ",nu

end subroutine hybrid_initstats

!Print the final stats to info.200.
subroutine hybrid_finalstats(ic, counter, failcount, badfieldcounter, &
		& errorcount, printing,  infounit)
	implicit none

	integer, intent(inout) :: ic, counter, failcount, badfieldcounter, errorcount
	integer, optional, intent(in) :: infounit
	integer :: u
	real(dp) :: ratio
	logical, intent(in) :: printing

	if(present(infounit)) then
		u=infounit
	else
		u=200
	end if

	ratio=DBLE(counter)/DBLE(counter+failcount)
	if (ic == 1) then
		if(printing) print*,"ZERO VELOCITY SLICE"
		write(unit=u,fmt=*) "ZERO VELOCITY SLICE"
	else if (ic == 2) then
		if(printing) print*,"EQUAL ENERGY SLICE"
		write(unit=u,fmt=*) "EQUAL ENERGY SLICE"
	else if (ic==3) then
		if(printing) print*, "2D SLICING OF EQEN SLICE."
		write(unit=u,fmt=*) "2D SLICING OF EQEN SLICE"
	else if (ic == 4) then
		if(printing) print*,"IC FROM METROPOLIS SAMPLING"
    write(unit=u,fmt=*) "IC FROM METROPOLIS SAMPLING"
	else if (ic==5) then
		if(printing) print*,"IC FROM FILE"
    write(unit=u,fmt=*) "IC FROM FILE"
	else if (ic==6) then
    if (printing) print*, "IC FROM ZOOMING IN ON POINT"
    write(unit=u,fmt=*) "IC FROM ZOOMING IN ON POINT"
	end if

	if(printing) print*,"Number of succ points ",counter,&
		&"Number of fail points ",failcount," ratio ",ratio
	write(unit=u,fmt=*) "Number of succ points ",counter,&
		&"Number of fail points ",failcount," ratio ",ratio
	if(printing) print*,"Number of IC that gave bad field values ",badfieldcounter
	write(unit=u,fmt=*)"Number of IC that gave bad field values",&
		&badfieldcounter
	if(printing) print*,"Number that didn't fall into minimum ",errorcount
	write(unit=u,fmt=*) "Number that didn't fall into minimum ",errorcount


end subroutine hybrid_finalstats

!Get a new point with method specified by "IC."
subroutine new_point(y0,iccounter,sample_table,ic, ic_table, yref, toler)
	implicit none

	real(dp), dimension(:), intent(inout) :: y0
	real(dp), dimension(:), optional, intent(inout) :: yref
	real(dp), optional, intent(in) :: toler
	integer, intent(inout) :: iccounter
	integer, intent(in) :: ic
	real(dp), dimension(:,:), allocatable, intent(inout) :: sample_table
	real(dp), dimension(:,:), optional, intent(in) :: ic_table


	if (ic == 1) then
		call d_ic_zerov(y0)
	else if (ic == 2) then
		call d_ic_eqen(y0,iccounter)
	else if (ic==3) then
		call eqen_slicing(y0)
	else if (ic==4) then
		call ic_metr(y0,sample_table,iccounter)
	else if (ic==5) then
		call ic_fromarray(y0,ic_table,iccounter+1)
	else if (ic==6) then
		call ic_eqen_pert(yref,y0,iccounter,euclidean,toler)
  else if (ic==7) then
    stop
  end if

end subroutine new_point

!Determine if an iteration of the integrator has reached a success or failure condition yet.  If it has, then this routine returns "leave=.true." -- an exit condition for the loop in the main program.
subroutine succ_or_fail(Y, success, successlocal, faillocal, &
		&sucunit, failunit, ic, leave, printing)
	implicit none

	integer, intent(inout) :: success, successlocal, faillocal, sucunit, failunit, ic
	real(dp), dimension(:), intent(inout) :: Y
	logical, intent(inout) :: leave
	real(dp) :: V, check
	logical, intent(in) :: printing

	!Check if  fields fell into minima (if engy dnsty <
	! dnsty at infl end)
	check = (.5_dp*( (Y(4)*Y(4)) + (Y(5)*Y(5))) &
		&+V_h(Y) - (lambda*lambda*lambda*lambda))

	!Gives a 1 if inflation is successful.
	leave = .false.
	if (Y(1)>65_dp) then
		success = 1
		successlocal = successlocal + 1
		if (ic==1) then
			write(unit=sucunit),psi_0,phi_0
		else
			write(unit=sucunit), psi_0, phi_0,&
			& psi_dot_0, phi_dot_0
		end if
		if(printing .and. mod(successlocal,1000)==0) then
			print*,successlocal,success
		end if
		leave = .true.
	elseif (check<0_dp) then
		success = 0
		faillocal = faillocal + 1
		if (ic==1) then
			write(unit=failunit),psi_0,phi_0
		else
			write(unit=failunit), psi_0, phi_0,&
			& psi_dot_0, phi_dot_0
		end if
		if(printing .and. mod(successlocal,1000)==0) then
			print*,successlocal,success
		end if
		leave = .true.
	end if

end subroutine succ_or_fail

!Set some necessary parameters for the integrator FCVODE.
subroutine set_paramsFCVODE(rpar, neq, nglobal, numtasks, iatol, atol, rtol, &
		& meth, itmeth, t0, t, itask, tout, dt)
	implicit none

	!FCVODE params
	real(dp), intent(out) :: rpar(5)
	integer, intent(out) :: meth, itmeth
	integer(kind=8), intent(out) :: neq, nglobal
	integer, intent(out) :: iatol, itask
	real(dp), intent(out) :: t0, t, tout, rtol, atol(5)
	real(dp), intent(in) :: dt

	!Other params
	integer, intent(in) :: numtasks
	integer :: ier

	!Pass params to FCVODE subroutines.  Errors if don't do this...?
	rpar(1) = lambda
	rpar(2) = m
	rpar(3) = mu
	rpar(4) = nu
	rpar(5) = beta

	!DOF global and local.
	neq = 5
	nglobal = neq*numtasks

	!Set params for FCVODE.
	call fnvinits(1, neq, ier)
		if (ier .ne. 0) then
			print*, ' SUNDIALS_ERROR: FNVINITS returned IER = ', ier
			stop
		endif
	iatol = 2
	atol = 1e-10_dp
	rtol = 1e-10_dp
	meth = 2
	itmetH = 2
	t0=0_dp
	t=t0
	itask = 1
	tout = dt

end subroutine set_paramsFCVODE

!Subroutine to check whether the integrator is giving any success points.  If check is present, then this will also check to see if the ratio of succ to fail points is less than check.
subroutine all_fail_check(successlocal, faillocal, allfailcheck, printing, check)
	implicit none

	integer, intent(in) :: successlocal, faillocal
	logical, intent(out) :: allfailcheck
	real(dp), optional, intent(in) :: check
	real(dp) :: ratio
	logical, intent(in) :: printing

	allfailcheck=.false.

	if (present(check)) then
		ratio = dble(successlocal)/dble(successlocal+faillocal)
		if(successlocal>=1.AND. ratio<check) then
			if(printing) print*,"Few successful.  Ratio is ", ratio
			allfailcheck=.true.
		end if
	else
		if(successlocal==0 .AND. faillocal>10000) then
			if(printing) print*,"None successful and fail is 10000."
			allfailcheck=.true.
		end if
	end if

end subroutine all_fail_check


!Subroutine to record the trajectory in a linked list..
subroutine rec_traj(Y, list)
	use linked_list
	implicit none

	real(dp), dimension(:), intent(in) :: Y
	type(linkedlist)  :: list
	type(llnode), pointer :: new

	!Make a node out of the array Y
	call ll_make(new,Y)

	!Append new node to list.
	call ll_append(new,list)

end subroutine rec_traj

!Prints the linked list and then deletes it.
subroutine print_del_traj(list, numb)
	use linked_list
	implicit none

	type(linkedlist), intent(inout) :: list
	type(llnode), pointer :: move
	integer, intent(in) :: numb
	integer :: i

	!Write to file.
	!Check if list is empty.
	if (.not. associated(list%head)) then
		print*, "The list is empty."
	else
		move=>list%head
		do
			if (allocated(move%a)) then
				write(unit=numb),(move%a(i),i=1,size(move%a))
			end if
			move=>move%next
			if (.not. associated(move)) exit
		end do
	end if
	!Delete the list.
	call ll_del_all(list)
end subroutine print_del_traj

!Subroutine which will initialize the zoom-in technique.  Finds a reference point from namelist and returns one point on the equal energy surface, a distance "tol" away from yref.  Uses the Euclidean metric on field space.
subroutine ic_zoom_init(y0, yref, iccounter, toler)
  use features, only : newunit
	implicit none

	real(dp), dimension(:), intent(inout) :: y0, yref
  real(dp), dimension(5) :: yinit
	real(dp), intent(inout) :: toler
	integer, intent(inout) :: iccounter
  integer :: u

	namelist / zoom / yinit, toler

	!Read the reference point "yinit" and tolerance "toler" from file.
	open(unit=newunit(u), file="parameters_hybrid.txt", status="old",&
	& delim = "apostrophe")
	read(unit=u, nml=zoom)
	close(unit=u)

  !Change yinit to proper form (phi,psi,...)
  yref(2)=yinit(3)
  yref(3)=yinit(2)
  yref(4)=yinit(5)
  yref(5)=yinit(4)

	!Get a new point on eq en slice a dist "toler" away from yref.
	call ic_eqen_pert(yref,y0,iccounter,euclidean,toler)

end subroutine ic_zoom_init

!Euclidean metric.
pure real(dp) function euclidean(pt1,pt2)
	implicit none

	real(dp), dimension(:), intent(in) :: pt1, pt2

	euclidean=sqrt(sum((pt1-pt2)*(pt1-pt2)))

end function euclidean

subroutine ic_init(ic, y0, iccounter, sample_table, burnin, rank,&
  &numtasks, ic_table, yref, toler, printing)
  implicit none

  real(dp), dimension(:), intent(out) :: y0, yref
  integer, intent(inout) :: ic, iccounter, rank, numtasks, burnin
  real(dp), dimension(:,:), allocatable, intent(out) :: sample_table, ic_table
  real(dp), intent(inout) :: toler
  logical, intent(in) :: printing

	if (ic == 1) then
		!Zero vel slice.
		call d_ic_zerov(y0)
	else if (ic == 2) then
		!Eq energy slice.
		call d_ic_eqen(y0,iccounter)
	else if (ic==3) then
		!Subslice of eq en slice.
		call eqen_slicing(y0)
	else if (ic==4) then
		!Metropolis sample a dataset.
    if (printing) print*, "Building Metropolis sample."
		call ic_metr_init(y0, iccounter, sample_table, burnin)
    if (printing) print*, "Burn in done."
	else if (ic==5) then
		!Read IC from a file.
    if (printing) print*,"Getting ICs from file."
		call ic_file_init(y0, rank,numtasks,ic_table)
    if (printing) print*, "IC table built."
	else if (ic==6) then
		!Zoom in on one point on eq en surface.
		call ic_zoom_init(y0, yref, iccounter, toler)
  else if (ic==7) then
    !Get a fixed initial condition.
    call fixed_ic(y0)
  else
    print*,"ERROR: IC out of range."
	end if

end subroutine ic_init


end module hybrid_subroutines

