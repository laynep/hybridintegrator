module hybrid_subroutines
use d_hybrid_initialconditions
implicit none

contains

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
	if (IC == 1) then
		if(printing) print*,"ZERO VELOCITY SLICE"
		write(unit=u,fmt=*) "ZERO VELOCITY SLICE"
	else if (IC == 2) then
		if(printing) print*,"EQUAL ENERGY SLICE"
		write(unit=u,fmt=*) "EQUAL ENERGY SLICE"
	else if (IC == 3) then
		if(printing) print*,"2D SLICING OF EQEN SLICE"
                write(unit=u,fmt=*) "2D SLICING OF EQEN SLICE"
	else if (IC == 4) then
		if(printing) print*,"IC FROM METROPOLIS SAMPLING"
                write(unit=u,fmt=*) "IC FROM METROPOLIS SAMPLING"
	end if

	if(printing) print*,"Planck mass is ",m_planck
	write(unit=u,fmt=*) "Planck mass is ",m_planck
	if(printing) print*,"E=",energy_scale
	write(unit=u,fmt=*) "E=",energy_scale
	if(printing) print*,"Lambda is ",lambda
	write(unit=u,fmt=*) "Lambda is ",lambda
       	if(printing) print*,"M is ",M
	write(unit=u,fmt=*) "M is ",M
       	if(printing) print*,"mu is ",mu
	write(unit=u,fmt=*) "mu is ",mu
       	if(printing) print*,"nu is ",nu
	write(unit=u,fmt=*) "nu is ",nu

end subroutine hybrid_initstats

subroutine hybrid_finalstats(ic, counter, failcount, badfieldcounter, &
		& errorcount, printing,  infounit)
implicit none

	integer, intent(inout) :: ic, counter, failcount, badfieldcounter, errorcount
	integer, optional, intent(in) :: infounit
	integer :: u
	double precision :: ratio
	logical, intent(in) :: printing

	if(present(infounit)) then
		u=infounit
	else
		u=200
	end if

	ratio=DBLE(counter)/DBLE(counter+failcount)
	if (IC == 1) then
		if(printing) print*,"ZERO VELOCITY SLICE"
		write(unit=u,fmt=*) "ZERO VELOCITY SLICE"
	else if (IC == 2) then
		if(printing) print*,"EQUAL ENERGY SLICE"
		write(unit=u,fmt=*) "EQUAL ENERGY SLICE"
	else if (IC==3) then
		if(printing) print*, "2D SLICING OF EQEN SLICE."
		write(unit=u,fmt=*) "2D SLICING OF EQEN SLICE"
	else if (IC == 4) then
		if(printing) print*,"IC FROM METROPOLIS SAMPLING"
                write(unit=u,fmt=*) "IC FROM METROPOLIS SAMPLING"
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


subroutine new_point(y0,iccounter,sample_table,ic)
implicit none

	double precision, dimension(:), intent(inout) :: y0
	integer, intent(inout) :: iccounter, ic
	double precision, dimension(:,:), allocatable, intent(inout) :: sample_table


	if (ic == 1) then
		call d_ic_zerov(y0)
	else if (ic == 2) then
		call d_ic_eqen(y0,iccounter)
	else if (ic==3) then
		call eqen_slicing(y0)
	else if (ic==4) then
		call ic_metr(y0,sample_table,iccounter)
	end if

end subroutine new_point




subroutine succ_or_fail(Y, success, successlocal, faillocal, &
		&sucunit, failunit, ic, leave, printing)
implicit none

	integer, intent(inout) :: success, successlocal, faillocal, sucunit, failunit, ic
	double precision, dimension(:), intent(inout) :: Y
	logical, intent(inout) :: leave
	double precision :: V, check
	logical, intent(in) :: printing


	!Potential function.
	V = V_h(Y)
	
	!Check if  fields fell into minima (if engy dnsty <
	! dnsty at infl end)
	check = (.5D0*( (Y(4)*Y(4)) + (Y(5)*Y(5))) &
		&+V - (lambda*lambda*lambda*lambda))

	!Gives a 1 if inflation is successful.
	leave = .false.
	if (Y(1)>65) then
		success = 1
		successlocal = successlocal + 1
		if (IC==1) then
			write(unit=sucunit),psi_0,phi_0
		else
			write(unit=sucunit), psi_0, phi_0,&
			& psi_dot_0, phi_dot_0
		end if
		if(printing .and. mod(successlocal,1000)==0) then
			print*,successlocal,success
		end if
		leave = .true.
	elseif (check<0D0) then
		success = 0
		faillocal = faillocal + 1
		if (IC==1) then
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


subroutine set_paramsFCVODE(rpar, neq, nglobal, numtasks, iatol, atol, rtol, &
		& meth, itmeth, t0, t, itask, tout)
implicit none

	!FCVODE PARAMS
	double precision, intent(out) :: RPAR(5)
	integer, intent(out) :: METH, ITMETH
	integer(kind=8), intent(out) :: NEQ, NGLOBAL
	integer, intent(out) :: IATOL, ITASK
	double precision, intent(out) :: T0, T, TOUT, RTOL, ATOL(5)

	!Other params
	integer, intent(in) :: numtasks	
	integer :: ier

	!Pass params to FCVODE subroutines.  Errors if don't do this...?
	RPAR(1) = lambda
	RPAR(2) = m
	RPAR(3) = mu
	RPAR(4) = nu
	RPAR(5) = beta

	!DOF global and local.
	NEQ = 5
	NGLOBAL = NEQ*numtasks

	!Set params for FCVODE.
	CALL FNVINITS(1, NEQ, IER)
		if (IER .NE. 0) then
			print*, ' SUNDIALS_ERROR: FNVINITS returned IER = ', IER
			STOP
		ENDIF
	IATOL = 2
	ATOL = 1D-10
	RTOL = 1D-10
	METH = 2
	ITMETH = 2
	T0=0D0
	T=T0
	ITASK = 1
	TOUT = 1D1

end subroutine set_paramsFCVODE

!Subroutine to check whether the integrator is giving any success points.  If check is present, then this will also check to see if the ratio of succ to fail points is less than check.
subroutine all_fail_check(successlocal, faillocal, allfailcheck, printing, check)
implicit none

	integer, intent(in) :: successlocal, faillocal
	logical, intent(inout) :: allfailcheck
	double precision, optional, intent(in) :: check
	double precision :: ratio
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


end module hybrid_subroutines


























