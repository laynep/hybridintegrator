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
		OPEN(UNIT=sucunit,status='new',file=sucname, form=formt)
		OPEN(UNIT=failunit,status='new',file=failname, form=formt)
	else
		OPEN(UNIT=sucunit,status='new',file=sucname, form='unformatted')
		OPEN(UNIT=failunit,status='new',file=failname, form='unformatted')
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


	OPEN(unit=u,status='new',file='info.200')
	if(printing) PRINT*, "Hybrid Inflation -- Homogeneous"
	WRITE(unit=u,fmt=*) "Hybrid Inflation -- Homogeneous"
	IF (IC == 1) THEN
		if(printing) PRINT*,"ZERO VELOCITY SLICE"
		WRITE(unit=u,fmt=*) "ZERO VELOCITY SLICE"
	ELSE IF (IC == 2) THEN
		if(printing) PRINT*,"EQUAL ENERGY SLICE"
		WRITE(unit=u,fmt=*) "EQUAL ENERGY SLICE"
	ELSE IF (IC == 3) THEN
		if(printing) PRINT*,"2D SLICING OF EQEN SLICE"
                WRITE(unit=u,fmt=*) "2D SLICING OF EQEN SLICE"
	ELSE IF (IC == 4) THEN
		if(printing) PRINT*,"IC FROM METROPOLIS SAMPLING"
                WRITE(unit=u,fmt=*) "IC FROM METROPOLIS SAMPLING"
	END IF

	if(printing) PRINT*,"Planck mass is ",m_planck
	WRITE(unit=u,fmt=*) "Planck mass is ",m_planck
	if(printing) PRINT*,"E=",energy_scale
	WRITE(unit=u,fmt=*) "E=",energy_scale
	if(printing) PRINT*,"Lambda is ",lambda
	WRITE(unit=u,fmt=*) "Lambda is ",lambda
       	if(printing) PRINT*,"M is ",M
	WRITE(unit=u,fmt=*) "M is ",M
       	if(printing) PRINT*,"mu is ",mu
	WRITE(unit=u,fmt=*) "mu is ",mu
       	if(printing) PRINT*,"nu is ",nu
	WRITE(unit=u,fmt=*) "nu is ",nu

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
	IF (IC == 1) THEN
		if(printing) PRINT*,"ZERO VELOCITY SLICE"
		WRITE(unit=u,fmt=*) "ZERO VELOCITY SLICE"
	ELSE IF (IC == 2) THEN
		if(printing) PRINT*,"EQUAL ENERGY SLICE"
		WRITE(unit=u,fmt=*) "EQUAL ENERGY SLICE"
	ELSE IF (IC==3) THEN
		if(printing) PRINT*, "2D SLICING OF EQEN SLICE."
		WRITE(unit=u,fmt=*) "2D SLICING OF EQEN SLICE"
	ELSE IF (IC == 4) THEN
		if(printing) PRINT*,"IC FROM METROPOLIS SAMPLING"
                WRITE(unit=u,fmt=*) "IC FROM METROPOLIS SAMPLING"
	END IF
	if(printing) PRINT*,"Number of succ points ",counter,&
		&"Number of fail points ",failcount," ratio ",ratio
	WRITE(unit=u,fmt=*) "Number of succ points ",counter,&
		&"Number of fail points ",failcount," ratio ",ratio
	if(printing) PRINT*,"Number of IC that gave bad field values ",badfieldcounter
	WRITE(unit=u,fmt=*)"Number of IC that gave bad field values",&
		&badfieldcounter
	if(printing) PRINT*,"Number that didn't fall into minimum ",errorcount
	WRITE(unit=u,fmt=*) "Number that didn't fall into minimum ",errorcount


end subroutine hybrid_finalstats


subroutine new_point(y0,iccounter,sample_table,ic)
implicit none

	double precision, dimension(:), intent(inout) :: y0
	integer, intent(inout) :: iccounter, ic
	double precision, dimension(:,:), allocatable, intent(inout) :: sample_table


	IF (IC == 1) THEN
		CALL D_IC_ZEROV(Y0)
	ELSE IF (IC == 2) THEN
		CALL D_IC_EQEN(Y0,iccounter)
	ELSE IF (IC==3) THEN
		CALL EQEN_SLICING(Y0)
	ELSE IF (IC==4) THEN
		CALL IC_METR(Y0,sample_table,iccounter)
	END IF

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
		&+V - (lambda**4D0))

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
		IF (IER .NE. 0) THEN
			PRINT*, ' SUNDIALS_ERROR: FNVINITS returned IER = ', IER
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


























