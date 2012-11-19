!**************************************************************************************
!Layne Price, University of Auckland, May 5, 2012

!A series of subroutines which set the initial conditions for hybrid inflation subject 
!to a set of constraints, dependent on each subroutine.

!Note that the parameters_hybrid() subroutine must be declared first, before the other 
!routines can be used.

!subroutine parameters_hybrid() ................... Sets the particle physics 
!						parameters for hybrid inflation.
!function V_h(y)................................... Gives the hybrid potential.
!subroutine D_IC_EQEN(y,iccounter) ................ Sets the ICs according to equal
!						energy slice. Field values are set
!						randomly and field velocities are 
!						alternately set on subsequent calls
!						by energy constraint.
!subroutine D_IC_ZEROV(y) ......................... Sets IC with zero velocity.  Clesse
!						original.

!y(1)=N, y(2)=phi, y(3)=psi, y(4)=phi_dot, y(5)=psi_dot.

!**************************************************************************************

module d_hybrid_initialconditions
  use types, only : dp, pi
  implicit none

	!Global data.
	real(dp) :: phi_0, psi_0, phi_dot_0, psi_dot_0
	real(dp), parameter :: m_planck=1000e0_dp
	real(dp), parameter :: beta=8.37758040957278e0_dp/(m_planck*m_planck)
	real(dp) :: m
	real(dp) :: mu

!	real(dp), parameter :: phi_min=44.9e0_dp
!	real(dp), parameter :: phi_max=45.1e0_dp
!	real(dp), parameter :: psi_min=8.9e0_dp
!	real(dp), parameter :: psi_max=9.1e0_dp
	real(dp), parameter :: phi_min=0e0_dp
	real(dp), parameter :: phi_max=.2e0_dp*m_planck
	real(dp), parameter :: psi_min=0e0_dp
	real(dp), parameter :: psi_max=.2e0_dp*m_planck
	real(dp) :: nu
	real(dp) :: energy_scale
	real(dp) :: lambda

	!For the nearest neighbor sampling.
	integer, dimension(:,:), allocatable :: numbpoints_sample
	integer :: xxglobal_fail, xxglobal_succ

	namelist /parameters/ m, mu, nu, lambda, energy_scale

contains

!*************************************************************************************
!Subroutine which declares potential parameters for hybrid inflation.

subroutine parameters_hybrid()
  use features, only : newunit
  implicit none

  integer :: u
  !Use the function newunit to be compatible with older versions of gfortran.
	open(unit=newunit(u), file="parameters_hybrid.txt", status="old", delim = "apostrophe")
	read(unit=u, nml=parameters)
	close(unit=u)

end subroutine parameters_hybrid

!*************************************************************************************
!function which describes the potential.  y must be (N,phi,psi,phi_dot,psi_dot)

pure real(dp) function V_h(y)
	implicit none

	real(dp), dimension(:), intent(in) :: y

	V_h = (lambda*lambda*lambda*lambda)*((1e0_dp - ((y(3)*y(3))/(m*m)))**2e0_dp +&
		& ((y(2)*y(2))/(mu*mu)) + ((y(2)*y(2)*y(3)*y(3))/ (nu*nu*nu*nu)))

end function V_h


!******************************************************************************************


!This subroutine will set the initial conditions on an equal energy slice at the start of inflation that corresponds to the value energy_scale which is specified in the parameter declaration in the main program.

subroutine d_ic_eqen(y,iccounter)
	implicit none

  real(dp), intent(out) :: y(5)
	integer, intent(in) :: iccounter
	real(dp) :: rand, rho_kinetic, dot_min, dot_max
	integer :: param_constr, a, b

	!Gives parameter to set by energy constraint: 0~phi_dot,1~psi_dot.  The IC for iccounter equals number of parameters to oscillate between.
	param_constr = MOD(iccounter,2)
	if (param_constr==0) then
		a=5
		b=4
	else
		a=4
		b=5
	end if

	do
		!Set phi from energy constraint.
		call random_number(rand)
		y(2)= (rand*(phi_max-phi_min)) + phi_min

		!IC for y(3)~psi randomly in range psi_min to psi_max.
		call random_number(rand)
		y(3)= (rand*(psi_max-psi_min)) + psi_min
	
		!Energy density remaining in kinetic term.
		rho_kinetic = (energy_scale**4e0_dp) - V_h(y)
		if (rho_kinetic<0) cycle
		!Set so that psi_dot is chosen with flat prior.
		dot_max = SQRT(2e0_dp*rho_kinetic)
		dot_min = -1e0_dp*SQRT(2e0_dp*rho_kinetic)	
		!Sets the psi_dot IC to the range psi_dot_max to psi_dot_min
		call random_number(rand)
		y(a) = (rand*(dot_max-dot_min)) + dot_min
	
		if(2e0_dp*rho_kinetic - (y(a)*y(a)) > 0) exit	
	end do
		
	!Set the phi_dot IC by the total energy density constraint.
	call random_number(rand)
	if(rand < .5) then
		y(b) = SQRT(2e0_dp*rho_kinetic - (y(a)*y(a)))
	else 
		y(b) = -1e0_dp*SQRT(2E0*rho_kinetic - (y(a)*y(a)))
	end if	

	!Set initial conditions.
	phi_0 = y(2)
	psi_0 = y(3)
	psi_dot_0 = y(5)
	phi_dot_0 = y(4)

	!Reinitialize the e-fold value: y(1)~N.
	y(1)=0e0_dp 		

end subroutine d_ic_eqen

!*********************************************************************************
!Subroutine which will take one pt y0 and give another point y1 that is very close to it, but also on the equal energy slice.  

subroutine ic_eqen_pert(y0,y1,iccounter,metric,sig,en)
	implicit none

	real(dp), dimension(5), intent(in) :: y0
	real(dp), dimension(5), intent(out) :: y1
	interface
		pure function metric(pt1,pt2)
      use types, only : dp
			implicit none
			real(dp), dimension(:), intent(in) :: pt1, pt2
			real(dp) :: metric
		end function metric
	end interface
	real(dp), optional, intent(in) :: en, sig
	integer, intent(in) :: iccounter
	real(dp) :: rand, rho_kinetic, e, tol, sgn
	integer :: param_constr, x, y, i
	real(dp), dimension(5) :: maxim, minim

	!Set energy.
	if (present(en)) then
		e=en
	else
		e=energy_scale
	end if
	!Set tolerance
	if (present(sig)) then
		tol=sig
	else
		tol=1e-12_dp
	end if

	!Set the max and min params.
	maxim=(y0+tol/2e0_dp)
	minim=(y0-tol/2e0_dp)

	!Initialize y1.
	y1=y0
	!Gives parameter to set by energy constraint: 0~phi_dot,1~psi_dot.  
	param_constr = MOD(iccounter,2)
	if (param_constr==0) then
		x=4
		y=5
	else
		x=5
		y=4
	end if

  do1:	do
    do2:	do
			!set phi from energy constraint.
			call random_number(rand)
			y1(2)= (rand*(maxim(2)-minim(2))) + minim(2)
			!ic for y(3)~psi randomly in range psi_min to psi_max.
			call random_number(rand)
			y1(3)= (rand*(maxim(3)-minim(3))) + minim(3)
			!energy density remaining in kinetic term.
			rho_kinetic = (e**4) - v_h(y1)
			if (rho_kinetic<0) cycle 
			!sets the psi_dot ic to the range psi_dot_max to psi_dot_min
			call random_number(rand)
			y1(y) = (rand*(maxim(y)-minim(y))) + minim(y)
			if(2e0_dp*rho_kinetic - (y1(y)*y1(y)) > 0) exit do2
		end do do2
		!Set the phi_dot IC by the total energy density constraint.
		sgn=y0(x)/abs(y0(x))
		y1(x) = sgn*sqrt(2e0_dp*rho_kinetic - (y1(y)*y1(y)))
		!Exit condition.
		if (metric(y0(2:5),y1(2:5)) .le. tol ) exit do1
	end do do1

	!Set initial conditions.
	phi_0 = y1(2)
	psi_0 = y1(3)
	psi_dot_0 = y1(5)
	phi_dot_0 = y1(4)

	!Reinitialize the e-fold value: y(1)~N.
	y1(1)=0e0_dp 		

end subroutine ic_eqen_pert


!*********************************************************************************


!Program subroutine: This subroutine will follow the results of Clesse and set the initial values of the velocity parameters to zero and choose the initial conditions of the fields at random.

subroutine d_ic_zerov(y)
	implicit none

  real(dp), intent(out) :: y(5)
	real(dp) :: rand_1, rand_2

	!Set IC for y(2)~phi randomly in range phi_min to phi_max.
	call random_number(rand_1)
	y(2)= (rand_1*(phi_max-phi_min)) + phi_min
	phi_0 = y(2)

	!Set IC for y(3)~psi randomly in range psi_min to psi_max.
	call random_number(rand_2)
	y(3)= (rand_2*(psi_max-psi_min)) + psi_min
	psi_0 = y(3)

	!Sets velocities to zero initially.
	y(4)=0e0_dp
	phi_dot_0=y(4)
	y(5)=0e0_dp
	psi_dot_0=y(5)

	!Reinitialize the cosmology equation values.
	y(1)=0e0_dp 
		
	
end subroutine d_ic_zerov

!*******************************************************************************************
!Fixed initial conditions set here.

subroutine fixed_ic(y)
  implicit none

	real(dp), dimension(5), intent(out) :: y

	phi_0 = 29.092534874291754_dp
	psi_0 = 64.663923706055400_dp
	phi_dot_0 =46.619723821534521_dp
	psi_dot_0 =133.51629619827270_dp
	y(1)=0e0_dp
	y(2)=phi_0
	y(3)=psi_0
	y(4)=phi_dot_0
	y(5)=psi_dot_0

end subroutine fixed_ic

!****************************************************************************
!Subroutine that slices the eqen surface.  This particular choice sets psi_0=0 so that all the initial conditions are set in (if the velocities were then zero) the inflationary valley.  psi_dot_0 is set by the energy constraint so that we can plot phi_0 vs phi_dot_0 as a two-dimensional representation of the three-dimensional slice.

subroutine eqen_slicing(y)
	implicit none

	real(dp), dimension(:), intent(out) :: y
	real(dp) :: rand_1, rand_2, rand_3, chi
	real(dp) :: phi_dot_max, phi_dot_min, rho_kinetic, psi_dot_max, psi_dot_min

	!Initialize.
	y=0e0_dp
	do
   	!Set phi_0 randomly.
    call random_number(rand_1)
    y(2)=(rand_1*(phi_max-phi_min)) + phi_min
    phi_0 = y(2)

    !Set psi_0 randomly.
    call random_number(rand_1)
    y(3)=(rand_1*(psi_max-psi_min)) + psi_min
    psi_0 = y(3)

    !Find remaining energy.
    rho_kinetic = (energy_scale**4e0_dp) - V_h(y)

    if (rho_kinetic<0) cycle

    chi = -1e0_dp*SQRT(rho_kinetic)
    y(4)=chi
    phi_dot_0=chi
    y(5)=chi
    psi_dot_0=chi

    exit

   end do


end subroutine eqen_slicing

subroutine eqen_eqvel(y)
	implicit none

	real(dp), dimension(:), intent(out) :: y
	real(dp) :: rand_1, rho_kinetic, chi
	
	!Initialize
	y=0e0_dp

	do
		!Set phi_0 randomly.
    call random_number(rand_1)
    y(2)=(rand_1*(phi_max-phi_min)) + phi_min
    phi_0 = y(2)

		!Set psi_0 randomly.
    call random_number(rand_1)
    y(3)=(rand_1*(psi_max-psi_min)) + psi_min
    psi_0 = y(3)

		!Find remaining energy.
    rho_kinetic = (energy_scale**4e0_dp) - V_h(y)

		if (rho_kinetic<0) cycle

		chi = sqrt(rho_kinetic)
		y(4)=chi
		phi_dot_0=chi
		y(5)=chi
		psi_dot_0=chi

		exit

	end do

end subroutine eqen_eqvel

!Gives uniform sample of IC from 0-1.  GOOD FOR TESTING.
subroutine ic_test(y)
	implicit none

	real(dp), dimension(:), intent(out) :: y
	real(dp) :: rand_1
	integer :: i

	do i=1,size(y)
		call random_number(rand_1)
		y(i)=rand_1
	end do

end subroutine ic_test

!Subroutine to open a file to read ICs from.  File has "length" x "dimn"-dimensional points, which are both passed as arguments to this routine.  Reads into temporary array from rank=0 (master), then scatters pieces to other threads in variable ic_table.
subroutine readdist_icfromfile(rank, numtasks, ic_table, fname, formt, length, dimn)
	use features, only : newunit
	use mpi
	implicit none

	real(dp), dimension(:,:), allocatable :: mastertable, masttransp
	real(dp), dimension(:,:), allocatable, intent(out) :: ic_table
	real(dp), dimension(:,:), allocatable :: ic_transp
	integer, intent(in) :: rank, numtasks, length, dimn
	character(len=*), intent(in) :: fname, formt
	integer :: i, j, n, low, receiv
	integer :: ierr, u

	!Make the ic_table for each thread.
	!Find how many ICs to give to each thread.  Note, some will be lost bc
	!I want to give each thread the same number of ICs.
	n=length/numtasks	!Note int div.
	allocate(ic_table(n,dimn),ic_transp(dimn,n))
	allocate(mastertable(length,dimn),masttransp(dimn,length))

	!Read ICs into mastertable by master thread.
	!NOTE: we have to do the kludgy transpose mess because the file we read from
	!will be in row-major order and scatter works in column-major order for Fortran.
	if (rank==0) then
		!Read.
		open(unit=newunit(u),file=fname,form=formt,status="old")
		do i=1,length
			read(u), (mastertable(i,j),j=1,dimn)
		end do
		close(u)
		masttransp=transpose(mastertable)
		deallocate(mastertable)
	end if
	!Scatter mastertransp to other threads.
	call mpi_scatter(masttransp,n*dimn,mpi_double_precision,&
	&ic_transp,n*dimn,mpi_double_precision,0, mpi_comm_world, ierr)

	!Halt processors until master sends them ic_table.
 	call mpi_barrier(mpi_comm_world,ierr)
	!Get transpose of ictable.
	ic_table=transpose(ic_transp)

	!Deallocate my arrays.
	if(allocated(ic_transp)) deallocate(ic_transp)

end subroutine readdist_icfromfile

!Subroutine to get the xth initial condition from an array ic_table.
subroutine ic_fromarray(y,ic_table,x)
	implicit none

	real(dp), dimension(:,:), intent(in) :: ic_table
	real(dp), dimension(5), intent(out) :: y
	integer, intent(in) :: x

	y(2)=ic_table(x,2)
	y(3)=ic_table(x,1)
	y(4)=ic_table(x,4)
	y(5)=ic_table(x,3)	
	y(1)=0e0_dp
	phi_0 = y(2)
	psi_0 = y(3)
	psi_dot_0 = y(5)
	phi_dot_0 = y(4)


end subroutine ic_fromarray


!Subroutine that initializes the file to read from.
subroutine ic_file_init(y0, rank,numtasks,ic_table)
	use features, only : newunit
  implicit none

	real(dp), dimension(:), intent(out) :: y0
	integer, intent(in) :: rank, numtasks
	real(dp), dimension(:,:), allocatable, intent(out) :: ic_table
	character(len=100) :: fname
	character(len=100) :: fform
	integer :: fleng, fwid, u

	namelist /filetoread/ fname, fleng, fwid, fform

	!Open parameter file and read namelist.
	open(unit=newunit(u),file="parameters_hybrid.txt", status="old", delim = "apostrophe")
	read(unit=u, nml=filetoread)
	close(unit=u)
	fname=trim(fname)
	fform=trim(fform)

	!Read ics into local ictable.
	call readdist_icfromfile(rank, numtasks, ic_table, fname, &
		& fform, fleng, fwid)

	!Get first ic.
	call ic_fromarray(y0,ic_table,1)


end subroutine ic_file_init

!Calculate the energy of a vector --- for testing.
pure real(dp) function energy_vector(y)
  implicit none

  real(dp), dimension(:), intent(in) :: y

  energy_vector=(.5e0_dp*(y(4)*y(4)+y(5)*y(5))+v_h(y))**(.25e0_dp)

end function energy_vector


end module d_hybrid_initialconditions

