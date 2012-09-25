!**************************************************************************************
!Layne Price, University of Auckland, May 5, 2012

!A series of subroutines which set the initial conditions for hybrid inflation subject to a set of constraints, dependent on each subroutine.

!Note that the parameters_hybrid() subroutine must be declared first, before the other routines can be used.

!subroutine parameters_hybrid() ................... Sets the particle physics 
!						parameters for hybrid inflation.
!function V_h(Y)................................... Gives the hybrid potential.
!subroutine D_IC_EQEN(Y,iccounter) ................ Sets the ICs according to equal
!						energy slice. Field values are set
!						randomly and field velocities are 
!						alternately set on subsequent calls
!						by energy constraint.
!subroutine D_IC_ZEROV(Y) ......................... Sets IC with zero velocity.  Clesse
!						original.

!Y(1)=N, Y(2)=phi, Y(3)=psi, Y(4)=phi_dot, Y(5)=psi_dot.

!**************************************************************************************

module d_hybrid_initialconditions
  use types, only : dp, pi
  implicit none

	!global data.
	real(dp) :: phi_0, psi_0, phi_dot_0, psi_dot_0
	real(dp), parameter :: m_planck=1000_dp
	real(dp), parameter :: beta=8.37758040957278_dp/(m_planck*m_planck)
	real(dp) :: m
	real(dp) :: mu
	real(dp), parameter :: phi_min=0_dp
	real(dp), parameter :: phi_max=.2_dp*m_planck
	real(dp), parameter :: psi_min=0_dp
	real(dp), parameter :: psi_max=.2_dp*m_planck
	real(dp) :: nu
	real(dp) :: energy_scale
	real(dp) :: lambda

	!for the nearest neighbor sampling.
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
!function which describes the potential.

pure real(dp) function V_h(Y)
	implicit none

	real(dp), dimension(:), intent(in) :: Y

	V_h = (lambda*lambda*lambda*lambda)*((1_dp - ((Y(3)*Y(3))/(m*m)))**2_dp +&
		& ((Y(2)*Y(2))/(mu*mu)) + ((Y(2)*Y(2)*Y(3)*Y(3))/ (nu*nu*nu*nu)))

end function V_h


!******************************************************************************************


!This subroutine will set the initial conditions on an equal energy slice at the start of inflation that corresponds to the value energy_scale which is specified in the parameter declaration in the main program.

subroutine D_IC_EQEN(Y,iccounter)
	implicit none

  real(dp), intent(out) :: Y(5)
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
		Y(2)= (rand*(phi_max-phi_min)) + phi_min

		!IC for Y(3)~psi randomly in range psi_min to psi_max.
		call random_number(rand)
		Y(3)= (rand*(psi_max-psi_min)) + psi_min
	
		!Energy density remaining in kinetic term.
		rho_kinetic = (energy_scale**4_dp) - V_h(Y)
		if (rho_kinetic<0) cycle
		!Set so that psi_dot is chosen with flat prior.
		dot_max = SQRT(2_dp*rho_kinetic)
		dot_min = -1_dp*SQRT(2_dp*rho_kinetic)	
		!Sets the psi_dot IC to the range psi_dot_max to psi_dot_min
		call random_number(rand)
		Y(a) = (rand*(dot_max-dot_min)) + dot_min
	
		if(2_dp*rho_kinetic - (Y(a)*Y(a)) > 0) exit	
	end do
		
	!Set the phi_dot IC by the total energy density constraint.
	call random_number(rand)
	if(rand < .5) then
		Y(b) = SQRT(2_dp*rho_kinetic - (Y(a)*Y(a)))
	else 
		Y(b) = -1_dp*SQRT(2E0*rho_kinetic - (Y(a)*Y(a)))
	end if	

	!Set initial conditions.
	phi_0 = Y(2)
	psi_0 = Y(3)
	psi_dot_0 = Y(5)
	phi_dot_0 = Y(4)

	!Reinitialize the e-fold value: Y(1)~N.
	Y(1)=0_dp 		

end subroutine D_IC_EQEN

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
	maxim=(y0+tol/2_dp)
	minim=(y0-tol/2_dp)

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
			if(2_dp*rho_kinetic - (y1(y)*y1(y)) > 0) exit do2
		end do do2
		!Set the phi_dot IC by the total energy density constraint.
		sgn=y0(x)/abs(y0(x))
		Y1(x) = sgn*sqrt(2_dp*rho_kinetic - (Y1(y)*Y1(y)))
		!Exit condition.
		if (metric(y0(2:5),y1(2:5)) .le. tol ) exit do1
	end do do1

	!Set initial conditions.
	phi_0 = Y1(2)
	psi_0 = Y1(3)
	psi_dot_0 = Y1(5)
	phi_dot_0 = Y1(4)

	!Reinitialize the e-fold value: Y(1)~N.
	Y1(1)=0_dp 		

end subroutine ic_eqen_pert




!*********************************************************************************


!Program subroutine: This subroutine will follow the results of Clesse and set the initial values of the velocity parameters to zero and choose the initial conditions of the fields at random.

subroutine D_IC_ZEROV(Y)
	implicit none

  real(dp), intent(out) :: Y(5)
	real(dp) :: rand_1, rand_2

	!Set IC for Y(2)~phi randomly in range phi_min to phi_max.
	call random_number(rand_1)
	Y(2)= (rand_1*(phi_max-phi_min)) + phi_min
	phi_0 = Y(2)

	!Set IC for Y(3)~psi randomly in range psi_min to psi_max.
	call random_number(rand_2)
	Y(3)= (rand_2*(psi_max-psi_min)) + psi_min
	psi_0 = Y(3)

	!Sets velocities to zero initially.
	Y(4)=0_dp
	phi_dot_0=Y(4)
	Y(5)=0_dp
	psi_dot_0=Y(5)

	!Reinitialize the cosmology equation values.
	Y(1)=0_dp 
		
	
end subroutine D_IC_ZEROV

!*******************************************************************************************
!Fixed initial conditions set here.

subroutine FIXED_IC(Y)
  implicit none

	real(dp), dimension(5), intent(out) :: Y

	phi_0 = 0_dp
	psi_0 = 4_dp*m_planck
	phi_dot_0 =0_dp
	psi_dot_0 =0_dp
	Y(1)=0_dp
	Y(2)=phi_0
	Y(3)=psi_0
	Y(4)=phi_dot_0
	Y(5)=psi_dot_0

end subroutine FIXED_IC

!****************************************************************************
!Subroutine that slices the eqen surface.  This particular choice sets psi_0=0 so that all the initial conditions are set in (if the velocities were then zero) the inflationary valley.  psi_dot_0 is set by the energy constraint so that we can plot phi_0 vs phi_dot_0 as a two-dimensional representation of the three-dimensional slice.

subroutine EQEN_SLICING(Y)
	implicit none

	real(dp), dimension(:), intent(out) :: Y
	real(dp) :: rand_1, rand_2, rand_3, chi
	real(dp) :: phi_dot_max, phi_dot_min, rho_kinetic, psi_dot_max, psi_dot_min

	!Initialize.
	Y=0_dp
	do
   	!Set phi_0 randomly.
    call random_number(rand_1)
    Y(2)=(rand_1*(phi_max-phi_min)) + phi_min
    phi_0 = Y(2)

    !Set psi_0 randomly.
    call random_number(rand_1)
    Y(3)=(rand_1*(psi_max-psi_min)) + psi_min
    psi_0 = Y(3)

    !Find remaining energy.
     rho_kinetic = (energy_scale**4_dp) - V_h(Y)

    if (rho_kinetic<0) then
			cycle
		end if

    chi = -1_dp*SQRT(rho_kinetic)
    Y(4)=chi
    phi_dot_0=chi
    Y(5)=chi
    psi_dot_0=chi

    exit

   end do


end subroutine EQEN_SLICING

subroutine EQEN_EQVEL(Y)
	implicit none

	real(dp), dimension(:), intent(out) :: Y
	real(dp) :: rand_1, rho_kinetic, chi
	
	!Initialize
	Y=0_dp

	do
		!Set phi_0 randomly.
    call random_number(rand_1)
    Y(2)=(rand_1*(phi_max-phi_min)) + phi_min
    phi_0 = Y(2)

		!Set psi_0 randomly.
    call random_number(rand_1)
    Y(3)=(rand_1*(psi_max-psi_min)) + psi_min
    psi_0 = Y(3)

		!Find remaining energy.
    rho_kinetic = (energy_scale**4_dp) - V_h(Y)

		if (rho_kinetic<0) cycle

		chi = SQRT(rho_kinetic)
		Y(4)=chi
		phi_dot_0=chi
		Y(5)=chi
		psi_dot_0=chi

		exit

	end do


end subroutine EQEN_EQVEL



!subroutine that will take as input a point Y0 and a sample table, sample_table.  This subroutine will perform a nearest neighbor interpolation on sample_table with respect to the density of points, within a given box size, eps.  By default this will not give duplicates, i.e. Y_init != Y_fin.  If you want to override this, specify dup=1 .

subroutine ic_metr(y,sample_table,iccounter,eps,dup,test)
	use sorters, only : locate, heapsorttotal
	implicit none

	real(dp), dimension(:), intent(inout) :: y
	real(dp), dimension(:,:), allocatable, intent(inout) :: sample_table
	real(dp), optional, intent(inout) :: eps
	integer, intent(in) :: iccounter
	integer, optional, intent(in) :: test
	integer, optional, intent(in) :: dup
	real(dp) :: accept, rand
	real(dp), dimension(5) :: yprop
	integer, dimension(:,:), allocatable :: boxcover
	integer :: i,j,k, n, start
	logical :: check

	!Load eps if not provided.
	if (.not. present(eps)) eps=1_dp

	!Load the density calculation if not already done.
	if (.not. allocated(numbpoints_sample)) then
		!box cover, of size eps.
		allocate(boxcover(size(sample_table,1),size(sample_table,2)))
		boxcover=ceiling(sample_table/eps)
		if (allocated(sample_table)) deallocate(sample_table)
		!count numb of boxes.
		call heapsorttotal(boxcover)
		n=1
		do i=2,size(boxcover,1)
			do j=1,size(boxcover,2)
				check=.true.
				if (boxcover(i,j).ne.boxcover(i-1,j)) then
					check=.false.
					exit
				end if
			end do
			if (check) cycle
			n=n+1			
		end do
		!Allocate numbpoints_sample.  Each row will be Y(2),...,Y(5),# of Boxes.
		allocate(numbpoints_sample(n,5))
		!Load first elt of numbpoints_sample
		do i=1,size(boxcover,2)
			numbpoints_sample(1,i)=boxcover(1,i)*eps
		end do
		numbpoints_sample(1,5)=1
		!Count elts per box. Copy to numbpoints_sample.
		start=2
doi:		do i=1,n
	doj:		do j=start,size(boxcover,1)
		dok1:		do k=1,size(boxcover,2)
					check=.true.
					if (boxcover(j,k).ne.boxcover(j-1,k)) then
						check=.false.
						exit dok1
					end if
				end do dok1
				if (check) then
					numbpoints_sample(i,5)=numbpoints_sample(i,5)+1
				else
		dok2:			do k=1,size(boxcover,2)
						numbpoints_sample(i+1,k)=boxcover(j,k)
					end do dok2
					numbpoints_sample(i+1,5)=1
					start=j+1
					exit doj
				end if
			end do doj
		end do doi
		
		!Deallocate boxcover & sample_table
		if (allocated(boxcover)) deallocate(boxcover)
		if (allocated(sample_table)) deallocate(sample_table)
	end if
	do
		!GET A NEW POINT ON EQEN SLICE.
		if (.not. present(test)) then 
			call D_IC_EQEN(YPROP,iccounter)
		else
			!TEST!!!!
			call EQEN_SLICING(YPROP)
		end if

		!CALC ACCEPT RATIO
		accept=accept_ratio(Y,YPROP,eps)

		!GEN RAND NUMB
		call random_number(rand)
	
		!MOVE TO NEW POINT if RAND<A
		if (rand<accept) then
			xxglobal_succ=xxglobal_succ+1
			Y=YPROP
			exit
		end if
		xxglobal_fail=xxglobal_fail+1

		!if NO DUPLICATES, then cycle UNTIL GET UNIQUE NUMB, OTHERWISE RETURN
		if (present(dup)) exit		
	end do

end subroutine ic_metr


!Subroutine that initializes the IC_METR routine.  This loads the sample_table from a file and does a burn in period.  Eps is the size of the coarse-graining for the interpolant.
subroutine ic_metr_init(y, iccounter, sample_table, bperiod, eps)
  use features, only : newunit
	implicit none

	real(dp), dimension(:), intent(out) :: y
	real(dp), optional, intent(inout) :: eps
	integer, intent(inout) :: iccounter
	integer, optional, intent(in) :: bperiod
	real(dp), dimension(:,:), allocatable, intent(out) :: sample_table
	character(len=100) :: datafile
	integer :: samp_len, samp_wid
	integer :: i, j, iend, ierr, u

	namelist /sample/ samp_len, samp_wid, datafile

	!Set IC at random on EQEN slice.
	call D_IC_EQEN(Y,iccounter)

	!Get info on sample table from namelist.
	open(unit=newunit(u), file="parameters_hybrid.txt", status="old",&
	& delim = "apostrophe")
	read(unit=u, nml=sample)
	close(unit=u)
	datafile=trim(datafile)
		
	!Load the data to sample via nearest neighbor interpolation.
	allocate(sample_table(samp_len,samp_wid))
	sample_table=0_dp
	!Data *must* be in form Y(2),...,Y(5), where Y(1)=0_dp assumed.
	open(unit=newunit(u), file=datafile, status="old", form="unformatted")
	do i=1,size(sample_table,1)		
		read(unit=u,iostat=ierr) (sample_table(i,j),j=2,5)
		if (is_iostat_end(ierr)) exit
	end do
	close(unit=u)

	!Burn in.
	if (present(bperiod)) then
		iend = bperiod
	else
		iend=10000
	end if
	do i=1, iend
		!Get new IC from sample_table.
		iccounter=iccounter+1
		if (present(eps)) then
			call IC_METR(Y,sample_table,iccounter,eps)
		else
			call IC_METR(Y,sample_table,iccounter)
		end if
	end do

end subroutine ic_metr_init



!function to calculate the acceptance ratio for the Metropolis algorithm for the density function obtained in IC_METR.  Y1 is old point, Y2 is new point.

real(dp) function accept_ratio(Y1,Y2,eps)
	use sorters, only : locate
	implicit none

	real(dp), dimension(:), intent(in) :: Y1, Y2
	integer :: start, ending, i, j
	real(dp), intent(in) :: eps
	real(dp) :: a, p1,p2
	integer, dimension(5) :: test
	LOGICAL :: check

	!Calc probabilities.
	!Find where to start in table. Returns first value below Y1(2).
	call locate(numbpoints_sample,Y1(2)/eps,start)
	call locate(numbpoints_sample,Y1(2)/eps+1_dp,ending)
	p1=0_dp

	do i=start,ending
		test=(/0,numbpoints_sample(i,1),numbpoints_sample(i,2),&
			&numbpoints_sample(i,3),numbpoints_sample(i,4) /)

	doj1:	do j=1,size(Y1)
			check=.TRUE.
			if (CEILING(Y1(j)/eps).ne.test(j)) then
				check=.FALSE.
				exit doj1
			end if
		end do doj1
		if (check) then
			p1=numbpoints_sample(i,5)*1_dp
			exit
		end if
	end do
	call locate(numbpoints_sample,Y2(2)/eps,start)
	call locate(numbpoints_sample,Y2(2)/eps+1_dp,ending)
	p2=0_dp
	do i=start,ending
		test=(/0, numbpoints_sample(i,1),numbpoints_sample(i,2),&
			&numbpoints_sample(i,3),numbpoints_sample(i,4) /)
	doj2:	do j=1,size(Y2)
			check=.TRUE.
			if (CEILING(Y2(j)/eps) .NE. test(j)) then
				check=.FALSE.
				exit doj2
			end if

		end do doj2

		if (check) then
			p2=numbpoints_sample(i,5)*1_dp
			exit
		end if
	end do

	if (p2<=1.1_dp) then
		a=0_dp
	!Avoid div by zero
	else if (p1<=0e-10_dp) then
		a=1_dp
	else
		a=p2/p1
	end if
  !MIN(1_dp,a)
  if (a<1_dp) then
    accept_ratio=a
  else
    accept_ratio=1_dp
  end if

end function accept_ratio



!Gives uniform sample of IC from 0-1.  GOOD FOR TESTING.
subroutine IC_TEST(Y)
	implicit none

	real(dp), dimension(:), intent(out) :: Y
	real(dp) :: rand_1
	integer :: i

	do i=1,size(Y)
		call random_number(rand_1)
		Y(i)=rand_1
	end do

end subroutine IC_TEST

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
	!NOTE: we have to do the kludgy transpose mess, because the file we read from
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
	&ic_transp,n*dimn,mpi_double_precision,0, MPI_COMM_WORLD, ierr)

	!Halt processors until master sends them ic_table.
 	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
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
	y(1)=0_dp
	phi_0 = Y(2)
	psi_0 = Y(3)
	psi_dot_0 = Y(5)
	phi_dot_0 = Y(4)


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



!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!***************************************************************************
!*********************
!*********************
!****************************************************************************
!NOT WORKING!!!

!A routine to generate an MCMC Gibbs-like sample, ICtable(points,9), from the whole parameter space (phi_0,psi_0,phi_dot_0,psi_dot_0,Lambda,m,mu,nu) subject to the equal energy constraint.  This involves embedding the 7-D constraint space in the 8-D parameter space, doing Gibbs sampling in the 8-D space, then projecting this perturbation to the constraint surface.  As the flat version of this problem is statistically represented by a uniform distribution, it will not be necessary to reject any points unless the projection does not lie on the equal energy surface.

!Needs call to random_seed per processor before using.  Needs to declare ICtable in main program with dimension #points_desired x 9.


!subroutine mcmc_eqen(ICtable)
!	implicit none

!	real(dp), dimension(:,:), intent(inout) :: ICtable
!	real(dp) :: Y(5), Y2(5), a(9), norm(8), pert(8), params(4), params2(4)
!	real(dp) :: lam_4, V, mag, KE, en_left, check
!	real(dp) :: sigma(8)
!	real(dp) :: tol
!	integer :: iccounter, i, j, k, k_1, k_2, l, i_1, counter
!	LOGICAL :: onsurf, bounds

	!Load tolerance and normal's sigma.
!	tol = 1e-6_dp
!	sigma(1) = 10_dp			!Phi
!	sigma(2) = 10_dp			!Psi
!	sigma(3) = 100_dp		!Phi_dot
!	sigma(4) = 100_dp		!Psi_dot
!	sigma(5) = 100_dp		!Lambda
!	sigma(6) = 10_dp			!M
!	sigma(7) = 10_dp			!mu
!	sigma(8) = 10_dp			!nu

	!First pick a point on eq en surface.
!	iccounter = 0
!	counter = 0
!	call parameters_hybrid()
	!Load params vect.
!	params(1)=Lambda
!	params(2)=M
!	params(3)=mu
!	params(4)=nu	

!	call SEED_IC_EQEN(Y)
!print*,"Y0",Y
!print*,"PARAMS0",params

!	do WHILE (counter < size(ICtable,1)+10000)

		!Define the normal vector at Y.
		!a_1*(f_dot)+a_2*(s_dot)+a_3*(f)+a_4*(s)+
		!a_5(l)+a_6*(m)+a_7*(z)+a_8*(n)+a_9=0
!		lam_4 = params(1)**4
!		V = lam_4*((1_dp-((Y(3)*Y(3))/(params(2)*params(2))))**2_dp+&
!			&((Y(2)*Y(2))/(params(3)*params(3)))+&
!			&((Y(2)*Y(2)*Y(3)*Y(3))/(params(4)**4_dp)))
!		a(1) = Y(4)
!		a(2) = Y(5)
!		a(3) = 2_dp*lam_4*Y(2)*(1_dp/(mu*mu)+(Y(3)*Y(3))/(nu*nu*nu*nu))
!		a(4) = 2_dp*lam_4*Y(3)*((Y(2)*Y(2))/(nu*nu*nu*nu)-&
!			&(2_dp*Y(3)*Y(3))/(M*M*M*M)-2_dp/(M*M))
!		a(5) = (4_dp/lambda)*V
!		a(6) = ((4_dp*lam_4*Y(2)*Y(2))/(M*M*M))*(1_dp + (Y(2)*Y(2))/(M*M))
!		a(7) = -2_dp*lam_4*Y(2)*Y(2)/(mu*mu*mu)
!		a(8) = -4_dp*lam_4*Y(2)*Y(2)*Y(3)*Y(3)/(nu*nu*nu*nu*nu)
!		a(9) = .5_dp*(Y(4)*Y(4)+Y(5)*Y(5)) + V - energy_scale**4
		!Normalize.
!		mag = SQRT(a(1)*a(1)+a(2)*a(2)+a(3)*a(3)+a(4)*a(4)+a(5)*a(5)+&
!			&a(6)*a(6)+a(7)*a(7)+a(8)*a(8))
!		do i=1,8
!			norm(i) = a(i)/mag
!		end do

		!Perturb each parameter by Gaussian with mean=0.
!		do j=1,8
!			pert(j) = normal(0_dp, sigma(j))
!		end do

		!Create new point.
!		do k=2,5	
!			Y2(k) = Y(k) + pert(k-1)
!		end do
!		do l=1,4
!			params2(l) = params(l)+pert(l+4)
!		end do

		!Project onto the eq en slice.
!print*,"Calling proj."
!		call proj_eqen(Y2,Y,params2,params,norm,tol)

		!Double check if on eq en surface
!		V = (params2(1)**4)*((1_dp-((Y2(3)*Y2(3))/(params2(2)*params2(2))))**2_dp +&
!			&((Y2(2)*Y2(2))/(params2(3)*params2(3))) + &
!			&((Y2(2)*Y2(2)*Y2(3)*Y2(3))/(params2(4)**4)))
!		KE = .5_dp*(Y2(4)*Y2(4)+Y2(5)*Y2(5))
!		en_left = energy_scale**4 - KE - V
!		if(ABS(en_left) .LE. tol) then
!			onsurf = .TRUE.
!		else
!			onsurf = .FALSE.
!print*,"FAIL not on surf"
!			cycle
!		end if
		!Test if within bounds. 
!		if(Y2(2)<phi_min .OR. Y2(2)>phi_max .OR. Y2(3)<phi_min .OR. Y2(3)>phi_max &
!			& .OR. params2(1)<lambda_min .OR. params2(1)>lambda_max &
!			&.OR. params2(2)<m_min .OR. params2(2)>m_max .OR. &
!			&params2(3)<mu_min .OR. params2(3)>mu_max .OR.&
!			&params2(4)<nu_min .OR. params2(4)>nu_max) then
!			bounds = .FALSE.
!print*,"FAIL not in bounds"
!			cycle
!		else
!			bounds = .TRUE.
!		end if
		!If new point on surf, upgrade to new point.
!print*,"SUCCESS"
!		counter = counter + 1
!		Y=Y2
!		params = params2
		!Burn in period.
!		if(counter>10000) then
!			do k_1=1,5
!				ICtable(counter-10000,k_1)=Y(k_1)
!			end do
!			do k_2=1,4
!				ICtable(counter-10000,k_2+5)=params(k_2)
!			end do
!		end if
		
!	end do
	


!end subroutine mcmc_eqen




!**********************************************************************
!Function that takes an m-D vector and gives the projection of that vector onto the (m-1)-D tangent space with normal, norm.  If the normal is a unit coordinate, then the normal direction will have a value of zero.  Euclidean.

function projection(vect, norm)
	implicit none

	real(dp), dimension(:) :: vect, norm
	real(dp), dimension(size(vect)) :: projection
	
	projection = vect - norm*dot_product(vect,norm)

end function projection



!*************************************************************************
!NOTE: Not working properly...

!Subroutine which takes the perturbed vector Y=Y0+pert and params=params0+pert and projects them down the normal direction until they lie within tol of the equal energy slice.

!subroutine proj_eqen(Y,Y0,params,params0,normal,tol)
!	implicit none
!
!	real(dp), dimension(:), intent(inout) :: Y, params
!	real(dp), dimension(:), intent(in) :: normal, Y0, params0
!	real(dp), dimension(size(normal)) :: norm
!	real(dp), intent(in) :: tol
!	real(dp), dimension(size(Y)-1+size(params)) :: vect, vect_norm
!	real(dp) :: alpha, V0, KE, en_left, en_left2, check, checknorm
!	integer :: i_1, i_2, i_3, i_4, i_5, i_6, i_7, counter
!	
!	norm = normal
!	counter = 0
!
!	!Load perturbed vector.
!	do i_1=1,4
!		vect(i_1)=Y(i_1+1)
!	end do
!	do i_2=1,4
!		vect(i_2+4)=params(i_2)
!	end do
!	!Load initial vector.
!	do i_6=1,4
!		vect_norm(i_6) = Y0(i_6+1)
!	end do
!	do i_7=1,4
!		vect_norm(i_7+4) = params0(i_7)
!	end do
!
!	!Check if on eq en slice.
!	en_left = en_diff(vect)
!print*,"Projection en left1 ", en_left
!	if(ABS(en_left)<=tol) RETURN
!
!	!Determine which way normal points and adjust to move onto eq en slice.
!	vect_norm = vect_norm + norm
!	checknorm = en_diff(vect_norm)
!print*,"Checknorm ",checknorm
!	if(checknorm*en_left<0)then
!		norm = (-1_dp)*norm
!	end if
!
!	!Step size.
!	alpha = en_left*.001_dp
!
!	!Move to eq en slice.
!	do i_5=1,500000
!		!Translate down norm by alpha.
!		vect = vect - alpha*norm
!		
!		en_left2 = en_diff(vect)
!
!!print*,"Projection en left 2", en_left2
!		
!		!Check for overshoot.
!		if(en_left2*en_left<0)then
!			en_left = en_left2
!			norm = -1_dp*norm
!			alpha = alpha*.001_dp
!print*,"Overshoot test"
!		end if
!
!		if(ABS(en_left2)<=tol) exit 
!		
!	end do
!
!	!Rebuild Y and params.
!	do i_3=1,4
!		Y(i_3+1)=vect(i_3)
!	end do
!	do i_4=1,4
!		params(i_4)=vect(i_4+4)
!	end do
!print*,"Y ",Y
!print*,"Params",params	
!
!end subroutine proj_eqen





!******************************************************************************************
!Subroutine which seeds the MCMC routine somewhere on the equal energy slice.

subroutine SEED_IC_EQEN(Y)
	implicit none

	real(dp), dimension(:), intent(inout) :: Y
	real(dp) ::  rand_1, rand_2, V, left

	do
	Y(1)=0_dp
	Y(4)=0_dp
	call random_number(rand_1)
	call random_number(rand_2)
	Y(2) = rand_2 *(phi_max-phi_min)+phi_min
	Y(3) = rand_1 *(psi_max-psi_min)+psi_min
	V = (lambda**4_dp)*((1_dp-((Y(3)*Y(3))/(m*m)))**2_dp +((Y(2)*Y(2))/(mu*mu)) + &
		&((Y(2)*Y(2)*Y(3)*Y(3))/ (nu**4_dp)))
	left = 2_dp*(energy_scale**4 - V)
	if(left<0) then
		cycle
	else
		Y(5)=SQRT(left)
		exit
	end if
	end do

end subroutine SEED_IC_EQEN


!*****************************************************************************************
!Function which will return the amount of difference in energy between a point in R^n and the surface E^4.

!real(dp) function en_diff(vect)
!	implicit none
!
!	real(dp), dimension(8) :: vect
!	real(dp) :: V0, KE
!
!	V0 = (vect(5)**4)*((1_dp-((vect(2)*vect(2))/(vect(6)*vect(6))))**2_dp +&
!		&((vect(1)*vect(1))/(vect(7)*vect(7))) + &
!		&((vect(1)*vect(1)*vect(2)*vect(2))/(vect(8)**4)))
!	KE = .5_dp*(vect(3)*vect(3)+vect(4)*vect(4))
!	en_diff = energy_scale**4 - KE -V0
!
!end function en_diff 

!Calculates the energy missing.  Send the function all of the given variables.
!The energy remaining for the variable which isn't sent to the function will be
!returned.
!pure real(dp) function en_diff(phi,psi,phiv,psiv)
!  	implicit none
!
!  real(dp), optional, intent(in) :: phi, psi, phiv, psiv
!  real(dp), dimension(5) :: y
!
!  if (.not. present(phiv)) then
!    en_diff=energy_scale**4-.5_dp*(psiv**2)-v_h(phi,psi)
!  else if (.not. present(psiv)) then
!    en_diff=energy_scale**4-.5_dp*(phiv**2)-v_h(phi,psi)
!  else if (.not. present(phi)) then
!    en_diff=energy_scale**4-.5_dp*(phiv**2+psiv**2)-v_h(0_dp,psi)
!  else if (.not. present(psi)) then
!    en_diff=energy_scale**4-.5_dp*(phiv**2+psiv**2)-v_h(phi,m)
!  end if
!
!end function en_diff








end module d_hybrid_initialconditions
