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

	phi_0 = 0e0_dp
	psi_0 = 4e0_dp*m_planck
	phi_dot_0 =0e0_dp
	psi_dot_0 =0e0_dp
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

    if (rho_kinetic<0) then
			cycle
		end if

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

!#############################################################################
!Metropolis Subroutines.

!Subroutine that will take as input a point y0 and a sample table, sample_table.  This subroutine will perform a nearest neighbor interpolation on sample_table with respect to the density of points, within a given box size, eps.  By default this will not give duplicates, i.e. y_init != y_fin.  If you want to override this, specify dup=1.  The interpolated density function is stored in numbpoints_sample.

subroutine ic_metr(y,sample_table,iccounter,epsil,dup)
	use sorters, only : locate, heapsort
	implicit none

	real(dp), dimension(:), intent(inout) :: y
	real(dp), dimension(:,:), allocatable, intent(inout) :: sample_table
	real(dp), optional, intent(inout) :: epsil
	integer, intent(in) :: iccounter
	integer :: test
	integer, optional, intent(in) :: dup
	real(dp) :: accept, rand, eps
	real(dp), dimension(5) :: yprop
	integer :: i,j,k

	!Load eps if not provided.  Default to value of Hubble parameter.
	if (present(epsil)) then
    eps=epsil
  else
    eps=energy_scale**2/m_planck
  end if

	!Load the density calculation if not already done.
  !numbpoints_sample global var defined in this mod.
	if (.not. allocated(numbpoints_sample)) then
    call build_nearest_neighbor(sample_table, eps)
  end if

  !Metropolis sampling.
	do
		!Get a new point on eqen slice.
		call d_ic_eqen(yprop,iccounter)

		!Calc accept ratio
		accept=accept_ratio(y,yprop,eps)

		!Gen rand numb
		call random_number(rand)
	
		!Move to new point if rand<a
		if (rand<accept) then
			xxglobal_succ=xxglobal_succ+1
			y=yprop
			exit
		end if
		xxglobal_fail=xxglobal_fail+1

		!If no duplicates, then cycle until get unique numb, otherwise return
    !dup is present iff we will accept duplicates, i.e. a proper MCMC alg.
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
	integer :: i, j, iend
  logical :: switch

	!Set IC at random on EQEN slice.
	!call d_ic_eqen(y,iccounter)

  !Load the sample table from file where file data stored in namelist in
  !parameter file.  Make sure to switch the data so stored as
  !(phi,psi,phi_dot,psi_dot)
  switch=.true.
  call load_sample_table(sample_table,switch)

  !Set first IC as first element in sample_table.
  y(:)=sample_table(1,:)

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
			call ic_metr(y,sample_table,iccounter,eps)
		else
			call ic_metr(y,sample_table,iccounter)
		end if
	end do

end subroutine ic_metr_init

!Subroutine that builds the nearest neighbor interpolation of the sample_table.
!Loads into the global variable numbpoints_sample.
subroutine build_nearest_neighbor(sample_table, eps)
  use sorters, only : heapsort
  implicit none

	integer, dimension(:,:), allocatable :: boxcover
  real(dp), intent(in) :: eps
	real(dp), dimension(:,:), allocatable, intent(inout) :: sample_table
  integer :: i, j, k, n, start, width
  logical :: check

		!Box cover, of size eps.
		allocate(boxcover(size(sample_table,1),size(sample_table,2)))
    do i=1, size(boxcover,1)
		  boxcover(i,:)=ceiling(sample_table(i,:)/eps)
    end do
		if (allocated(sample_table)) deallocate(sample_table)

		!Count numb of boxes it takes to cover sample_table.
		call heapsort(boxcover)
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

		!Each row in numbpoints_sample will be the n-D position of the
    !upper-right-hand-corner of the box expressed in the first n entries of the
    !row and the (n+1)st entry will be the number of points in that box.
    width=size(boxcover,2)+1
		allocate(numbpoints_sample(n,width))
    !Initialize.
    numbpoints_sample=0
    do i=1,size(boxcover,2)
		  numbpoints_sample(1,i)=boxcover(1,i)
    end do
		numbpoints_sample(1,size(numbpoints_sample,2))=1

		!Count elts per box. Copy to final entry in corresponding row of
    !numbpoints_sample.
		start=2
doi:		do i=1,n
doj:  		do j=start,size(boxcover,1)
dok1:   		do k=1,size(boxcover,2)
			    		check=.true.
				    	if (boxcover(j,k).ne.boxcover(j-1,k)) then
					    	check=.false.
					    	exit dok1
				    	end if
			    	end do dok1
				if (check) then
					numbpoints_sample(i,width)=numbpoints_sample(i,width)+1
				else
dok2:			do k=1,size(boxcover,2)
						numbpoints_sample(i+1,k)=boxcover(j,k)
					end do dok2
					numbpoints_sample(i+1,width)=1
					start=j+1
					exit doj
				end if
			end do doj
		end do doi
		!Deallocate boxcover & sample_table
		if (allocated(boxcover)) deallocate(boxcover)
		if (allocated(sample_table)) deallocate(sample_table)

end subroutine build_nearest_neighbor

!Subroutine to load the sample table from file.  File data stored in the
!namelist given in parameters_hybrid.txt as &sample.  The sample_table will be
!loaded in as (0,data1,data2,...,dataN) for N-dimensional data.
!NOTE: Optional argument switch is invoked when we want to switch a file with
!vectors stored in the form (psi,phi,psi_dot,phi_dot) to vectors we will
!manipulate stored as (phi,psi,phi_dot,psi_dot).
subroutine load_sample_table(sample_table,switch)
  use features, only : newunit
  implicit none

  real(dp), dimension(:,:), allocatable, intent(out) :: sample_table
  real(dp), dimension(:), allocatable :: interm_vect
  integer :: u, i, j, ierr
  integer :: samp_len, samp_wid
  character(len=100) :: datafile
  logical, optional, intent(in) :: switch

  namelist /sample/ samp_len, samp_wid, datafile

  !Get info on sample table from namelist.
 	open(unit=newunit(u), file="parameters_hybrid.txt", status="old",&
 	& delim = "apostrophe")
 	read(unit=u, nml=sample)
 	close(unit=u)
 	datafile=trim(datafile)

  allocate(sample_table(samp_len,samp_wid+1))
 	!Data *must* be in form y(2),...,y(n), where y(1)=0_dp assumed.
  sample_table=0e0_dp
 	open(unit=newunit(u), file=datafile, status="old", form="unformatted")
 	do i=1,samp_len
  	read(unit=u,iostat=ierr) (sample_table(i,j),j=2,samp_wid+1)
  	if (is_iostat_end(ierr)) exit
 	end do
 	close(unit=u)

  !If switch is both present and set to .true. then swap the columns 2&3 and
  !columns 4&5.  Requires an intermediary vector.
  allocate(interm_vect(size(sample_table,1)))
  if (present(switch) .and. switch) then
    if (samp_wid<4) then
      print*,"ERROR: can't perform sample_table switch. Not enough columns."
      stop
    end if
    interm_vect(:)=sample_table(:,2)
    sample_table(:,2)=sample_table(:,3)
    sample_table(:,3)=interm_vect(:)
    interm_vect(:)=sample_table(:,4)
    sample_table(:,4)=sample_table(:,5)
    sample_table(:,5)=interm_vect(:)
  end if
  deallocate(interm_vect)

end subroutine load_sample_table


!function to calculate the acceptance ratio for the Metropolis algorithm for the density function obtained in IC_METR.  y1 is old point, y2 is new point.

real(dp) function accept_ratio(y1,y2,eps)
	use sorters, only : locate
	implicit none

	real(dp), dimension(:), intent(in) :: y1, y2
	real(dp), intent(in) :: eps
	real(dp) :: a, p1,p2

	!Calc probabilities.
  p1= get_prob(y1, eps)
  p2= get_prob(y2, eps)

!print*,"y1",y1
!print*,"y2",y2
!print*,"from accept",p2,p1

  if (p2<0e0_dp .or. p1<0e0_dp) then
    print*,"ERROR: Negative probability in accept ratio."
	!Avoid div by almost-zero
	else if (p1<=1e-10_dp) then
		a=1e0_dp
	else
		a=p2/p1
	end if

  !MIN(1e0_dp,a)
  if (a<1e0_dp) then
    accept_ratio=a
  else
    accept_ratio=1e0_dp
  end if

  contains

    !Find where to start in table. We use locate in the second column of
    !numbpoints_sample, since the first column should be all 0 bc it's the
    !e-folds.
    subroutine find(start,ending,y)
      implicit none

      integer, intent(out) :: start, ending
      real(dp), intent(in) :: y(:)
      integer :: width, i, low, high

      width=size(numbpoints_sample,2)
      low =ceiling(y(2)/eps)
      high =ceiling(y(2)/eps)+1

      call locate(numbpoints_sample(:,2:width),low,start)
	    call locate(numbpoints_sample(:,2:width),high,ending)

    end subroutine find

    !Compares box for vect y vs rows in numbpoints_sample to see if they're the same.  If
    !so, then take the number of points in the box as value for p1, else assume
    !that there are no points in that box and leave p1=0.
    real(dp) function get_prob(y, eps)
      implicit none

      integer :: start, ending, i, j
      real(dp), intent(in) :: y(:), eps
      real(dp), dimension(size(y)) :: test
      logical :: check

      !Find start and ending search points.
      call find(start,ending,y)

	    get_prob=0e0_dp

doi:  do i=start,ending
        !Load comparison vector with first N elements in ith row of
        !array numbpoints_sample
        test(:)=numbpoints_sample(i,1:size(y))

        !Check if y is same as test.
  doj1:	do j=1,size(y)
    			check=.true.
    			if (ceiling(y(j)/eps).ne.test(j)) then
    				check=.false.
    				exit doj1
    			end if
    		end do doj1

    		if (check) then
!print*,"check is same",y
!print*,numbpoints_sample(i,2:5)

    			get_prob=numbpoints_sample(i,size(numbpoints_sample,2))*1e0_dp
          exit doi
    		end if
      end do doi

    end function get_prob

end function accept_ratio


end module d_hybrid_initialconditions

