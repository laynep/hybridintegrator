!#############################################################################
!Subroutines for building an MCMC Metropolis sampling of initial conditions from
!a given data set for hybrid inflation.

!USE: The ic_metr routine will build a nearest neighbor sample from a data set
!stored in sample_table.  It will then perform a Metropolis sampling of that
!set, where the area over which the sample's interpolation has been done is
!specified in the variable epsil.  The first time this routine would be called,
!instead call ic_metr_init so that all relevant variables can be declared, etc.
!The proposal distribution is a semi-flat sampler over the allowed range of
!initial conditions (d_ic_eqen).  New points are accepted from old points at the
!rate min(1,p_new/p_old) where p_new is the number of points in the
!interpolation region for the new point and p_old is the same for the old point.
!By option we can turn off the ability of the routine to take the same point for
!subsequent links in the Markov chain.  Turning this off makes the sampler
!pseudo-Metropolis, giving a non-representative sample of the interpolation
!function, but it allows us to get unique points to test off of.
!#############################################################################

module hybrid_metropolis
  use d_hybrid_initialconditions, only : numbpoints_sample, energy_scale,&
    &m_planck, xxglobal_fail, xxglobal_succ, d_ic_eqen
  use types, only : dp
  implicit none

contains

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

  !Load the sample table from file where file data stored in namelist in
  !parameter file.  Make sure to switch the data so stored as
  !(phi,psi,phi_dot,psi_dot)
  switch=.true.
  call load_sample_table(sample_table,switch)

  !Set first IC as first element in sample_table.
  y(:)=sample_table(1,:)
	!Set IC at random on EQEN slice.
	!call d_ic_eqen(y,iccounter)

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
    doi: do i=1,n
      doj: do j=start,size(boxcover,1)
        dok1:	do k=1,size(boxcover,2)
			  	check=.true.
				 	if (boxcover(j,k).ne.boxcover(j-1,k)) then
				   	check=.false.
				   	exit dok1
				 	end if
			  end do dok1
				if (check) then
					numbpoints_sample(i,width)=numbpoints_sample(i,width)+1
				else
          dok2:	do k=1,size(boxcover,2)
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

      doi: do i=start,ending
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
    			get_prob=numbpoints_sample(i,size(numbpoints_sample,2))*1e0_dp
          exit doi
    		end if
      end do doi

    end function get_prob

end function accept_ratio



end module hybrid_metropolis
