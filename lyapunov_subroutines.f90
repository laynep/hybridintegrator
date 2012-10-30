!Module that contains necessary subroutines for the calculation of Lyapunov
!exponents in conjunction with the program file lyapunov_integrator.f90.

!Subroutines for the lesnls/leslis software.

!NOTE: y(1)=n, y(2)=phi, y(3)=psi, y(4)=phi_dot, y(5)=psi_dot.

module lyapunov_subroutines
  use types, only : dp
  use numerical_routines, only : cspline
  implicit none

!  	real(dp), dimension(5) :: y_1, y_2, v_1, v_2
!  	real(dp) :: x_1, x_2

  !Do we want to include friction in the calculation? Var set in main.
  logical, public :: friction

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !subroutine to set the leslis parameter values.

  subroutine leslis_parameters(ipar,tolt,tolq,toll,work)
  	implicit none

  	integer, dimension(:), intent(out) :: ipar
  	real(dp), intent(out) :: tolt, tolq, toll(:)
  	integer, intent(in) :: work
	
  	ipar(1) = 0 !variable step size.
  	ipar(2) = 0 !integr from t0->te
  	ipar(3) = 1 !set 1 before init and before new call to lesnls
  	ipar(4) = 1 !yes ic on y
  	ipar(5) = work
    ipar(6) = 0 !1 integrate solution only, 0 approx LEs, too.
    ipar(7) = 0
  	ipar(8) = 0 !default method to calc lexps
  	ipar(9) = 0 !default method.
  	ipar(10) = 10 !error control on exps and q

  	tolt = 1e-6_dp	!tolerances
  	tolq = 1e-6_dp
  	toll = 1e-6_dp

  end subroutine leslis_parameters

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !reinitialize leslis parameter values for a different run.

  subroutine reinit_leslis(t0,te,ipar,le_dt)
  	implicit none

  	integer, dimension(:), intent(out) :: ipar
  	real(dp), intent(out) :: t0, te
  	real(dp), intent(in) :: le_dt

  	t0=0_dp
  	te=le_dt
  	ipar(3) = 1	!New call to LESLIS integrator.

  end subroutine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !subroutine to set the lesnls parameter values.

  pure subroutine lesnls_parameters(ipar,rearr,tolt,tolq,toll,work,x0,lyap)
    use d_hybrid_initialconditions
  	implicit none

  	integer, dimension(:), intent(out) :: ipar
  	real(dp), intent(out) :: tolt, tolq, toll(:), x0(:,:)
    real(dp), dimension(:), intent(out) :: rearr
    logical, intent(in) :: lyap
  	integer, intent(in) :: work
    integer :: i, j
	
  	ipar(1) = 0 !variable step size.
    if (lyap) then
      ipar(2) = 1 !integr one time step.
    else
  	  ipar(2) = 0 !integr from t0->te
    end if
  	ipar(3) = 1 !set 1 before init and before new call to lesnls
  	ipar(4) = 0 !no ic on x
  	ipar(5) = work
  	ipar(6) = 0 !calc lexps. 1 else.
    ipar(7) = 0
  	ipar(8) = 0 !default method to calc lexps
  	ipar(9) = 1 !default method.
  	ipar(10) = 20 !error control on traj and exps

    !Pass params to LESNLS subroutines.
  	rearr(1) = lambda
  	rearr(2) = m
  	rearr(3) = mu
  	rearr(4) = nu
  	rearr(5) = beta

    !Tolerances
  	tolt = 1e-9_dp
  	tolq = 1e-9_dp
  	toll = 1e-9_dp

    !IC on fundamental matrix.
    do i=1,size(x0,1)
      do j=1,size(x0,2)
        if(i==j) then
          x0=1e0_dp
        else
          x0=0e0_dp
        end if
      end do
    end do

  end subroutine lesnls_parameters

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !reinitialize lesnls parameter values for a different run.

  pure subroutine reinit_lesnls(t0,te,dt,ipar)
  	implicit none

  	integer, dimension(:), intent(out) :: ipar
  	real(dp), intent(inout) :: t0, te, dt

  	t0=0e0_dp
  	te=t0+dt
  	ipar(3) = 1

  end subroutine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine le_stats(ipar, iflag, t0)
  	implicit none

  	integer, intent(in) :: ipar(*), iflag
  	real(dp), intent(in) :: t0

        if (ipar(8).eq.0) then
           print *, 'method used for les is cont qr with proj-dp5 '
          else if (ipar(8).eq.1) then
           print *, 'method used for les is cont qr with hybrid-dp5 '
          else if (ipar(8).eq.2) then
           print *, 'method used for les is cont qr with proj-rk38 '
          else if (ipar(8).eq.3) then
           print *, 'method used for les is cont qr with hybrid-rk38 '
          else if (ipar(8).eq.4) then
           print *, 'method used for les is disc qr with dp5 '
          else if (ipar(8).eq.5) then
           print *, 'method used for les is disc qr with rk38 '
        endif
        if (ipar(8).le.3) then
           if (ipar(9).eq.0) then
              print *, 'les by nu-integration '
             else if (ipar(9).eq.1) then
              print *, 'les by comp-trap-rule '
           endif
        endif
        if (ipar(10).eq.0) then
           print *, 'error control on trajectory '
          else if (ipar(10).eq.1) then
           print *, 'error control on les '
          else if (ipar(10).eq.2) then
           print *, 'error control on q '
          else if (ipar(10).eq.10) then
           print *, 'error control on trajectory and les '
          else if (ipar(10).eq.20) then
           print *, 'error control on trajectory and q '
          else if (ipar(10).eq.21) then
           print *, 'error control on les and q'
          else if (ipar(10).eq.210) then
           print *, 'error control on trajectory, les and q'
        end if
        print *, 'iflag = ', iflag, ' final t = ', t0
        print *, 'steps rejected = ', ipar(12)
        print *, 'number of steps = ', ipar(13)

  end subroutine le_stats

  !Subroutine to record time, Y vector, and lyapunov exponents in a linked list.
  subroutine rec_LE(lyap_exp,t,Y,le_list)
    use linked_list
    implicit none

    type(linkedlist), intent(inout) :: le_list
    real(dp), dimension(:), intent(in) :: lyap_exp
    real(dp), dimension(5), intent(in) :: y
    real(dp), intent(in) :: t
    type(llnode), pointer :: new
    integer :: i
    real(dp), dimension(size(lyap_exp)+6) :: arr

    !Make an array.
    arr(1)=t
    arr(2)=y(1)
    arr(3)=y(3)
    arr(4)=y(2)
    arr(5)=y(5)
    arr(6)=y(4)
    do i=1, size(lyap_exp)
      arr(i+6) = lyap_exp(i)
    end do

	  !Make a node out of the array.
	  call ll_make(new,arr)

	  !Append new node to list.
	  call ll_append(new,le_list)

  end subroutine rec_LE



end module lyapunov_subroutines

