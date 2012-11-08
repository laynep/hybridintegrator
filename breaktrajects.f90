!A program that will read in a file and break it into smaller files
!based on when the values in the first column start to decrease.  This
!is useful for breaking the trajectory files into smaller units.  In
!these files I have multiple trajectories being recorded with
!monotonically increasing e-folding stored in the first column.
!Consequently, I am able to break up the file containing many
!trajectories into many files containing only one trajectory.

!Only works for text files.  Does not work if the first column in the
!first row is less than zero.

program breaktrajects
  use types, only : dp
  use features, only : newunit
  implicit none

  real(dp), dimension(5) :: now
  real(dp) :: prev
  integer :: i, j, u, w, counter
	character(len=32) :: arg1
  character(len=13) :: newfilename

  !Trajectory counter.
  counter = 1

	!Get the bin file name from command line
	call get_command_argument(1,arg1)
  open(unit=newunit(u), file=arg1)

  !Make file name and see if it exists.
  call makenewfilename(newfilename,counter,w)

  i=1
  prev=0e0_dp
  do
    read(u,fmt=*) (now(j),j=1,size(now))
    if (now(1)<prev) then
      !Open new file.
      close(w)
      call makenewfilename(newfilename,counter,w)
      write(w,*) (now(j),j=1,size(now))
    else
      write(w,*) (now(j),j=1,size(now))
    end if
    prev = now(1)
  end do

end program breaktrajects

subroutine makenewfilename(newfilename,counter,newfilenumb)
  use features, only : newunit
  implicit none

  logical :: file_exists
  character(len=13) :: newfilename
  integer :: counter, newfilenumb

  do
    write(newfilename,'(a,i4.4,a)')'traj_',counter,".txt"
    inquire(file=newfilename, exist=file_exists)
    if (.not. file_exists) then
      open(unit=newunit(newfilenumb), file=newfilename, status="new")
      exit
    end if
    counter=counter+1
  end do
end subroutine makenewfilename
