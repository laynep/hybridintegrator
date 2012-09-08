!**************************************************************************************
!Layne Price, University of Auckland, May 5, 2012

!A series of subroutines which set the initial conditions for hybrid inflation subject to a set of constraints, dependent on each subroutine.

!Note that the parameters_hybrid() subroutine must be declared first, before the other routines can be used.

!SUBROUTINE parameters_hybrid() ................... Sets the particle physics 
!						parameters for hybrid inflation.
!FUNCTION V_h(Y)................................... Gives the hybrid potential.
!SUBROUTINE D_IC_EQEN(Y,iccounter) ................ Sets the ICs according to equal
!						energy slice. Field values are set
!						randomly and field velocities are 
!						alternately set on subsequent calls
!						by energy constraint.
!SUBROUTINE D_IC_ZEROV(Y) ......................... Sets IC with zero velocity.  Clesse
!						original.

!Y(1)=N, Y(2)=phi, Y(3)=psi, Y(4)=phi_dot, Y(5)=psi_dot.

!**************************************************************************************

MODULE d_hybrid_initialconditions
IMPLICIT NONE
 
	!Global data.
	DOUBLE PRECISION :: phi_0, psi_0, phi_dot_0, psi_dot_0
	DOUBLE PRECISION, PARAMETER :: pi= 3.14159265358979323846D0
	DOUBLE PRECISION, PARAMETER :: m_planck=1000D0
	DOUBLE PRECISION, PARAMETER :: beta=8.37758040957278D0/(m_planck*m_planck)
	DOUBLE PRECISION :: m
	DOUBLE PRECISION :: mu
	DOUBLE PRECISION, PARAMETER :: phi_min=0D0
	DOUBLE PRECISION, PARAMETER :: phi_max=.2D0*m_planck
	DOUBLE PRECISION, PARAMETER :: psi_min=0D0
	DOUBLE PRECISION, PARAMETER :: psi_max=.2D0*m_planck
	DOUBLE PRECISION :: nu
	DOUBLE PRECISION :: energy_scale
	DOUBLE PRECISION :: lambda

	!For the nearest neighbor sampling.
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: numbpoints_sample
	INTEGER :: xxglobal_fail, xxglobal_succ

	!From Clesse.
	DOUBLE PRECISION, PARAMETER :: m_min = .049787068D0*m_planck
	DOUBLE PRECISION, PARAMETER :: m_max = .932239382D0*m_planck
	DOUBLE PRECISION, PARAMETER :: mu_min = 0.367879441D0*m_planck
	DOUBLE PRECISION, PARAMETER :: mu_max = 54.598150033D0*m_planck
	DOUBLE PRECISION, PARAMETER :: nu_min = 0.22313016*m_planck
	DOUBLE PRECISION, PARAMETER :: nu_max = 2.718281828*m_planck
	DOUBLE PRECISION, PARAMETER :: lambda_min = -2D0*m_planck
	DOUBLE PRECISION, PARAMETER :: lambda_max = 2D0*m_planck

	NAMELIST /parameters/ m, mu, nu, lambda, energy_scale

	

CONTAINS

!*************************************************************************************
!Subroutine which declares potential parameters for hybrid inflation.

SUBROUTINE parameters_hybrid()
IMPLICIT NONE

	OPEN(unit=10000, file="parameters_hybrid.txt", status="old", delim = "apostrophe")
	READ(unit=10000, nml=parameters)
	CLOSE(unit=10000)

END SUBROUTINE parameters_hybrid

!*************************************************************************************
!FUNCTION which describes the potential.

pure DOUBLE PRECISION FUNCTION V_h(Y)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:), intent(in) :: Y

	V_h = (lambda*lambda*lambda*lambda)*((1D0 - ((Y(3)*Y(3))/(m*m)))**2D0 +&
		& ((Y(2)*Y(2))/(mu*mu)) + ((Y(2)*Y(2)*Y(3)*Y(3))/ (nu*nu*nu*nu)))

END FUNCTION V_h


!******************************************************************************************


!This subroutine will set the initial conditions on an equal energy slice at the start of inflation that corresponds to the value energy_scale which is specified in the parameter declaration in the main program.

SUBROUTINE D_IC_EQEN(Y,iccounter)
IMPLICIT NONE
	DOUBLE PRECISION, INTENT(OUT) :: Y(5)
	INTEGER, INTENT(IN) :: iccounter
	DOUBLE PRECISION :: rand_1, rand_2, rand_3, rand_4,rand_5, V_0, rho_kinetic, phi_dot_min, phi_dot_max, psi_dot_min, psi_dot_max
	INTEGER :: param_constr

	!Gives parameter to set by energy constraint: 0~phi_dot,1~psi_dot.  The IC for iccounter equals number of parameters to oscillate between.
	param_constr = MOD(iccounter,2)

	IF (param_constr==0) THEN
		DO
			!Set phi from energy constraint.
			CALL random_number(rand_1)
			Y(2)= (rand_1*(phi_max-phi_min)) + phi_min
			phi_0 = Y(2)
	
			!IC for Y(3)~psi randomly in range psi_min to psi_max.
			CALL random_number(rand_2)
			Y(3)= (rand_2*(psi_max-psi_min)) + psi_min
			psi_0 = Y(3)
	
			!Initial value of the potential.
			V_0 = V_h(Y)
	
			!Energy density remaining in kinetic term.
			rho_kinetic = (energy_scale**4D0) - V_0

			IF (rho_kinetic<0) CYCLE 
			
			!Set so that psi_dot is chosen with flat prior.
			psi_dot_max = SQRT(2D0*rho_kinetic)
			psi_dot_min = -1D0*SQRT(2D0*rho_kinetic)
	
			!Sets the psi_dot IC to the range psi_dot_max to psi_dot_min
			CALL random_number(rand_3)
			Y(5) = (rand_3*(psi_dot_max-psi_dot_min)) + psi_dot_min
			psi_dot_0 = Y(5)
	
			IF(2D0*rho_kinetic - (Y(5)*Y(5)) > 0) EXIT
	
		END DO
	
		!Set the phi_dot IC by the total energy density constraint.
		CALL random_number(rand_4)
		IF(rand_4 < .5) THEN
			Y(4) = SQRT(2D0*rho_kinetic - (Y(5)*Y(5)))
			phi_dot_0 = Y(4)
		ELSE 
			Y(4) = -1D0*SQRT(2E0*rho_kinetic - (Y(5)*Y(5)))
			phi_dot_0 = Y(4)
		END IF
		
	ELSEIF (param_constr==1) THEN
		DO
			!IC for Y(2)~phi randomly in range phi_min to phi_max.
			CALL random_number(rand_1)
			Y(2)= (rand_1*(phi_max-phi_min)) + phi_min
			phi_0 = Y(2)
	
			!IC for Y(3)~psi randomly in  range psi_min to psi_max.
			CALL random_number(rand_2)
			Y(3)= (rand_2*(psi_max-psi_min)) + psi_min
			psi_0 = Y(3)
	
			!Initial value of the potential.
			V_0 = V_h(Y)
	
			!Energy density remaining in kinetic term.
			rho_kinetic = (energy_scale**4D0) - V_0

			IF (rho_kinetic<0) CYCLE 
	
			!Set so that phi_dot is chosen with flat prior.
			phi_dot_max = SQRT(2D0*rho_kinetic)
			phi_dot_min = -1D0*SQRT(2D0*rho_kinetic)
	
			!Sets the phi_dot IC to the range phi_dot_max to phi_dot_min
			CALL random_number(rand_3)
			Y(4) = (rand_3*(phi_dot_max-phi_dot_min)) + phi_dot_min
			phi_dot_0 = Y(4)
	
		IF(2D0*rho_kinetic - (Y(4)*Y(4)) > 0) EXIT
	
		END DO
	
		CALL random_number(rand_5)
		IF(rand_5<.5)THEN
			!Set the psi_dot IC by the total energy density constraint.
			Y(5) = SQRT(2D0*rho_kinetic - (Y(4)*Y(4)))
			psi_dot_0 = Y(5)
		ELSE
			!Set the psi_dot IC by the total energy density constraint.
			Y(5) = -SQRT(2D0*rho_kinetic - (Y(4)*Y(4)))
			psi_dot_0 = Y(5)
		END IF
	END IF

	!Reinitialize the e-fold value: Y(1)~N.
	Y(1)=0D0 
		

END SUBROUTINE D_IC_EQEN

!*********************************************************************************


!Program subroutine: This subroutine will follow the results of Clesse and set the initial values of the velocity parameters to zero and choose the initial conditions of the fields at random.

SUBROUTINE D_IC_ZEROV(Y)
IMPLICIT NONE
	DOUBLE PRECISION, INTENT(OUT) :: Y(5)
	DOUBLE PRECISION :: rand_1, rand_2

	!Set IC for Y(2)~phi randomly in range phi_min to phi_max.
	CALL random_number(rand_1)
	Y(2)= (rand_1*(phi_max-phi_min)) + phi_min
	phi_0 = Y(2)

	!Set IC for Y(3)~psi randomly in range psi_min to psi_max.
	CALL random_number(rand_2)
	Y(3)= (rand_2*(psi_max-psi_min)) + psi_min
	psi_0 = Y(3)

	!Sets velocities to zero initially.
	Y(4)=0D0
	phi_dot_0=Y(4)
	Y(5)=0D0
	psi_dot_0=Y(5)

	!Reinitialize the cosmology equation values.
	Y(1)=0D0 
		
	
END SUBROUTINE D_IC_ZEROV

!*******************************************************************************************
!Fixed initial conditions set here.

SUBROUTINE FIXED_IC(Y)

	DOUBLE PRECISION, DIMENSION(5), INTENT(OUT) :: Y

	phi_0 = 80.837797631740003D0
	psi_0 =.1D-6
	phi_dot_0 =0D0
	psi_dot_0 =0D0
	Y(1)=0D0
	Y(2)=phi_0
	Y(3)=psi_0
	Y(4)=phi_dot_0
	Y(5)=psi_dot_0

END SUBROUTINE FIXED_IC

!****************************************************************************
!Subroutine that slices the eqen surface.  This particular choice sets psi_0=0 so that all the initial conditions are set in (if the velocities were then zero) the inflationary valley.  psi_dot_0 is set by the energy constraint so that we can plot phi_0 vs phi_dot_0 as a two-dimensional representation of the three-dimensional slice.

SUBROUTINE EQEN_SLICING(Y)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: Y
	DOUBLE PRECISION :: rand_1, rand_2, rand_3, chi
	DOUBLE PRECISION :: phi_dot_max, phi_dot_min, rho_kinetic, psi_dot_max, psi_dot_min

	!Initialize.
	Y=0D0

	!No initial e-foldings.
	Y(1)=0D0

!	DO
!		!Set phi_0 randomly.
!		CALL random_number(rand_1)
!		Y(2)=(rand_1*(phi_max-phi_min)) + phi_min
!		phi_0 = Y(2)

		!Constrain psi_0 so in inflationary valley.
!		Y(3)=0D0
!		psi_0 = Y(3)

		!Choose phi_dot_0 uniformly from remaining energy.
!		rho_kinetic = (energy_scale**4D0) - V_h(Y)

!		IF (rho_kinetic<0) CYCLE

!		phi_dot_max = SQRT(2D0*rho_kinetic)
 !       	phi_dot_min = -1D0*SQRT(2D0*rho_kinetic)
!		CALL random_number(rand_2)
!		Y(4)=(rand_3*(phi_dot_max-phi_dot_min)) + phi_dot_min
!		phi_dot_0= Y(4)

		!Set psi_dot_0 by energy constraint
!		psi_dot_max = SQRT(2D0*(energy_scale**4D0 - V_h(Y) - .5D0*Y(4)*Y(4)))
!		psi_dot_min = 0D0 	! -1D0*psi_dot_max
!		CALL random_number(rand_3)
!		Y(5)=(rand_3*(psi_dot_max-psi_dot_min)) + psi_dot_min
!		psi_dot_0=Y(5)

!		EXIT
!	END DO


	!Initialize
        Y=0D0

	DO
          	!Set phi_0 randomly.
                CALL random_number(rand_1)
                Y(2)=(rand_1*(phi_max-phi_min)) + phi_min
                phi_0 = Y(2)

                !Set psi_0 randomly.
                CALL random_number(rand_1)
                Y(3)=(rand_1*(psi_max-psi_min)) + psi_min
                psi_0 = Y(3)

                !Find remaining energy.
                rho_kinetic = (energy_scale**4D0) - V_h(Y)

                IF (rho_kinetic<0) THEN
			CYCLE
		END IF

                chi = -1D0*SQRT(rho_kinetic)
                Y(4)=chi
                phi_dot_0=chi
                Y(5)=chi
                psi_dot_0=chi

                EXIT

        END DO


END SUBROUTINE EQEN_SLICING

SUBROUTINE EQEN_EQVEL(Y)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: Y
	DOUBLE PRECISION :: rand_1, rho_kinetic, chi
	
	!Initialize
	Y=0D0

	DO
		!Set phi_0 randomly.
                CALL random_number(rand_1)
                Y(2)=(rand_1*(phi_max-phi_min)) + phi_min
                phi_0 = Y(2)

		!Set psi_0 randomly.
                CALL random_number(rand_1)
                Y(3)=(rand_1*(psi_max-psi_min)) + psi_min
                psi_0 = Y(3)

		!Find remaining energy.
                rho_kinetic = (energy_scale**4D0) - V_h(Y)

		IF (rho_kinetic<0) CYCLE

		chi = SQRT(rho_kinetic)
		Y(4)=chi
		phi_dot_0=chi
		Y(5)=chi
		psi_dot_0=chi

		EXIT

	END DO


END SUBROUTINE EQEN_EQVEL



!SUBROUTINE that will take as input a point Y0 and a sample table, sample_table.  This subroutine will perform a nearest neighbor interpolation on sample_table with respect to the density of points, within a given box size, eps.  By default this will not give duplicates, i.e. Y_init != Y_fin.  If you want to override this, specify dup=1 .

SUBROUTINE IC_METR(Y,sample_table,iccounter,eps,dup,test)
USE sorters, ONLY : locate, heapsorttotal
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: Y
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: sample_table
	DOUBLE PRECISION, OPTIONAL, INTENT(INOUT) :: eps
	INTEGER, INTENT(IN) :: iccounter
	integer, optional, intent(in) :: test
	INTEGER, OPTIONAL, INTENT(IN) :: dup
	DOUBLE PRECISION :: accept, rand
	DOUBLE PRECISION, DIMENSION(5) :: YPROP
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: boxcover
	INTEGER :: i,j,k, n, start
	LOGICAL :: check

	!Load eps if not provided.
	IF (PRESENT(eps) .EQV. .FALSE.) eps=1D0

	!Load the density calculation if not already done.
	IF (ALLOCATED(numbpoints_sample) .EQV. .FALSE.) THEN
		!Box cover, of size eps.
		ALLOCATE(boxcover(SIZE(sample_table,1),SIZE(sample_table,2)))
		boxcover=CEILING(sample_table/eps)
		IF (ALLOCATED(sample_table)) DEALLOCATE(sample_table)
		!Count numb of boxes.
		CALL heapsorttotal(boxcover)
		n=1
		DO i=2,SIZE(boxcover,1)
			DO j=1,SIZE(boxcover,2)
				check=.TRUE.
				IF (boxcover(i,j).NE.boxcover(i-1,j)) THEN
					check=.FALSE.
					EXIT
				END IF
			END DO
			IF (check) CYCLE
			n=n+1			
		END DO
		!Allocate numbpoints_sample.  Each row will be Y(2),...,Y(5),# of Boxes.
		ALLOCATE(numbpoints_sample(n,5))
		!Load first elt of numbpoints_sample
		DO i=1,SIZE(boxcover,2)
			numbpoints_sample(1,i)=boxcover(1,i)*eps
		END DO
		numbpoints_sample(1,5)=1
		!Count elts per box. Copy to numbpoints_sample.
		start=2
doi:		DO i=1,n
	doj:		DO j=start,SIZE(boxcover,1)
		dok1:		DO k=1,SIZE(boxcover,2)
					check=.TRUE.
					IF (boxcover(j,k).NE.boxcover(j-1,k)) THEN
						check=.FALSE.
						EXIT dok1
					END IF
				END DO dok1
				IF (check) THEN
					numbpoints_sample(i,5)=numbpoints_sample(i,5)+1
				ELSE
		dok2:			DO k=1,SIZE(boxcover,2)
						numbpoints_sample(i+1,k)=boxcover(j,k)
					END DO dok2
					numbpoints_sample(i+1,5)=1
					start=j+1
					EXIT doj
				END IF
			END DO doj
		END DO doi
		
		!Deallocate boxcover & sample_table
		IF (ALLOCATED(boxcover)) DEALLOCATE(boxcover)
		IF (ALLOCATED(sample_table)) DEALLOCATE(sample_table)
	END IF
	DO
		!GET A NEW POINT
		if (.not. present(test)) then 
			call D_IC_EQEN(YPROP,iccounter)
		else
			!TEST!!!!
			call EQEN_SLICING(YPROP)
		end if

		!CALC ACCEPT RATIO
		accept=accept_ratio(Y,YPROP,eps)

		!GEN RAND NUMB
		CALL random_number(rand)
	
		!MOVE TO NEW POINT IF RAND<A
		IF (rand<accept) THEN
			xxglobal_succ=xxglobal_succ+1
			Y=YPROP
			EXIT
		END IF
		xxglobal_fail=xxglobal_fail+1

		!IF NO DUPLICATES, THEN CYCLE UNTIL GET UNIQUE NUMB, OTHERWISE RETURN
		IF (PRESENT(dup)) EXIT		
	END DO

END SUBROUTINE IC_METR



!FUNCTION to calculate the acceptance ratio for the Metropolis algorithm for the density function obtained in IC_METR.  Y1 is old point, Y2 is new point.

DOUBLE PRECISION FUNCTION accept_ratio(Y1,Y2,eps)
USE sorters, ONLY : locate
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: Y1, Y2
	INTEGER :: start, ending, i, j
	DOUBLE PRECISION, INTENT(IN) :: eps
	DOUBLE PRECISION :: a, p1,p2
	INTEGER, DIMENSION(5) :: test
	LOGICAL :: check

	!Calc probabilities.
	!Find where to start in table. Returns first value below Y1(2).
	CALL locate(numbpoints_sample,Y1(2)/eps,start)
	CALL locate(numbpoints_sample,Y1(2)/eps+1D0,ending)
	p1=0D0

	DO i=start,ending
		test=(/0,numbpoints_sample(i,1),numbpoints_sample(i,2),&
			&numbpoints_sample(i,3),numbpoints_sample(i,4) /)

	doj1:	DO j=1,SIZE(Y1)
			check=.TRUE.
			IF (CEILING(Y1(j)/eps).ne.test(j)) THEN
				check=.FALSE.
				EXIT doj1
			END IF
		END DO doj1
		IF (check) THEN
			p1=numbpoints_sample(i,5)*1D0
			EXIT
		END IF
	END DO
	CALL locate(numbpoints_sample,Y2(2)/eps,start)
	CALL locate(numbpoints_sample,Y2(2)/eps+1D0,ending)
	p2=0D0
	DO i=start,ending
		test=(/0, numbpoints_sample(i,1),numbpoints_sample(i,2),&
			&numbpoints_sample(i,3),numbpoints_sample(i,4) /)
	doj2:	DO j=1,SIZE(Y2)
			check=.TRUE.
			IF (CEILING(Y2(j)/eps) .NE. test(j)) THEN
				check=.FALSE.
				EXIT doj2
			END IF

		END DO doj2

		IF (check) THEN
			p2=numbpoints_sample(i,5)*1D0
			EXIT
		END IF
	END DO

	IF (p2<=1.1D0) THEN
		a=0D0
	!Avoid div by zero
	ELSE IF (p1<=0D-10) THEN
		a=1D0
	ELSE
		a=p2/p1
	END IF
	accept_ratio=MIN(1D0,a)
	

END FUNCTION accept_ratio



!Gives uniform sample of IC from 0-1.  GOOD FOR TESTING.
SUBROUTINE IC_TEST(Y)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: Y
	DOUBLE PRECISION :: rand_1
	INTEGER :: i

	DO i=1,SIZE(Y)
		CALL random_number(rand_1)
		Y(i)=rand_1
	END DO

END SUBROUTINE IC_TEST




!************************************************************************
!This returns a normal distribution.  Taken from http://www.sdsc.edu/~tkaiser/f90.html .

function normal(mean,sigma) 
        implicit none
	DOUBLE PRECISION, parameter :: pi = 3.141592653589793239D0
        DOUBLE PRECISION :: normal,tmp, ran1, ran2
        DOUBLE PRECISION ::  mean,sigma
        INTEGER :: flag
        DOUBLE PRECISION :: fac,gsave,rsq,r1,r2
        save flag,gsave
        data flag /0/
        if (flag.eq.0) then
        rsq=2.0D0
        do while(rsq.ge. 1.0D0 .or. rsq.eq. 0.0D0)
		CALL random_number(ran1)
		CALL random_number(ran2)
                r1=2.0D0*ran1-1.0D0
                r2=2.0D0*ran2-1.0D0
                rsq=r1*r1+r2*r2
        enddo
            fac=sqrt(-2.0D0*log(rsq)/rsq)
            gsave=r1*fac
            tmp=r2*fac
            flag=1
        else
            tmp=gsave
            flag=0
        endif
        normal=tmp*sigma+mean
      	RETURN
end function normal


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


SUBROUTINE mcmc_eqen(ICtable)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: ICtable
	DOUBLE PRECISION :: Y(5), Y2(5), a(9), norm(8), pert(8), params(4), params2(4)
	DOUBLE PRECISION :: lam_4, V, mag, KE, en_left, check
	DOUBLE PRECISION :: sigma(8)
	DOUBLE PRECISION :: tol
	INTEGER :: iccounter, i, j, k, k_1, k_2, l, i_1, counter
	LOGICAL :: onsurf, bounds

	!Load tolerance and normal's sigma.
	tol = 1D-6
	sigma(1) = 10D0			!Phi
	sigma(2) = 10D0			!Psi
	sigma(3) = 100D0		!Phi_dot
	sigma(4) = 100D0		!Psi_dot
	sigma(5) = 100D0		!Lambda
	sigma(6) = 10D0			!M
	sigma(7) = 10D0			!mu
	sigma(8) = 10D0			!nu

	!First pick a point on eq en surface.
	iccounter = 0
	counter = 0
	CALL parameters_hybrid()
	!Load params vect.
	params(1)=Lambda
	params(2)=M
	params(3)=mu
	params(4)=nu	

	CALL SEED_IC_EQEN(Y)
PRINT*,"Y0",Y
PRINT*,"PARAMS0",params

	DO WHILE (counter < SIZE(ICtable,1)+10000)

		!Define the normal vector at Y.
		!a_1*(f_dot)+a_2*(s_dot)+a_3*(f)+a_4*(s)+
		!a_5(l)+a_6*(m)+a_7*(z)+a_8*(n)+a_9=0
		lam_4 = params(1)**4
		V = lam_4*((1D0-((Y(3)*Y(3))/(params(2)*params(2))))**2D0+&
			&((Y(2)*Y(2))/(params(3)*params(3)))+&
			&((Y(2)*Y(2)*Y(3)*Y(3))/(params(4)**4D0)))
		a(1) = Y(4)
		a(2) = Y(5)
		a(3) = 2D0*lam_4*Y(2)*(1D0/(mu*mu)+(Y(3)*Y(3))/(nu*nu*nu*nu))
		a(4) = 2D0*lam_4*Y(3)*((Y(2)*Y(2))/(nu*nu*nu*nu)-&
			&(2D0*Y(3)*Y(3))/(M*M*M*M)-2D0/(M*M))
		a(5) = (4D0/lambda)*V
		a(6) = ((4D0*lam_4*Y(2)*Y(2))/(M*M*M))*(1D0 + (Y(2)*Y(2))/(M*M))
		a(7) = -2D0*lam_4*Y(2)*Y(2)/(mu*mu*mu)
		a(8) = -4D0*lam_4*Y(2)*Y(2)*Y(3)*Y(3)/(nu*nu*nu*nu*nu)
		a(9) = .5D0*(Y(4)*Y(4)+Y(5)*Y(5)) + V - energy_scale**4
		!Normalize.
		mag = SQRT(a(1)*a(1)+a(2)*a(2)+a(3)*a(3)+a(4)*a(4)+a(5)*a(5)+&
			&a(6)*a(6)+a(7)*a(7)+a(8)*a(8))
		DO i=1,8
			norm(i) = a(i)/mag
		END DO

		!Perturb each parameter by Gaussian with mean=0.
		DO j=1,8
			pert(j) = normal(0D0, sigma(j))
		END DO

		!Create new point.
		DO k=2,5	
			Y2(k) = Y(k) + pert(k-1)
		END DO
		DO l=1,4
			params2(l) = params(l)+pert(l+4)
		END DO

		!Project onto the eq en slice.
PRINT*,"Calling proj."
		CALL proj_eqen(Y2,Y,params2,params,norm,tol)

		!Double check if on eq en surface
		V = (params2(1)**4)*((1D0-((Y2(3)*Y2(3))/(params2(2)*params2(2))))**2D0 +&
			&((Y2(2)*Y2(2))/(params2(3)*params2(3))) + &
			&((Y2(2)*Y2(2)*Y2(3)*Y2(3))/(params2(4)**4)))
		KE = .5D0*(Y2(4)*Y2(4)+Y2(5)*Y2(5))
		en_left = energy_scale**4 - KE - V
		IF(ABS(en_left) .LE. tol) THEN
			onsurf = .TRUE.
		ELSE
			onsurf = .FALSE.
PRINT*,"FAIL not on surf"
			CYCLE
		END IF
		!Test if within bounds. 
		IF(Y2(2)<phi_min .OR. Y2(2)>phi_max .OR. Y2(3)<phi_min .OR. Y2(3)>phi_max &
			& .OR. params2(1)<lambda_min .OR. params2(1)>lambda_max &
			&.OR. params2(2)<m_min .OR. params2(2)>m_max .OR. &
			&params2(3)<mu_min .OR. params2(3)>mu_max .OR.&
			&params2(4)<nu_min .OR. params2(4)>nu_max) THEN
			bounds = .FALSE.
PRINT*,"FAIL not in bounds"
			CYCLE
		ELSE
			bounds = .TRUE.
		END IF
		!If new point on surf, upgrade to new point.
PRINT*,"SUCCESS"
		counter = counter + 1
		Y=Y2
		params = params2
		!Burn in period.
		IF(counter>10000) THEN
			DO k_1=1,5
				ICtable(counter-10000,k_1)=Y(k_1)
			END DO
			DO k_2=1,4
				ICtable(counter-10000,k_2+5)=params(k_2)
			END DO
		END IF
		
	END DO
	


END SUBROUTINE mcmc_eqen



!Necessary subroutine for FKINSOL.

SUBROUTINE FKFUN (U, FVAL, IER)
IMPLICIT NONE

	DOUBLE PRECISION :: U(*), FVAL(*)
	INTEGER :: IER
	
	FVAL(1) = U(1) - 1D0 +SIN(U(1))*COS(U(1)) - 2D0*COS(U(1))

END SUBROUTINE FKFUN


!**********************************************************************
!Function that takes an m-D vector and gives the projection of that vector onto the (m-1)-D tangent space with normal, norm.  If the normal is a unit coordinate, then the normal direction will have a value of zero.  Euclidean.

FUNCTION projection(vect, norm)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:) :: vect, norm
	DOUBLE PRECISION, DIMENSION(SIZE(vect)) :: projection
	
	projection = vect - norm*dot_product(vect,norm)

END FUNCTION projection



!*************************************************************************
!NOTE: Not working properly...

!Subroutine which takes the perturbed vector Y=Y0+pert and params=params0+pert and projects them down the normal direction until they lie within tol of the equal energy slice.

SUBROUTINE proj_eqen(Y,Y0,params,params0,normal,tol)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: Y, params
	DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: normal, Y0, params0
	DOUBLE PRECISION, DIMENSION(SIZE(normal)) :: norm
	DOUBLE PRECISION, INTENT(IN) :: tol
	DOUBLE PRECISION, DIMENSION(SIZE(Y)-1+SIZE(params)) :: vect, vect_norm
	DOUBLE PRECISION :: alpha, V0, KE, en_left, en_left2, check, checknorm
	INTEGER :: i_1, i_2, i_3, i_4, i_5, i_6, i_7, counter
	
	norm = normal
	counter = 0

	!Load perturbed vector.
	DO i_1=1,4
		vect(i_1)=Y(i_1+1)
	END DO
	DO i_2=1,4
		vect(i_2+4)=params(i_2)
	END DO
	!Load initial vector.
	DO i_6=1,4
		vect_norm(i_6) = Y0(i_6+1)
	END DO
	DO i_7=1,4
		vect_norm(i_7+4) = params0(i_7)
	END DO

	!Check if on eq en slice.
	en_left = en_diff(vect)
PRINT*,"Projection en left1 ", en_left
	IF(ABS(en_left)<=tol) RETURN

	!Determine which way normal points and adjust to move onto eq en slice.
	vect_norm = vect_norm + norm
	checknorm = en_diff(vect_norm)
PRINT*,"Checknorm ",checknorm
	IF(checknorm*en_left<0)THEN
		norm = (-1D0)*norm
	END IF

	!Step size.
	alpha = en_left*.001D0

	!Move to eq en slice.
	DO i_5=1,500000
		!Translate down norm by alpha.
		vect = vect - alpha*norm
		
		en_left2 = en_diff(vect)

!PRINT*,"Projection en left 2", en_left2
		
		!Check for overshoot.
		IF(en_left2*en_left<0)THEN
			en_left = en_left2
			norm = -1D0*norm
			alpha = alpha*.001D0
PRINT*,"Overshoot test"
		END IF

		IF(ABS(en_left2)<=tol) EXIT 
		
	END DO

	!Rebuild Y and params.
	DO i_3=1,4
		Y(i_3+1)=vect(i_3)
	END DO
	DO i_4=1,4
		params(i_4)=vect(i_4+4)
	END DO
PRINT*,"Y ",Y
PRINT*,"Params",params	

END SUBROUTINE proj_eqen





!******************************************************************************************
!Subroutine which seeds the MCMC routine somewhere on the equal energy slice.

SUBROUTINE SEED_IC_EQEN(Y)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: Y
	DOUBLE PRECISION ::  rand_1, rand_2, V, left

	DO
	Y(1)=0D0
	Y(4)=0D0
	CALL random_number(rand_1)
	CALL random_number(rand_2)
	Y(2) = rand_2 *(phi_max-phi_min)+phi_min
	Y(3) = rand_1 *(psi_max-psi_min)+psi_min
	V = (lambda**4D0)*((1D0-((Y(3)*Y(3))/(m*m)))**2D0 +((Y(2)*Y(2))/(mu*mu)) + &
		&((Y(2)*Y(2)*Y(3)*Y(3))/ (nu**4D0)))
	left = 2D0*(energy_scale**4 - V)
	IF(left<0) THEN
		CYCLE
	ELSE
		Y(5)=SQRT(left)
		EXIT
	END IF
	END DO

END SUBROUTINE SEED_IC_EQEN


!*****************************************************************************************
!Function which will return the amount of difference in energy between a point in R^n and the surface E^4.

DOUBLE PRECISION FUNCTION en_diff(vect)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(8) :: vect
	DOUBLE PRECISION :: V0, KE

	V0 = (vect(5)**4)*((1D0-((vect(2)*vect(2))/(vect(6)*vect(6))))**2D0 +&
		&((vect(1)*vect(1))/(vect(7)*vect(7))) + &
		&((vect(1)*vect(1)*vect(2)*vect(2))/(vect(8)**4)))
	KE = .5D0*(vect(3)*vect(3)+vect(4)*vect(4))
	en_diff = energy_scale**4 - KE -V0

END FUNCTION en_diff 






END MODULE d_hybrid_initialconditions
