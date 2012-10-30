SHELL = /bin/bash

prefix       = /usr/local
exec_prefix  = ${prefix}
includedir   = ${prefix}/include
libdir       = ${exec_prefix}/lib
genmoddir = /home/lpri691/Cosmology_Research/F_projects/libraries/genmodules

sourcedir= /home/lpri691/Cosmology_Research/Fortran/clesse_copy/main_double/cvodeintegrator

F77         = mpif90 #gfortran
FFLAGS      = -O2 -pg
F77_LNKR    = mpif90 #gfortran
F77_LDFLAGS = -O2 -pg
F77_LIBS    = #  -L/usr/lib/gcc/i686-linux-gnu/4.6 -L/usr/lib/gcc/i686-linux-gnu/4.6/../../../i386-linux-gnu -L/usr/lib/gcc/i686-linux-gnu/4.6/../../../../lib -L/lib/i386-linux-gnu -L/lib/../lib -L/usr/lib/i386-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/i686-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath -lm 

LIBRARIES =  -lsundials_fcvode -lsundials_cvode -lsundials_fnvecserial -lsundials_nvecserial # -lsundials_fnvecparallel -lsundials_nvecparallel ${LIBS}
LIBRARIES_BL = -L/usr/lib64  -llapack -lblas -lm #  -L/usr/lib/gcc/i686-linux-gnu/4.6 -L/usr/lib/gcc/i686-linux-gnu/4.6/../../../i386-linux-gnu -L/usr/lib/gcc/i686-linux-gnu/4.6/../../../../lib -L/lib/i386-linux-gnu -L/lib/../lib -L/usr/lib/i386-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/i686-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath

LD_FLAGS = -llapack -lblas -lm

EXAMPLES = hybrid_integrator_d

OBJECTS = ${EXAMPLES:=.o} d_hybrid_initialconditions.o

# -----------------------------------------------------------------------------------------

.SUFFIXES : .o .f

.f.o :
	${F77} ${FFLAGS} -c $^

%.o:%.mod


hybrid_integrator_d.exe: hybrid_integrator_d.o d_hybrid_initialconditions.o hybrid_metropolis.o ${genmoddir}/*.o hybrid_subroutines.o 
	${F77_LNKR} -L${sourcedir} -o hybrid_integrator_d hybrid_integrator_d.o hybrid_subroutines.o d_hybrid_initialconditions.o  hybrid_metropolis.o ${genmoddir}/*.o ${F77_LDFLAGS} ${F77_LIBS} -L${libdir} ${LIBRARIES} ${LIBRARIES_BL} 


hybrid_integrator_d.o: hybrid_integrator_d.f90 d_hybrid_initialconditions.o hybrid_metropolis.o hybrid_subroutines.o ${genmoddir}/*.o 
	${F77_LNKR} -I$(genmoddir) -c hybrid_integrator_d.f90 ${F77_LDFLAGS} ${F77_LIBS} -L${libdir} ${LIBRARIES} ${LIBRARIES_BL} -L${sourcedir}

d_hybrid_initialconditions.o: d_hybrid_initialconditions.f90 ${genmoddir}/*.o
	${F77_LNKR} -I$(genmoddir) -c d_hybrid_initialconditions.f90 ${F77_LDFLAGS} ${F77_LIBS} -L${libdir} ${LIBRARIES} ${LIBRARIES_BL} -L${sourcedir}

hybrid_subroutines.o: hybrid_subroutines.f90 d_hybrid_initialconditions.o hybrid_metropolis.o
	${F77_LNKR} -I$(genmoddir) -c hybrid_subroutines.f90 ${F77_LDFLAGS} ${F77_LIBS} -L${libdir} ${LIBRARIES} ${LIBRARIES_BL} -L${sourcedir}

hybrid_metropolis.o: hybrid_metropolis.f90 d_hybrid_initialconditions.o ${genmoddir}/*.o
	${F77_LNKR} -I$(genmoddir) -c hybrid_metropolis.f90 ${F77_LDFLAGS} ${F77_LIBS} -L${libdir} ${LIBRARIES} ${LIBRARIES_BL} -L${sourcedir}


clean:
	rm -f ${OBJECTS}
	rm -f ${EXAMPLES}
	rm -f *~
	rm -f fort.*
	rm -f fail* succ*
	rm -f info.*
	rm -f *.o *.mod
	rm -f traj*.bin
	rm -f lyapunov_integrator

clean2:
	rm -f fort.*
	rm -f fail* succ*
	rm -f info.*
	rm -f traj*.bin

# -----------------------------------------------------------------------------------------

