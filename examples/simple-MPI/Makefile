KGEN_HOME := ../..
KGEN := ${KGEN_HOME}/bin/kgen

#PAPI_DIR := /path/to/PAPI

SRC_DIR := ${PWD}/src
SRC := ${SRC_DIR}/update_mod.F90
CALLSITE := update_mod:update:calc

test:
	${KGEN} \
		--mpi enable \
		--cmd-clean "cd ${SRC_DIR}; make -f Makefile.mpirun clean" \
		--cmd-build "cd ${SRC_DIR}; make -f Makefile.mpirun build" \
        --cmd-run "cd ${SRC_DIR}; make -f Makefile.mpirun run" \
        --data condition=and:"i == 2" \
		${SRC}:${CALLSITE}

# additiional flags
#        --repr-code enable \
#        --repr-etime nbins=5,ndata=20 \
#        --repr-papi header=${PAPI_DIR}/include/f90papi.h,static=${PAPI_DIR}/lib/libpapi.a,event=PAPI_TOT_INS \
#        --source alias=/path/a:path/b # add this option if your system requires alias paths
#        --rebuild all \

clean:
	${MAKE} -f Makefile.mpirun clean -C src
	rm -rf kernel state kgen.log _kgen_compflag_cmdwrapper.sh strace.log include.ini model model.ini elapsedtime coverage papi
