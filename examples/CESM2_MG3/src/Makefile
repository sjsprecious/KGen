include ./Macros

PCOLS := 384
PVER  := 32

# uncommand followings to build with PAPI library
#PAPI_DIR := /path/to/PAPI
#PAPI_LIB := ${PAPI_DIR}/lib/libpapi.a

FC := pgf90 
FC_FLAGS += -DPCOLS=${PCOLS} -DPVER=${PVER}

run: build
	./mg3_kernel.exe

build: mg3_kernel.exe

mg3_kernel.exe: program.o micro_mg_cam.o micro_mg3_0.o \
		micro_mg_utils.o shr_spfn_mod.o \
		shr_const_mod.o shr_kind_mod.o wv_sat_methods.o
	${FC} ${LD_FLAGS} -o $@ $^ ${PAPI_LIB}

program.o: program.F90 micro_mg_cam.o
	${FC} ${FC_FLAGS} -c -o $@ $<

micro_mg_cam.o: micro_mg_cam.F90 micro_mg3_0.o shr_kind_mod.o
	${FC} ${FC_FLAGS} -c -o $@ $<

micro_mg3_0.o: micro_mg3_0.F90 micro_mg_utils.o \
		shr_spfn_mod.o wv_sat_methods.o 
	${FC} ${FC_FLAGS} -c -o $@ $<

micro_mg_utils.o: micro_mg_utils.F90 shr_spfn_mod.o
	${FC} ${FC_FLAGS} -c -o $@ $<

wv_sat_methods.o: wv_sat_methods.F90 
	${FC} ${FC_FLAGS} -c -o $@ $<

shr_spfn_mod.o: shr_spfn_mod.F90 shr_const_mod.o shr_kind_mod.o
	${FC} ${FC_FLAGS} -c -o $@ $<

shr_const_mod.o: shr_const_mod.F90 shr_kind_mod.o
	${FC} ${FC_FLAGS} -c -o $@ $<

shr_kind_mod.o: shr_kind_mod.F90
	${FC} ${FC_FLAGS} -c -o $@ $<

clean:
	rm -f *.exe *.o *.mod
