ROW := 16
COL := 16

FC := gfortran
FC_FLAGS := -DROW=${ROW} -DCOL=${COL}

run: build
	export OMP_NUM_THREADS=8; ./demo.exe

build: demo.exe

demo.exe: program.o update_mod.o calc_mod.o
	${FC} -fopenmp -o $@ $^

program.o: program.F90 update_mod.o
	${FC} -fopenmp -c -o $@ $<

update_mod.o: update_mod.F90 calc_mod.o
	${FC} -fopenmp ${FC_FLAGS} -c -o $@ $<

calc_mod.o: calc_mod.F90
	${FC} -fopenmp -c -o $@ $<

clean:
	rm -f *.exe *.o *.mod
