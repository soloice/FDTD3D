cc=mpicxx
cflags=
FDTD_MPI : main.cpp updateEH.cpp input.cpp output.cpp 
	$(cc) $(cflags) -o FDTD_MPI main.cpp updateEH.cpp input.cpp output.cpp

