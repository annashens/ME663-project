FC = gfortran
FFLAGS = -ffixed-form -ffixed-line-length-none -fcheck=bounds -g

# OBJS = global_data.o main.o calcu.o calcv.o calcp.o grid.o init.o
OBJS = global_data.o main.o nas2d.o io.o boundary.o grid_init.o	pressure_solver.o
cavity: $(OBJS)
	$(FC) $(OBJS) -o cavity

%.o: %.f
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod cavity