FC = gfortran
FFLAGS = -ffixed-form -ffixed-line-length-none -fcheck=bounds -g

# OBJS = global_data.o main.o calcu.o calcv.o calcp.o grid.o init.o
OBJS = global_data.o main.o fs.o io.o boundary.o grid_init.o	fs_sor.o fs_calc.o simple.o simple_calc.o
cavity: $(OBJS)
	$(FC) $(OBJS) -o cavity

%.o: %.f
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod cavity
	