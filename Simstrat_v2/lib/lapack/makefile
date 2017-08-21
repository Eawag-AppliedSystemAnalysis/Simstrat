CC = x86_64-w64-mingw32-gfortran

LAPACK_MODULES := xerbla.f \
									dgtsv.f

%.o: %.f
	$(CC) -Ofast -c $^ -o $@ 

all: $(LAPACK_MODULES:.f=.o)
	@ar -qc liblapack.a $^
