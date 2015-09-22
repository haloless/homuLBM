

# COORD := CART
# COORD := RZ


# BDIR := ./bin

# NDIM := 2
NDIM := 3

#
COPTS     = -O2
CINCLUDES = 
CDEFS     = 
#
CFLAGS    = $(COPTS) $(CINCLUDES) $(CDEFS)

#
F90OPTS   = -O2
#
F90FLAGS  = $(F90OPTS)

#
LIBS      = 
LFLAGS    = $(LIBS) -lstdc++ -lgfortran


EXESRC := ShanChen$(NDIM)D.F90
EXEOBJ := $(EXESRC:%.F90=%.o)

# F90
F90SRC := commons.F90 
#F90SRC += ft_color.F90 ft_tension.F90 ft_advance.F90
#F90SRC += ns_conv.F90 ns_proj.F90 ns_bndry.F90
#
F90OBJ := $(F90SRC:%.F90=%.o)

# CPP
CPPSRC := output.cpp 
CPPOBJ := $(CPPSRC:%.cpp=%.o)

all: fortpart cpart $(EXEOBJ)
	g++ $(F90OBJ) $(CPPOBJ) $(EXEOBJ) $(LFLAGS)



fortpart: fortecho $(F90OBJ)

fortecho:
	@echo "compiling <"$(F90SRC)"> to <"$(F90OBJ)">"




cpart: cecho $(CPPOBJ)

cecho:
	@echo "compiling <"$(CPPSRC)"> to <"$(CPPOBJ)">"


# rules
%.o: %.F90
	gfortran $(F90FLAGS) -c $< -o $@
%.o: %.cpp
	g++ $(CFLAGS) -c $< -o $@



.PHONY: clean
clean:
	rm -f a.exe *.o *.mod





