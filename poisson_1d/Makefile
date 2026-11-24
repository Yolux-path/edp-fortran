FC = gfortran

FFLAGS = -Wall -O2 -fcheck=all

# Source files
SRCS = donnees.f90 fonctions.f90 algebre.f90 schemas.f90 post_traitement.f90 DF.f90 

# Object files
OBJS = $(SRCS:.f90=.o)

# Executable
EXEC = exe

all: $(EXEC)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< 

$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) $^ -o $@

clean:
	rm -f $(OBJS) $(EXEC)