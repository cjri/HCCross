CC              = g++
CC_FLAGS        = -g3 -O3 -Wall -I  /usr/local/include/gsl/ -mmacosx-version-min=10.11
LD_FLAGS	= -L/usr/local/lib -lgsl -lgslcblas -lcblas -lm
SIM_OBJECTS	= simulation.o utilities.o

slocus: $(SIM_OBJECTS)
	$(CC) $(CC_FLAGS) $(SIM_OBJECTS) -o comp_sim $(LD_FLAGS)
simulation.o: simulation.cpp
	$(CC) $(CC_FLAGS) -c simulation.cpp
utilities.o: utilities.cpp
	$(CC) $(CC_FLAGS) -c utilities.cpp

