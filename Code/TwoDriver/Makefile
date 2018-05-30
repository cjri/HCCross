CC	      = g++
CC_FLAGS	= -g3 -O3 -Wall -I/usr/local/include/
LD_FLAGS	= -L/usr/local/lib -lgsl -lgslcblas -lm
SE_OBJECTS      = simulation_run_exp.o utilities_comp.o io.o
TD_OBJECTS	= simulation_run_exp_twodriver.o utilities_comp.o io.o
WRS_OBJECTS     = simulation_run_exp_wrs.o utilities_comp.o io_wrs.o
WT_OBJECTS	= simulation_run_exp_wrs_twodriver.o utilities_comp.o io_wrs.o

sdriver:  $(SE_OBJECTS)
	$(CC) $(CC_FLAGS) $(SE_OBJECTS) -o comp_sim_exp $(LD_FLAGS)
simulation_run_exp.o: simulation_run_exp.cpp
	$(CC) $(CC_FLAGS) -c simulation_run_exp.cpp
utilities_comp.o: utilities_comp.cpp
	$(CC) $(CC_FLAGS) -c utilities_comp.cpp
io.o: io.cpp
	$(CC) $(CC_FLAGS) -c io.cpp

tdriver: $(TD_OBJECTS)
	$(CC) $(CC_FLAGS) $(TD_OBJECTS) -o comp_sim_exp_td $(LD_FLAGS)
simulation_run_exp_twodriver.o: simulation_run_exp_twodriver.cpp
	$(CC) $(CC_FLAGS) -c simulation_run_exp_twodriver.cpp

sdriver_wrs:    $(WRS_OBJECTS)
	$(CC) $(CC_FLAGS) $(WRS_OBJECTS) -o comp_sim_exp_wrs $(LD_FLAGS)
simulation_run_exp_wrs.o: simulation_run_exp_wrs.cpp
	$(CC) $(CC_FLAGS) -c simulation_run_exp_wrs.cpp
io_wrs.o: io_wrs.cpp
	$(CC) $(CC_FLAGS) -c io_wrs.cpp

tdriver_wrs:	$(WT_OBJECTS)
	$(CC) $(CC_FLAGS) $(WT_OBJECTS) -o comp_sim_exp_wrs_td $(LD_FLAGS)
simulation_run_exp_wrs_twodriver.o: simulation_run_exp_wrs_twodriver.cpp
	$(CC) $(CC_FLAGS) -c simulation_run_exp_wrs_twodriver.cpp

