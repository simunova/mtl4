CXXFLAGS += -I${MTL_BOOST_ROOT} -I${BOOST_ROOT} 

# Other possible CXXFLAGS

# To see which power function was chosen
CXXFLAGS += -DMTL_TRACE_POWER_DISPATCHING

# For modular arithmetic to see which invertibility test was chosen
CXXFLAGS += -DMTL_TRACE_MOD_N_INVERTIBILITY_DISPATCHING 

# CXXFLAGS += -g
# CXXFLAGS += -O2

# CXXFLAGS += -Wno-long-double

# To use pre-compiled headers (can be put in conceptg++ script)
# CXXFLAGS += --include bits/stdc++.h                               

# For some strange problems on lunkwill
# CXXFLAGS += -L/u/pgottsch/software/conceptgcc/darwin/lib/ppc64/ -Winvalid-pch -m64


XTRA_SOURCES = algebraic_functions.hpp ${MTL_BOOST_ROOT}/boost/numeric/linear_algebra/concepts.hpp \
                 Makefile
 
CXX =  conceptg++
#CXX =  g++

TARGETS =	simple_test2 simple_test positive_real	positive_real_explicit positive_real_power \
		positive_real_power_explicit string_example \
		float_example int_example mod_n_example mod_n_power accumulation_example \
		concept_map_test      int_example  accumulation_measuring

it:	simple_test


all:	${TARGETS}


.cpp:	$(@:=.cpp) $(XTRA_SOURCES)
	${CXX} ${CXXFLAGS} -o $@ $?


clean:
	rm -rf ${TARGETS} *~