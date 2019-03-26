project(mtl4_binaries)
if (NOT MTL_IS_CONFIGURED)
  find_package(MTL REQUIRED)
endif()

macro (compile_all testing)
    # cycle through the sources
    # The directory assumes that each cpp file is a separate test
    # so, create a executable target and an associated test target
    foreach (source ${ARGN})
        get_filename_component (test ${source} NAME_WE)
		string(REPLACE " " ";" newSource ${source})
		# message(STATUS "Add test ${test} from source ${newSource}.")
        add_executable (${test} ${newSource})
		if(MTL_LIBRARIES)
			target_link_libraries(${test} ${MTL_LIBRARIES})
		endif(MTL_LIBRARIES)
		#if(MTL_LINK_FLAGS)
		#  get_target_property(CURLINK ${test} LINK_FLAGS)
		# message (status "Link flags are ${CURLINK}")
		#  if(CURLINK)
		#    set_target_properties(${test} PROPERTIES LINK_FLAGS "${MTL_LINK_FLAGS} ${CURLINK}")
		#  else()
		#    set_target_properties(${test} PROPERTIES LINK_FLAGS "${MTL_LINK_FLAGS}")
		#  endif()
		#endif(MTL_LINK_FLAGS)
		# add_test (${test} ${EXECUTABLE_OUTPUT_PATH}/${test})
		if (${testing} STREQUAL "true")
		  # message("testing: ${testing} ${source}")
		  #add_test (${test} ${EXECUTABLE_OUTPUT_PATH}/${test})
		  add_test(${test} ${test})
		endif()
    endforeach (source)
endmacro (compile_all)

# enable test generation
find_program(VALGRIND_COMMAND valgrind)
if(VALGRIND_COMMAND)
  set(CTEST_MEMORYCHECK_COMMAND ${VALGRIND_COMMAND})
endif(VALGRIND_COMMAND)
set(BUILDNAME "${TESTBUILDNAME}")
enable_testing()
include(CTest)

if(ENABLE_GCOV)
  add_definitions("-fprofile-arcs -ftest-coverage")
  list(APPEND MTL_LIBRARIES "gcov")
  set(COVERAGE_EXTRA_FLAGS "-p -l") # has no effect, set "-p" by hand
  # if (${COVERAGE_EXTRA_FLAGS})
  #   set(COVERAGE_EXTRA_FLAGS "${COVERAGE_EXTRA_FLAGS}" "-p")
  # else()
  #   set(COVERAGE_EXTRA_FLAGS "-p")
  # endif()
	set(COVERAGE_EXTRA_FLAGS "${COVERAGE_EXTRA_FLAGS} -p")
endif(ENABLE_GCOV)

if (MSVC)
    add_definitions(/wd4244) # narrowing, only suppress warning in tests instead of hiding it by conversion (so that users never see it)
    if (HAVE_UMFACK)
	# add_definitions(/NODEFAULTLIB:libcmt.lib /NODEFAULTLIB:msvcrt.lib /NODEFAULTLIB:libcd.lib /NODEFAULTLIB:libcmtd.lib /NODEFAULTLIB:msvcrtd.lib) # no impact
	add_definitions(/NODEFAULTLIB:msvcrt.lib) # no impact
    endif()
endif()

# Descent into the test program directory to build the test suite
# MTL4 test programs
add_subdirectory("numeric/mtl/test")
# as well as the example programs
add_subdirectory("numeric/mtl/examples")
# ITL test programs
add_subdirectory("numeric/itl/test")
# Experimental and timing are not compiled by default
add_subdirectory("numeric/mtl/experimental" EXCLUDE_FROM_ALL)
add_subdirectory("numeric/mtl/timing" EXCLUDE_FROM_ALL)

