find_package(UMFPACK QUIET)
if(UMFPACK_FOUND)
  set(HAVE_UMFPACK true)
else()
  #check manually
  find_file (UMFPACK_H umfpack.h HINTS /usr/include/suitesparse DOC "Path to the UMFPACK header file umfpack.h")
  find_library (UMFPACK_LIBRARY umfpack DOC "Path to the UMFPACK library.")
  set (HAVE_UMFPACK false)
  if (UMFPACK_H AND UMFPACK_LIBRARY)
    message (STATUS "Umfpack header and library found.")
    message (STATUS "  Search now for required headers and libraries.")
    get_filename_component (UMFPACK_INCLUDE_DIR_ONLY ${UMFPACK_H} PATH)
    set (UMFPACK_INCLUDE_DIRS ${UMFPACK_INCLUDE_DIR_ONLY})
    set (UMFPACK_LIBRARIES ${UMFPACK_LIBRARY})

    # check for UFConfig.h
    # find_file (UFCONFIG_H UFconfig.h 
    #            HINTS /usr/include/suitesparse ${UMFPACK_INCLUDE_DIR_ONLY}
    # 	       DOC "Path to the header file UFconfig.h")
    # if (UFCONFIG_H)
    #     get_filename_component (UFCONFIG_DIR ${UFCONFIG_H} PATH)
    # 	if (NOT (UFCONFIG_DIR STREQUAL UFCONFIG_DIR) )
    #         set (UMFPACK_INCLUDE_DIRS "${UMFPACK_INCLUDE_DIR};${UFCONFIG_DIR}")
    # 	    message (STATUS "UMFPACK_INCLUDE_DIRS are now ${UMFPACK_INCLUDE_DIRS}")
    # 	else()
    # 	    message (STATUS "UFconfig.h and umfpack.h in same directory")
    # 	endif()
    # endif()
    
    # set (UFCONFIG_HINT "" CACHED)
    # set (AMD_HINT "" CACHED)
    # Mandatory header files
    set (UMFPACK_NEEDS_HEADERS "UFconfig.h;amd.h")    
    foreach (FNAME ${UMFPACK_NEEDS_HEADERS})
        get_filename_component(FNAME_WE ${FNAME} NAME_WE)
	# message(STATUS "FNAME_WE = ${FNAME_WE}")
        string(TOUPPER "${FNAME_WE}" VARNAME)
	# set(${VARNAME}_HINT "" CACHE STRING "Path to the header file ${FNAME}")
	# message(STATUS "Search for ${FNAME}, hint [${VARNAME}_HINT] = ${${VARNAME}_HINT}.")
        find_file (${VARNAME}_H ${FNAME}
                   HINTS /usr/include/suitesparse ${UMFPACK_INCLUDE_DIR_ONLY} # ${VARNAME}_HINT
    		   DOC "Absolute file name for header file ${FNAME}.")
        if (${VARNAME}_H)
	    message (STATUS "${FNAME} found.")
            get_filename_component (DIR ${${VARNAME}_H} PATH)
	    # set(DIR "/dings/${VARNAME}/bums") # only for testing comparison
	    # message (STATUS "Found ${FNAME} in ${DIR}.")
    	    if (NOT (UMFPACK_INCLUDE_DIR_ONLY STREQUAL DIR) )
                set (UMFPACK_INCLUDE_DIRS "${UMFPACK_INCLUDE_DIRS};${DIR}")
    		# message (STATUS "UMFPACK_INCLUDE_DIRS are now ${UMFPACK_INCLUDE_DIRS}")
    	    else()
    	        # message (STATUS "${FNAME} and umfpack.h in same directory")
    	    endif()
	else()
	    message (STATUS "${FNAME} not found.")
	    message (STATUS "  You can set the absolute file name in ${VARNAME}_H.")
        endif()
    endforeach()

    # Libraries that are needed but may be already statically linked into libumfpack.
    set (UMFPACK_WANTS_LIBS "amd;blas")
    foreach (LNAME ${UMFPACK_WANTS_LIBS})
        string (TOUPPER "${LNAME}" VARNAME)
	find_library (${VARNAME}_LIBRARY ${LNAME} DOC "Absolute file name for library ${LNAME}.")
	if (${VARNAME}_LIBRARY)
	    message (STATUS "Library ${LNAME} found.")
	    set (UMFPACK_LIBRARIES "${UMFPACK_LIBRARIES};${${VARNAME}_LIBRARY}")
	else()
	    message (STATUS "Library ${LNAME} not found.")
	    message (STATUS "  You can set the absolute file name in ${VARNAME}_LIBRARY")
	endif()
    endforeach()
    # message (STATUS "UMFPACK_LIBRARIES are ${UMFPACK_LIBRARIES}.")

    try_compile(UMFPACK_RESULT ${CMAKE_BINARY_DIR}/CMakeTmp 
                ${MTL_DIR}/tools/cmake/umfpack_test.cpp
		CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${UMFPACK_INCLUDE_DIRS} " 
		CMAKE_FLAGS "-DLINK_LIBRARIES=${UMFPACK_LIBRARIES}"
		COMPILE_DEFINITIONS " " 
		OUTPUT_VARIABLE errors)
    if (UMFPACK_RESULT)
        set(HAVE_UMFPACK true)
    else()
        message(WARNING "Could not compile umfpack test because:\n${errors}")
    endif()



# 	get_filename_component(UMFPACK_INCLUDE_DIR ${UMFPACK_H} PATH)
# 	# message(STATUS "UMFPACK_H = ${UMFPACK_H}")
# 	# message(STATUS "UMFPACK_INCLUDE_DIRS = ${UMFPACK_INCLUDE_DIRS}")
# 	find_file(UFCONFIG_H UFconfig.h HINTS /usr/include/suitesparse DOC "Path to the UMFPACK header file umfpack.h")
	



# 	find_library(AMD_LIBRARY amd)
# 	find_library(BLAS_LIBRARY blas)
# 	if(AMD_LIBRARY AND BLAS_LIBRARY)
# 		set(UMFPACK_LIBRARIES "${UMFPACK_LIBRARY};${AMD_LIBRARY};${BLAS_LIBRARY}")
# 		set(HAVE_UMFPACK true)
# 	endif()
# 	message(STATUS "UMFPACK_LIBRARIES = ${UMFPACK_LIBRARIES}")
# 	try_compile(UMFPACK_RESULT ${CMAKE_BINARY_DIR}/CMakeTmp 
# 	            ${MTL_DIR}/tools/cmake/umfpack_test.cpp
# 		    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${UMFPACK_INCLUDE_DIRS} " 
# 		    CMAKE_FLAGS "-DLINK_LIBRARIES=${UMFPACK_LIBRARIES}"
# 		    COMPILE_DEFINITIONS " " 
# 		    OUTPUT_VARIABLE errors)
# 	if (NOT UMFPACK_RESULT)
# 	  message(WARNING "Could not compile umfpack test because: ${errors}")
# 	endif()
# 
  else()
     message(STATUS "Umfpack not found. You can set UMFPACK_H and UMFPACK_LIBRARY by hand.")
  endif()
endif()
