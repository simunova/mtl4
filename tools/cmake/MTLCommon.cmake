############################
#This configuration file defines some cmake variables:
#MTL_INCLUDE_DIRS: list of include directories for the mtl
#MTL_LIBRARIES: libraries needed for interfaces like umfpack and arprec, see below
#MTL_CXX_DEFINITIONS: definitions to enable the requested interfaces
#MTL_VERSION: version (current: 4)
#MTL_MINOR_VERSION: minor version 
#
#supported components:
#Umfpack, Arprec

option(ENABLE_OPENMP "switch on to enable OpenMP flags for mtl" OFF)
option(ENABLE_SHORT_ELE_PROD "enable short notation for element-wise product" OFF)
option(ENABLE_CXX_ELEVEN "enable C++11 features as far as compiler permits" ON)
option(USE_ASSERTS "Use assert instead of throwing exceptions" ON)


unset(MTL_LIBRARIES )
unset(MTL_CXX_DEFINITIONS )
unset(MTL_INCLUDE_DIRS )

find_package(Boost 1.40 REQUIRED)
if(Boost_FOUND)
	LIST(APPEND MTL_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
endif(Boost_FOUND)

include(${MTL_DIR}/tools/cmake/Vampir.cmake)
include(${MTL_DIR}/tools/cmake/UMFPACK.cmake)
include(${MTL_DIR}/tools/cmake/ARPREC.cmake)

find_package(Subversion)
set(MTL_VERSION "4")
if(Subversion_FOUND)
  if(EXISTS "${CMAKE_SOURCE_DIR}/.svn")
  	Subversion_WC_INFO(${MTL_DIR} mtlSubinfo)
  	set(MTL_MINOR_VERSION ${mtlSubinfo_WC_REVISION})
    #	message("current revision: ${mtlSubinfo_WC_REVISION}")
  else ()
    set(MTL_MINOR_VERSION "0")
  endif()
else(Subversion_FOUND)
	set(MTL_MINOR_VERSION "0")
endif(Subversion_FOUND)

if(EXISTS ${MTL_DIR}/tools/cmake/C++11Features.cmake)
  include(${MTL_DIR}/tools/cmake/C++11Features.cmake)
endif()

if (MSVC)
    add_definitions(/wd4522) # multiple assignment ops for single type, to be investigated further if avoidable
endif()

if (USE_ASSERTS)
  list(APPEND MTL_CXX_DEFINITIONS "-DMTL_ASSERT_FOR_THROW")
endif()

if(HAVE_UMFPACK)
	list(APPEND MTL_CXX_DEFINITIONS "-DMTL_HAS_UMFPACK")
	include_directories(${UMFPACK_INCLUDE_DIRS})
	list(APPEND MTL_LIBRARIES ${UMFPACK_LIBRARIES})
endif()
if(HAVE_ARPREC)
	list(APPEND MTL_CXX_DEFINITIONS "-DMTL_HAS_ARPREC")
	include_directories(${ARPREC_INCLUDE_DIRS})
	list(APPEND MTL_LIBRARIES ${ARPREC_LIBRARIES})
endif()
if(EXISTS ${MTL_DIR}/boost/numeric/mtl/mtl.hpp)
	list(APPEND MTL_INCLUDE_DIRS "${MTL_DIR}")
else()
	list(APPEND MTL_INCLUDE_DIRS "${MTL_DIR}/../../include")
endif(EXISTS ${MTL_DIR}/boost/numeric/mtl/mtl.hpp)

if(ENABLE_VAMPIR AND VAMPIR_FOUND)
#add_definitions("-DMTL_HAS_VPT -DVTRACE -vt:inst manual")
	list(APPEND MTL_CXX_DEFINITIONS "-DMTL_HAS_VPT" "-DVTRACE")
	list(APPEND MTL_CXX_DEFINITIONS ${VT_COMPILE_FLAGS})
	set(MTL_LINK_FLAGS "${MTL_LINK_FLAGS} ${VT_LINK_FLAGS}")
	if(EXISTS ${MTL_DIR}/boost/numeric/mtl/interface/vpt.cpp)
		set(VAMPIR_SRC ${MTL_DIR}/boost/numeric/mtl/interface/vpt.cpp)
	else()
		set(VAMPIR_SRC ${MTL_DIR}/vpt.cpp)
	endif()
	get_target_property(MTL_VAMPIR_ADDED mtl_vampir TYPE)
	if(NOT MTL_VAMPIR_ADDED)
		add_library(mtl_vampir ${VAMPIR_SRC})
	endif()
	list(APPEND MTL_LIBRARIES ${VT_LIBRARIES} mtl_vampir)
	set(HAVE_VAMPIR TRUE)
else()
	set(HAVE_VAMPIR FALSE)
endif(ENABLE_VAMPIR AND VAMPIR_FOUND)

if(ENABLE_OPENMP )
	find_package(OpenMP REQUIRED)
	if(OPENMP_FOUND)
		list(APPEND MTL_CXX_DEFINITIONS "-DMTL_WITH_OPENMP")
		list(APPEND MTL_CXX_DEFINITIONS ${OpenMP_CXX_FLAGS})
		list(APPEND MTL_LIBRARIES ${OpenMP_CXX_FLAGS})
	else()
		message(FATAL_ERROR "OpenMP not found")
	endif()
endif()
message(STATUS "MTL Find components: ${MTL_FIND_COMPONENTS}")
#we found nothing..
set(MTL_NOT_FOUND )
#remove?
foreach(CURCOMP ${MTL_FIND_COMPONENTS})
#look for a file called COMPONENT.cmake in the mtl-directory (/usr/share/mtl/)
	string(TOUPPER ${CURCOMP} CURCOMP_UPPER)
	set(curfile "${MTL_DIR}/tools/cmake/${CURCOMP_UPPER}.cmake")
	if(EXISTS ${curfile})
		include(${curfile})
		#look for component 
		#check if the component was correctly found
		if(HAVE_${CURCOMP_UPPER})
			list(APPEND MTL_INCLUDE_DIRS ${${CURCOMP_UPPER}_INCLUDE_DIRS})
			list(APPEND MTL_LIBRARIES ${${CURCOMP_UPPER}_LIBRARIES})
			list(APPEND MTL_CXX_DEFINITIONS "-DMTL_HAS_${CURCOMP_UPPER}")
		else()
			list(APPEND MTL_NOT_FOUND ${CURCOMP})
		endif()
	else()
		list(APPEND MTL_NOT_FOUND ${CURCOMP})
	endif()
endforeach()

if(MTL_FIND_REQUIRED AND MTL_NOT_FOUND)
	message(SEND_ERROR "could not find all components: ${MTL_NOT_FOUND}")
endif()
#include_directories(${MTL_INCLUDE_DIRS})

macro(mtl_check_cxx_compiler_flag FLAG RESULT)
  # counts entirely on compiler's return code, maybe better to combine it with check_cxx_compiler_flag
  file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx" "int main() { return 0;}\n")
  try_compile(${RESULT}
    ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx
    COMPILE_DEFINITIONS ${FLAG})  
endmacro()

