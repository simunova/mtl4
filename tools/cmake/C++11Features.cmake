include(CheckCXXCompilerFlag)

option(MTL_VERBOSE_TRYCOMPILE "Print error message when C++11 feature is not supported" OFF)

# Compiler flag for C++11: special case for VC only, everything else should be equal
# Might need later adaption, once compilers don't support this flag any longer
if(MSVC)
  set(CXX_ELEVEN_FLAG "")
else()
  set(CXX_ELEVEN_FLAG "-std=c++0x")
endif()

# only performed once (to reevaluated delete CMakeCache.txt)
check_cxx_compiler_flag("${CXX_ELEVEN_FLAG}" CXX_ELEVEN_FLAG_SUPPORTED) 

if (NOT CXX_ELEVEN_FLAG_SUPPORTED)
  if (ENABLE_CXX_ELEVEN)
    message("C++11 flag not supported by your compiler (probably too old).")
  endif()
  return()
endif()

if (NOT ENABLE_CXX_ELEVEN)
  return()
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" AND ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CXX_ELEVEN_FLAG "${CXX_ELEVEN_FLAG} -stdlib=libc++")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lc++")
endif()

message(STATUS "Add ${CXX_ELEVEN_FLAG}")
list(APPEND MTL_CXX_DEFINITIONS ${CXX_ELEVEN_FLAG})
#add_definitions("${CXX_ELEVEN_FLAG}")

# set (CXX_ELEVEN_FEATURE_LIST "MOVE" "AUTO" "RANGEDFOR" "INITLIST" "STATICASSERT" "DEFAULTIMPL")
file(GLOB CHECKS "${MTL_DIR}/tools/cmake/*_CHECK.cpp")

foreach (CURFILE ${CHECKS})
  get_filename_component(FNAME ${CURFILE} NAME_WE)
  string(REPLACE "_CHECK" "" feature ${FNAME})
  # message(STATUS "${feature}")
  try_compile(${feature}_RESULT ${CMAKE_BINARY_DIR} "${CURFILE}" COMPILE_DEFINITIONS "${CXX_ELEVEN_FLAG}" OUTPUT_VARIABLE errors)
  message(STATUS "Support C++11's ${feature} - ${${feature}_RESULT}")
  if (${feature}_RESULT)
    list(APPEND MTL_CXX_DEFINITIONS "-DMTL_WITH_${feature}")
  elseif (MTL_VERBOSE_TRYCOMPILE)
    message(STATUS "Failed because: ${errors}")
  endif()
endforeach()

# message(STATUS "C++11 flags: ${MTL_CXX_DEFINITIONS}")
