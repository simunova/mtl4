project (mtl4)

cmake_minimum_required(VERSION 2.6)

#let cmake look for the svn-version
find_package(Subversion)
if(Subversion_FOUND)
	Subversion_WC_INFO(${CMAKE_CURRENT_BINARY_DIR} mtlSubinfo)
	set(CurrentRevision ${mtlSubinfo_WC_REVISION})
#	message("current revision: ${mtlSubinfo_WC_REVISION}")
else(Subversion_FOUND)
	set(CurrentRevision "0")
endif(Subversion_FOUND)
#################################################################################
######################## THE INSTALL RULES   ####################################
#################################################################################
#set the cmake-policy to suppress the warning. 
#don't follow symbolic links when globbing files recursive

cmake_policy(SET CMP0009 NEW)
INSTALL(DIRECTORY ../../boost
	DESTINATION include/
	PATTERN ".svn" EXCLUDE
       )

INSTALL(FILES ../../MTLConfig.cmake ../../README ../../tools/license/license.mtl.txt
	../../boost/numeric/mtl/interface/vpt.cpp
	DESTINATION share/mtl/
       )

#install the helper scripts in tools/cmake
FILE(GLOB CMAKE_HELPERS ../../tools/cmake/*.cmake)
INSTALL(FILES ${CMAKE_HELPERS} 
	DESTINATION share/mtl/tools/cmake
	)

#install helper sources in tools/cmake
FILE(GLOB CMAKE_HELPERS ../../tools/cmake/*.cpp)
INSTALL(FILES ${CMAKE_HELPERS} 
	DESTINATION share/mtl/tools/cmake
	)

#configure the uninstall script (will remove each hpp-file in boost/*, but leaves the empty directories)
FILE(GLOB_RECURSE ALL_HPPS_LIST RELATIVE ${CMAKE_CURRENT_BINARY_DIR}/../../ ../../boost/*.hpp)

foreach(curFile IN LISTS ALL_HPPS_LIST)
	#get all directories for the debian packaging
	get_filename_component(curDir ${curFile} PATH)
	LIST(FIND HPPDIRS ${curDir} CURDIR_FOUND)
	IF(${CURDIR_FOUND} EQUAL -1)
		LIST(APPEND HPPDIRS ${curDir})
		SET(HPPDIRS_STR "${HPPDIRS_STR} ${curDir}")
	ENDIF(${CURDIR_FOUND} EQUAL -1)
endforeach(curFile IN LISTS ALL_HPPS_LIST)

#remove directories also?
#configure_file(uninstall_mtl.sh.in uninstall_mtl.sh @ONLY)
configure_file(preinst.in preinst @ONLY)
configure_file(postrm.in postrm @ONLY)

SET(PACKAGE_TITLE "Matrix Template Library 4")
SET(PACKAGE_DESCRIPTION "MTL4 is a generic C++ library for linear algebra including advanced Krylov sub-space solvers. This package contains all headers of the MTL. The headers are installed in /usr/includes/boost. For some more information consult the README in /usr/share/mtl/ or read http://www.mtl4.org")
SET(CPACK_PACKAGE_NAME "MTL")
SET(CPACK_PACKAGE_VENDOR "Peter Gottschling & Andrew Lumsdaine")
SET(CPACK_PACKAGE_VERSION "4.0.${CurrentRevision}")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PACKAGE_TITLE} header files")
SET(CPACK_PACKAGE_CONTACT "mtl@osl.iu.edu")
SET(CPACK_RESOURCE_FILE_LICENSE "../../tools/license/license.mtl.txt") #-> to usr/share/mtl/
SET(CPACK_RESOURCE_FILE_README "../../README") #-> to usr/share/mtl
SET(CPACK_PACKAGING_INSTALL_PREFIX "/usr/")

SET(CPACK_RPM_POST_UNINSTALL_SCRIPT_FILE "${CMAKE_CURRENT_BINARY_DIR}/postrm")
SET(CPACK_RPM_PACKAGE_ARCHITECTURE "noarch")
#DONT use empty lines. everything after an empty line is not displayed
SET(CPACK_RPM_PACKAGE_DESCRIPTION "${PACKAGE_DESCRIPTION}")
#TODO: rpm-boost dependence
SET(CPACK_RPM_PACKAGE_REQUIRES "boost-devel >= 1.40")

#the debian way..
SET(CPACK_GENERATOR "DEB" "RPM" "TGZ" "TBZ2" "ZIP")
SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "all")
SET(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "${CMAKE_CURRENT_SOURCE_DIR}/preinst" "${CMAKE_CURRENT_SOURCE_DIR}/postrm")
SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libboost-dev (>= 1.40)")
SET(CPACK_DEBIAN_PACKAGE_DESCRIPTION "${PACKAGE_TITLE}\n ${PACKAGE_DESCRIPTION}")
#SET(CPACK_DEBIAN_PACKAGE_CONFLICTS "pmtl")
INCLUDE(CPack)
#################################################################################
## handy for learning CMAKE and debugging possible problems
######################## record all the variables CMAKE defined #################
file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/AllVariables.txt "")
get_cmake_property(VARS VARIABLES)
foreach (var ${VARS})
    file (APPEND ${CMAKE_CURRENT_BINARY_DIR}/AllVariables.txt "${var} \"${${var}}\"\n")
endforeach (var ${VARS})
