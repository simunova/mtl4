project(mtl-package-all)
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

SET(PACKAGE_TITLE "Matrix Template Library 4 - all")
SET(PACKAGE_DESCRIPTION "This is only a dummy package, which depends on MTL, MTL-examples and MTL-tests")
SET(CPACK_PACKAGE_NAME "MTL-all")
SET(CPACK_PACKAGE_VENDOR "Peter Gottschling & Andrew Lumsdaine")
SET(CPACK_PACKAGE_VERSION "4.0.${CurrentRevision}")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PACKAGE_TITLE}")
SET(CPACK_PACKAGE_CONTACT "mtl@osl.iu.edu")
SET(CPACK_RESOURCE_FILE_LICENSE "../../tools/license/license.mtl.txt") #-> to usr/share/mtl/
SET(CPACK_RESOURCE_FILE_README "../../README") #-> to usr/share/mtl
SET(CPACK_PACKAGING_INSTALL_PREFIX "/usr/")

SET(CPACK_RPM_PACKAGE_ARCHITECTURE "noarch")
SET(CPACK_RPM_PACKAGE_DESCRIPTION "${PACKAGE_DESCRIPTION}")
#TODO: rpm-boost dependence
SET(CPACK_RPM_PACKAGE_REQUIRES "mtl==${CPACK_PACKAGE_VERSION} , mtl-examples==${CPACK_PACKAGE_VERSION} , mtl-tests==${CPACK_PACKAGE_VERSION} ")

#the debian way..
SET(CPACK_GENERATOR "DEB" "RPM" "TGZ" "TBZ2")
SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "all")
SET(CPACK_DEBIAN_PACKAGE_DEPENDS "mtl (==${CPACK_PACKAGE_VERSION}), mtl-examples (==${CPACK_PACKAGE_VERSION}), mtl-tests (==${CPACK_PACKAGE_VERSION})")
#DONT use empty lines. everything after an empty line is not displayed
SET(CPACK_DEBIAN_PACKAGE_DESCRIPTION "${PACKAGE_TITLE}\n ${PACKAGE_DESCRIPTION}")
#SET(CPACK_DEBIAN_PACKAGE_CONFLICTS "pmtl")
INCLUDE(CPack)
