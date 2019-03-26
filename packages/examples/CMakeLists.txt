project(MTL_examples)

cmake_minimum_required(VERSION 2.6)
#this subproject prepares the MTL_examples package for installing in /usr/share/mtl/examples

#let cmake look for the svn-version
find_package(Subversion)
if(Subversion_FOUND)
	Subversion_WC_INFO(${CMAKE_CURRENT_BINARY_DIR} mtlSubinfo)
	set(CurrentRevision ${mtlSubinfo_WC_REVISION})
#	message("current revision: ${mtlSubinfo_WC_REVISION}")
else(Subversion_FOUND)
	set(CurrentRevision "0")
endif(Subversion_FOUND)

INSTALL(DIRECTORY ../../libs/numeric/mtl/examples
	DESTINATION .
	FILES_MATCHING PATTERN "*.cpp"
	PATTERN "*.hpp"
	PATTERN "*.mtx"
	)

INSTALL(FILES CMakeLists_examples.txt
	DESTINATION examples/
	RENAME CMakeLists.txt
       )

#the documentation, build with doxygen before
#add_custom_command(OUTPUT documentation.tar.gz COMMAND tar -czf documentation.tar.gz ../libs/numeric/mtl/doc)
add_custom_target(documentation ALL ./create_doc.sh)
INSTALL(FILES documentation.tar.gz 
	DESTINATION .
       )
SET(PACKAGE_TITLE "Matrix Template Library 4 - examples")
SET(PACKAGE_DESCRIPTION "This package contains examples and the documentation for the MTL. They are stored in /usr/share/mtl/examples. You can simply check the examples by copying the directory /usr/share/mtl/examples. If you are familiar with cmake, you can use the file CMakeLists.txt from the same directory.")
SET(CPACK_PACKAGE_NAME "MTL-examples")
SET(CPACK_PACKAGE_VENDOR "Peter Gottschling & Andrew Lumsdaine")
SET(CPACK_PACKAGE_VERSION "4.0.${CurrentRevision}")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PACKAGE_TITLE}")
SET(CPACK_PACKAGE_CONTACT "mtl@osl.iu.edu")
SET(CPACK_RESOURCE_FILE_LICENSE "../../tools/license/license.mtl.txt")
SET(CPACK_RESOURCE_FILE_README "../../README")
SET(CPACK_PACKAGING_INSTALL_PREFIX "/usr/share/mtl/")

SET(CPACK_RPM_POST_UNINSTALL_SCRIPT_FILE "${CMAKE_CURRENT_BINARY_DIR}/postrm")
SET(CPACK_RPM_POST_INSTALL_SCRIPT_FILE "${CMAKE_CURRENT_BINARY_DIR}/postinst")
SET(CPACK_RPM_PACKAGE_ARCHITECTURE "noarch")
SET(CPACK_RPM_PACKAGE_DESCRIPTION "${PACKAGE_DESCRIPTION}")
SET(CPACK_RPM_PACKAGE_REQUIRES "mtl")

#the debian way..
SET(CPACK_GENERATOR "DEB" "RPM" "TGZ" "TBZ2" "ZIP")
SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "all")
SET(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "${CMAKE_CURRENT_SOURCE_DIR}/preinst" "${CMAKE_CURRENT_SOURCE_DIR}/postinst" "${CMAKE_CURRENT_SOURCE_DIR}/postrm")
SET(CPACK_DEBIAN_PACKAGE_DEPENDS "mtl")
#DONT use empty lines. everything after an empty line is not displayed
SET(CPACK_DEBIAN_PACKAGE_DESCRIPTION "${PACKAGE_TITLE}\n ${PACKAGE_DESCRIPTION}")
INCLUDE(CPack)
