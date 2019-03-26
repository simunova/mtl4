find_program(VAMPIR_CXX vtc++)
option(ENABLE_VAMPIR "enable or disable vampir trace" OFF)
if(VAMPIR_CXX)
	#function to compute the compile flags, libraries, linker directories and linker flags used by the vampir compiler
	#cflags: name of the global variable which will contain the compilerflags
	#libs: name of the global variable which will contain the libraries
	#ldirs: name of the global variable which will contain the linker directories
	#lflags: name of the global variable which will contain the linker flags
	#addArgs: additional arguments passed to the vampir compiler
	function(GetVampirData cflags libs ldirs lflags _addArgs)
		set(addArgs "${${_addArgs}}")
		#request only compile arguments
		execute_process(COMMAND ${VAMPIR_CXX} "-vt:show" "-vt:inst" "manual" ${addArgs} "-c" RESULT_VARIABLE VT_RES OUTPUT_VARIABLE VT_OUT)
		#message("vt_res: ${VT_RES}")
		#message("vt_out: ${VT_OUT}")
		string(REPLACE " " ";" asList ${VT_OUT} )
		list(LENGTH asList length)
		#message("asList: ${asList}")
		#message("length: ${length}")
		list(GET asList 0 VT_COMPILER)
		#remove the compiler
		list(REMOVE_AT asList 0)
		#add each flag except the last to the vampir compile flags
		while(length GREATER 1)
			list(GET asList 0 CURFLAG)
			#save the results in internal variable
			list(APPEND _VT_COMPILE_FLAGS ${CURFLAG})
			list(REMOVE_AT asList 0)
			list(LENGTH asList length)
		endwhile()
		#message("asList: ${asList}")
		#request all other flags
		execute_process(COMMAND ${VAMPIR_CXX} "-vt:show" "-vt:inst" "manual" ${addArgs} RESULT_VARIABLE VT_RES OUTPUT_VARIABLE VT_OUT)
		string(REPLACE " " ";" asList ${VT_OUT} )
		list(LENGTH asList length)
		#remove the compiler
		list(REMOVE_AT asList 0)
		while(length GREATER 1)
			list(GET asList 0 CURFLAG)
			#if the current flag is not a compile flag, add it
			list(FIND _VT_COMPILE_FLAGS ${CURFLAG} ISCOMPILEFLAG)
			#message("iscompileflag: ${ISCOMPILEFLAG}")
			if(ISCOMPILEFLAG EQUAL -1)
				string(SUBSTRING "${CURFLAG}" 0 2 FLAGBEG)

				string(REGEX MATCHALL "-L([^\" ]+|\"[^\"]+\")" CURLINKDIR "${CURFLAG}")
				set(CURLINKDIR ${CMAKE_MATCH_1})
				if(CURLINKDIR)
					list(APPEND _VT_LINKER_DIRECTORIES ${CURLINKDIR})
				endif()
				#save the results in internal variable
				#list(APPEND _VT_LINKER_DIRECTORIES ${CURLINKDIR})

				string(REGEX MATCHALL "-l([^\" ]+|\"[^\"]+\")" CURLIBRARY "${CURFLAG}")
				set(CURLIBRARY ${CMAKE_MATCH_1})
				if(CURLIBRARY)
					list(APPEND _VT_LIBRARIES ${CURLIBRARY})
				endif()
				#save the results in internal variable
				#list(APPEND _VT_LIBRARIES ${CURLIBRARY})
				#set(_VT_LINK_FLAGS "${_VT_LINK_FLAGS} ${CURFLAG}")
				if(NOT (CURLIBRARY OR CURLINKDIR)) 
					set(VT_LINK_FLAGS "${VT_LINK_FLAGS} ${CURFLAG}")
				endif()

			endif()
			list(REMOVE_AT asList 0)
			list(LENGTH asList length)
		endwhile()
		#propagate the values to parent scope
		set(${cflags} ${_VT_COMPILE_FLAGS} PARENT_SCOPE)
		set(${libs} ${_VT_LIBRARIES} PARENT_SCOPE)
		set(${ldirs} ${_VT_LINKER_DIRECTORIES} PARENT_SCOPE)
		set(${lflags} "${_VT_LINK_FLAGS}" PARENT_SCOPE)
	endfunction()
	set(VAMPIR_FOUND TRUE)
	#request only compile arguments
#	execute_process(COMMAND ${VAMPIR_CXX} "-vt:show" "-vt:mt" "-vt:noopari" "-vt:inst" "manual" "-c" RESULT_VARIABLE VT_RES OUTPUT_VARIABLE VT_OUT)
#	string(REPLACE " " ";" asList ${VT_OUT} )
#	list(LENGTH asList length)
#	list(GET asList 0 VT_COMPILER)
	#remove the compiler
#	list(REMOVE_AT asList 0)
	#add each flag except the last to the vampir compile flags
#	set(VT_COMPILE_FLAGS )
#	while(length GREATER 1)
#		list(GET asList 0 CURFLAG)
#		list(APPEND VT_COMPILE_FLAGS ${CURFLAG})
#      		list(REMOVE_AT asList 0)
#		list(LENGTH asList length)
#	endwhile()
	#request all other flags
#	execute_process(COMMAND ${VAMPIR_CXX} "-vt:show" "-vt:mt" "-vt:noopari" "-vt:inst" "manual" RESULT_VARIABLE VT_RES OUTPUT_VARIABLE VT_OUT)
#	string(REPLACE " " ";" asList ${VT_OUT} )
#	list(LENGTH asList length)
	#remove the compiler
#	list(REMOVE_AT asList 0)
#	set(VT_LIBRARIES )
#	set(VT_LINKER_DIRECTORIES )
#	set(VT_LINK_FLAGS "")
#	while(length GREATER 1)
#		list(GET asList 0 CURFLAG)
		#if the current flag is not a compile flag, add it
#		list(FIND VT_COMPILE_FLAGS ${CURFLAG} ISCOMPILEFLAG)
#		if(ISCOMPILEFLAG EQUAL -1)
#			string(SUBSTRING "${CURFLAG}" 0 2 FLAGBEG)

#			string(REGEX MATCHALL "-L([^\" ]+|\"[^\"]+\")" CURLINKDIR "${CURFLAG}")
#			set(CURLINKDIR ${CMAKE_MATCH_1})
#			if(CURLINKDIR)
#				list(APPEND VT_LINKER_DIRECTORIES ${CURLINKDIR})
#			endif()

#			string(REGEX MATCHALL "-l([^\" ]+|\"[^\"]+\")" CURLIBRARY "${CURFLAG}")
#			set(CURLIBRARY ${CMAKE_MATCH_1})
#			if(CURLIBRARY)
#				list(APPEND VT_LIBRARIES_REL ${CURLIBRARY})
#			endif()

#			if(NOT (CURLIBRARY OR CURLINKDIR)) 
#				set(VT_LINK_FLAGS "${VT_LINK_FLAGS} ${CURFLAG}")
#			endif()
#		endif()
#      		list(REMOVE_AT asList 0)
#		list(LENGTH asList length)
#	endwhile()

set(_ADD_VT_FLAGS "")
if(ENABLE_OPENMP)
	set(_ADD_VT_FLAGS "${VT_FLGAS};-vt:mt;-vt:noopari")
endif(ENABLE_OPENMP)
GetVampirData(VT_COMPILE_FLAGS VT_LIBRARIES_REL VT_LINKER_DIRECTORIES VT_LINK_FLAGS _ADD_VT_FLAGS)

#look for fullpath of vampir libraries
foreach(curLib IN LISTS VT_LIBRARIES_REL)
	find_library(CURRELLIB_${curLib} ${curLib} PATHS ${VT_LINKER_DIRECTORIES} NO_DEFAULT_PATH)
	if(CURRELLIB_${curLib})
	  list(APPEND VT_LIBRARIES ${CURRELLIB_${curLib}})
	else()
	  list(APPEND VT_LIBRARIES ${curLib})
	endif()
endforeach()
#message("vt_compiler: ${VT_COMPILER}")
#message("vt_compiler flags: ${VT_COMPILE_FLAGS}")
#message("vt_linker flags: ${VT_LINK_FLAGS}")
#message("vt_linker directories: ${VT_LINKER_DIRECTORIES}")
#message("vt_linker libraries: ${VT_LIBRARIES}")
else(VAMPIR_CXX)
	set(VAMPIR_FOUND FALSE)
endif()
