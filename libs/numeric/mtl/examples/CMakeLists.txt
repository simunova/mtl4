project (mtl4_examples)

# collect the source files
file (GLOB SOURCES *.cpp)
compile_all("true" ${SOURCES})

# Copy each matrix market file to build dir
file (GLOB MTX_FILES matrix_market/*.mtx)
foreach(FILE ${MTX_FILES})
  get_filename_component(FILENAME ${FILE} NAME)
  
  add_custom_target(
    CopyTestFile_${FILENAME} ALL
    COMMAND cmake -E copy_if_different ${FILE} ${CMAKE_CURRENT_BINARY_DIR}/matrix_market/${FILENAME}
    COMMENT "Copying ${FILENAME}"
  )
endforeach(FILE)