project (mtl4_tests)

# collect the source files
file (GLOB SOURCES "./*_test.cpp")
compile_all("true" "${SOURCES}")
compile_all("true" "vpt_test_1.cpp vpt_test_2.cpp")

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