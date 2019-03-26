project (mtl4_timing)

# collect the source files
file (GLOB SOURCES [a-zA-Z]*.cpp) # filenames starting with '#' screw build
compile_all("false" ${SOURCES})
