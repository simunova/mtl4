#!/bin/bash

extractDir="/tmp/testMtlPackage"
curDir=$(pwd)
curRevision=$(svn info . | grep Revision | cut -d ' ' -f 2)

cd standard
rm -rf CMakeCache.txt
cmake . > /dev/null
if test $? -ne 0; then
  echo "could not configure standard package!"
  exit 1
fi
rm -rf *.tar.gz
cpack -G TGZ > /dev/null

if test $? -ne 0; then
  echo "could not create standard package!"
  exit 1
fi
package_mtl="${curDir}/standard/MTL-4.0.${curRevision}-Linux.tar.gz"
cd ..

cd tests
rm -rf CMakeCache.txt
cmake . > /dev/null

if test $? -ne 0; then
  echo "could not configure test package!"
  exit 1
fi
rm -rf *.tar.gz
cpack -G TGZ > /dev/null

if test $? -ne 0; then
  echo "could not create test package!"
  exit 1
fi
package_test="${curDir}/tests/MTL-tests-4.0.${curRevision}-Linux.tar.gz"
cd ..

rm -rf ${extractDir}
if test -e ${extractDir}; then
  echo "extract exists and could not delete it...\n"
  exit 1
fi
mkdir -p ${extractDir}

cd ${extractDir}
tar xzf ${package_mtl}
tar xzf ${package_test}

curMTLDir="${extractDir}/MTL-4.0.${curRevision}-Linux/usr/share/mtl/"
cd ${extractDir}/MTL-tests-4.0.${curRevision}-Linux/usr/share/mtl/test
mkdir build
cd build
cmake -DMTL_DIR=${curMTLDir} .. > logfile_configure 2> errfile_configure
if test $? -ne 0; then
  echo "configuration went wrong"
  exit 1
fi

echo "compiling tests, this may take some time"
make -j16 > logfile_make 2> errfile_make
if test $? -ne 0; then
  echo "compilation went wrong"
  exit 1
fi

echo "package creation with TGZ works"
