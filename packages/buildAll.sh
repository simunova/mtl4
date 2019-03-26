./cleanAll.sh
cd standard
rm -rf CMakeCache.txt
cmake .
make package
cd ..

cd examples
rm -rf CMakeCache.txt
cmake .
make package
cd ..

cd tests/
rm -rf CMakeCache.txt
cmake .
make package
cd ..

cd all/
rm -rf CMakeCache.txt
cmake .
make package
./createArchive.sh
cd ..

./createDebRepo.sh
./createRpmRepo.sh
./createDownloadTable.sh
