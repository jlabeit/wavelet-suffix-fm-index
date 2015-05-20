include parallelDefs

# required files
SORT =  blockRadixSort.h transpose.h
OTHER = merge.h rangeMin.h
GLOBAL_REQUIRE = parallel.h sequence.h utils.h  $(SORT) $(OTHER)
LOCAL_REQUIRE = src/sssort.C src/trsort.C src/utils.C src/divsufsort.C
OBJS = pks.o src/sssort.o src/trsort.o src/utils.o src/divsufsort.o


include ../common/timeRequiredFiles
include MakeBench
