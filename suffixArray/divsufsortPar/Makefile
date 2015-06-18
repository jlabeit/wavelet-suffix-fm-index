include parallelDefs

# required files
SORT =  blockRadixSort.h transpose.h
OTHER = merge.h rangeMin.h
GLOBAL_REQUIRE = parallel.h sequence.h utils.h  $(SORT) $(OTHER)
LOCAL_REQUIRE =  divsufsort.C sssort.C trsort.C paralleltrsort.C 
OBJS =  divsufsort.o sssort.o trsort.o paralleltrsort.o 


include ../common/timeRequiredFiles
include MakeBench
