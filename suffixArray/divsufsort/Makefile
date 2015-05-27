include parallelDefs

# required files
SORT =  blockRadixSort.h transpose.h
OTHER = merge.h rangeMin.h
GLOBAL_REQUIRE = parallel.h sequence.h utils.h  $(SORT) $(OTHER)
LOCAL_REQUIRE =  divsufsort.C sssort.C trsort.C utils.C paralleltrsort.C 
OBJS =  paralleltrsort.o divsufsort.o sssort.o trsort.o utils.o 


include ../common/timeRequiredFiles
include MakeBench
