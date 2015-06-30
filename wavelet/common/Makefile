include parallelDefs

COMMON = IO.h parseCommandLine.h parallel.h runTests.py sequence.h utils.h

WTCheck: WTCheck.C $(COMMON)
	$(PCC) $(PCFLAGS) $(PLFLAGS) -o WTCheck WTCheck.C

$(COMMON) :
	ln -s ../../common/$@ .

clean :
	rm -f WTCheck

cleansrc :
	make -s clean
	rm -f $(COMMON) *.pyc
