all:
	ln -s EXPORT-15.1/DATA DATA
	(cd EXPORT-15.1/exec-v15.1; make)
	cp EXPORT-15.1/exec-v15.1/*_lu ../../bin/
clean:
	-rm *.o
