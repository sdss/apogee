all: 2020
15.1:
	ln -fs EXPORT-15.1/DATA DATA
	(cd EXPORT-15.1/exec-v15.1; make)
	cp EXPORT-15.1/exec-v15.1/*_lu ../../bin/
2020:
	ln -fs Turbospectrum2020/DATA DATA
	(cd Turbospectrum2020/exec; make)
	cp Turbospectrum2020/exec/*_lu ../../bin/
clean:
	-rm *.o
