OBJECTS = LinAl.o Driver_LinAl.o

MODULES = LinAl.mod

.PHONY: clean

output.txt: LinAl.exe
	./LinAl.exe > output.txt

LinAl.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) -o LinAl.exe

%.o: %.f90
	gfortran -c   $<

%.mod: %.f90
	gfortran -c   $<

clean:
	rm -f $(OBJECTS) $(MODULES) LinAl.exe
