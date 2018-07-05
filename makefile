option=

#LDFLAGS="-lrt" python setup.py build_ext -if $(option)
magcalc.so: magcalc.pyx sector_cext.c sector_cext.h
		python setup.py build_ext -if $(option)

clean:
	rm -rf build
	rm *.so
