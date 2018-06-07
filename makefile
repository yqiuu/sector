option=

#LDFLAGS="-lrt" python setup.py build_ext -if $(option)
magcalc.so: magcalc.pyx mag_calc_cext.c mag_calc_cext.h
		python setup.py build_ext -if $(option)

clean:
	rm -rf build
	rm *.so
