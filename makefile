option=

#LDFLAGS="-lrt" python setup.py build_ext -if $(option)
mag_calc.so: mag_calc.pyx mag_calc_cext.c mag_calc_cext.h
		python setup.py build_ext -if $(option)

clean:
	rm -rf build
	rm mag_calc.so
