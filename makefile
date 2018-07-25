option=

#LDFLAGS="-lrt" python setup.py build_ext -if $(option)
magcalc.so: magcalc.pyx
		python setup.py build_ext -if $(option)

clean:
	rm -rf build
	rm *.so
