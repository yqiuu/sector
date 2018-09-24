option=

sector.so: sector/sector.pyx
	python setup.py build_ext -if $(option)

install:
	@rm -f secotr/*.so
	@rm -f sector/sector.c
	@rm -rf build
	@rm -rf *.egg-info
	pip install -e . --user

clean:
	@rm -f secotr/*.so
	@rm -f sector/sector.c
	@rm -rf build
	@rm -rf *.egg-info
