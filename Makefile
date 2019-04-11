option=

sector.so: sector/sector.pyx
	python setup.py build_ext -if $(option)

install:
	@rm -f sector/*.so
	@rm -f sector/utils.c
	@rm -f sector/sector.c
	@rm -f sector/sfh.c
	@rm -f sector/dust/dust.c
	@rm -f sector/dust/*.so
	pip install -e . --user

clean:
	@rm -f sector/*.so
	@rm -f sector/utils.c
	@rm -f sector/sector.c
	@rm -f sector/sfh.c
	@rm -f sector/dust/dust.c
	@rm -f sector/dust/*.so
	@rm -rf build
	@rm -rf *.egg-info
