PYTHON=python3

all: build

build: phony
	${PYTHON} -c "from setuptools import setup; setup()" build_ext --inplace

clean: phony
	rm -f gcad_c.*.so gcad_c/gcad.c
	rm -rf __pycache__/ gcad/__pycache__/
	rm -rf build/ dist/ gcad.egg-info/

dist: phony
	${PYTHON} -m build
	rm -rf gcad.egg-info/

test: phony
	${PYTHON} -m pytest -m "not slow"

test-full: phony
	${PYTHON} -m pytest

phony:;
