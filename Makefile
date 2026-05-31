PYTHON=python3

all: build

build: phony
	${PYTHON} -c "from setuptools import setup; setup()" build_ext --inplace

clean: phony
	rm -f gcad_ext.*.so gcad_ext/gcad.c
	rm -rf __pycache__/ gcad/__pycache__/ tests/__pycache__/
	rm -rf build/ dist/ gcad.egg-info/

dist: phony
	${PYTHON} -m build
	rm -rf gcad.egg-info/

last-commit-dist: phony
	mkdir -p dist
	tmp=$$(mktemp -d) && \
	    git clone . "$${tmp}/" && \
	    ${PYTHON} -m build "$${tmp}/" && \
	    cp -a $${tmp}/dist/* dist/ && \
	    rm -rf "$${tmp}/"

test: build phony
	${PYTHON} -m pytest -x --full-trace -m "not slow"

test-full: build phony
	${PYTHON} -m pytest -x --full-trace

phony:;
