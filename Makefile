PYTHON=python3

all: build README.md

build: phony
	${PYTHON} -c "from setuptools import setup; setup()" build_ext --inplace

clean: phony
	rm -f gcad_ext.*.so gcad_ext/gcad_ext.cpp
	rm -rf __pycache__/ gcad/__pycache__/ tests/__pycache__/
	rm -rf build/ dist/ gcad.egg-info/
	rm -f README.md.tmp
	rm -f gcad/_version.py

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

README.md: phony
	sed '/## API/,$$d' README.md >README.md.tmp
	printf '## API\n\n' >>README.md.tmp
	$(PYTHON) mod2md.py >>README.md.tmp
	mv README.md.tmp README.md

phony:;
