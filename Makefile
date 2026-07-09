PYTHON?=python3
HEPWARE_COMMIT=39fa541a491f3ee4aaa485a305a334b35a6eb743 

all: build README.md

hepware/Makefile:
	git clone https://github.com/magv/hepware
	cd hepware && git checkout ${HEPWARE_COMMIT}

hepware/urlget: hepware/Makefile
	echo "import sys; from urllib.request import urlretrieve as get; get(sys.argv[2], filename=sys.argv[1])" >$@

print-hepware-id:
	@echo ${HEPWARE_COMMIT}

build-deps: hepware/Makefile hepware/urlget phony
	cd hepware && git checkout ${HEPWARE_COMMIT}
	+${MAKE} -C hepware flint.done mpfr.done gmp.done FETCH="${PYTHON} urlget"

build: build-deps phony
	${PYTHON} setup.py build_ext --inplace

test: build phony
	${PYTHON} -m pytest -x --full-trace -m "not slow"

test-full: build phony
	${PYTHON} -m pytest -x --full-trace

README.md: build phony
	sed '/## API/,$$d' README.md >README.md.tmp
	printf '## API\n\n' >>README.md.tmp
	${PYTHON} mod2md.py >>README.md.tmp
	mv README.md.tmp README.md

clean: phony
	rm -f gcad_ext.*.so gcad_ext/gcad_ext.cpp
	rm -rf __pycache__/ gcad/__pycache__/ tests/__pycache__/
	rm -rf build/ dist/ gcad.egg-info/
	rm -f README.md.tmp
	rm -f gcad/_version.py

sdist: phony
	${PYTHON} -m build --sdist --no-isolation

bdist: build phony
	${PYTHON} -m build --wheel --no-isolation

dist: sdist bdist

last-commit-sdist: phony
	${PYTHON} -m build --sdist

phony:;
