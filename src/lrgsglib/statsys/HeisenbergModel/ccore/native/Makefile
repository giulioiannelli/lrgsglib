.PHONY: all build clean

all: build

build:
	python setup.py build_ext --inplace
	@# Move built .so up to parent (ccore/) for import from HeisenbergModel
	@mv -f _heisenberg_native*.so ../ 2>/dev/null || true

clean:
	rm -rf build *.so ../_heisenberg_native*.so
