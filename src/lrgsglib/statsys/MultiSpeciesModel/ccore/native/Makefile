.PHONY: all build clean

all: build

build:
	python setup.py build_ext --inplace
	@# Move built .so up to parent (ccore/) for import from MultiSpeciesModel
	@mv -f _multispec_native*.so ../ 2>/dev/null || true

clean:
	rm -rf build *.so ../_multispec_native*.so
