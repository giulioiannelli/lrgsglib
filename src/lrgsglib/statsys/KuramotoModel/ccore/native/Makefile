.PHONY: all build clean

all: build

build:
	python setup.py build_ext --inplace
	@mv -f _kuramoto_native*.so ../ 2>/dev/null || true

clean:
	rm -rf build *.so ../_kuramoto_native*.so
