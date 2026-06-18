all: build

build:
	python setup.py build_ext --inplace
	@# Move built .so up to parent (ccore/) for import from VoterModel
	@mv -f _voter_native*.so ../ 2>/dev/null || true

clean:
	rm -rf build *.so ../_voter_native*.so
