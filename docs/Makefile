# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Custom targets for development

.PHONY: livehtml
livehtml:
	@echo "Starting live documentation server..."
	sphinx-autobuild "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O) \
		--watch ../src/lrgsglib \
		--open-browser \
		--delay 1

.PHONY: doctest
doctest:
	@echo "Running doctests..."
	@$(SPHINXBUILD) -b doctest "$(SOURCEDIR)" "$(BUILDDIR)/doctest" $(SPHINXOPTS) $(O)

.PHONY: clean-all
clean-all:
	@echo "Removing build directory and generated files..."
	rm -rf "$(BUILDDIR)"
	rm -rf "$(SOURCEDIR)/api/generated"
	find "$(SOURCEDIR)" -name "*.pyc" -delete
	find "$(SOURCEDIR)" -name "__pycache__" -delete

.PHONY: html-strict
html-strict:
	@echo "Building HTML documentation with warnings as errors..."
	@$(SPHINXBUILD) -W -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O)

.PHONY: linkcheck
linkcheck:
	@echo "Checking external links..."
	@$(SPHINXBUILD) -b linkcheck "$(SOURCEDIR)" "$(BUILDDIR)/linkcheck" $(SPHINXOPTS) $(O)

.PHONY: coverage
coverage:
	@echo "Checking documentation coverage..."
	@$(SPHINXBUILD) -b coverage "$(SOURCEDIR)" "$(BUILDDIR)/coverage" $(SPHINXOPTS) $(O)
	@cat "$(BUILDDIR)/coverage/python.txt"
