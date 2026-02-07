# HolmeKimGT C++ extension Makefile
# Builds graph-tool C++ extension for Holme-Kim powerlaw cluster graph

# Set default target before includes
.DEFAULT_GOAL := all

# Include shared graph-tool C++ build configuration
include ../../../../../../build/gt-cpp.mk

# Target
TARGETS := $(BUILD_DIR)/holme_kim.so

.PHONY: all clean

all: $(TARGETS)
	@echo "HolmeKimGT C++ extension built."
	@ls -la $(BUILD_DIR)/*.so 2>/dev/null || true

$(BUILD_DIR)/holme_kim.so: holme_kim.cc holme_kim.hh
	@echo "Compiling holme_kim.so..."
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< -o $@

clean:
	@echo "Cleaning HolmeKimGT C++ extension..."
	rm -f $(BUILD_DIR)/*.so $(BUILD_DIR)/*.o
