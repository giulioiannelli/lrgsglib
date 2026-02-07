# DualBarabasiAlbertGT C++ extension Makefile
# Builds graph-tool C++ extension for Dual Barabasi-Albert

# Set default target before includes
.DEFAULT_GOAL := all

# Include shared graph-tool C++ build configuration
include ../../../../../../build/gt-cpp.mk

# Target
TARGETS := $(BUILD_DIR)/dual_barabasi_albert.so

.PHONY: all clean

all: $(TARGETS)
	@echo "DualBarabasiAlbertGT C++ extension built."
	@ls -la $(BUILD_DIR)/*.so 2>/dev/null || true

$(BUILD_DIR)/dual_barabasi_albert.so: dual_barabasi_albert.cc dual_barabasi_albert.hh
	@echo "Compiling dual_barabasi_albert.so..."
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< -o $@

clean:
	@echo "Cleaning DualBarabasiAlbertGT C++ extension..."
	rm -f $(BUILD_DIR)/*.so $(BUILD_DIR)/*.o
