# ExtendedBarabasiAlbertGT C++ extension Makefile
# Builds graph-tool C++ extension for Extended Barabasi-Albert

# Set default target before includes
.DEFAULT_GOAL := all

# Include shared graph-tool C++ build configuration
include ../../../../../../build/gt-cpp.mk

# Target
TARGETS := $(BUILD_DIR)/extended_barabasi_albert.so

.PHONY: all clean

all: $(TARGETS)
	@echo "ExtendedBarabasiAlbertGT C++ extension built."
	@ls -la $(BUILD_DIR)/*.so 2>/dev/null || true

$(BUILD_DIR)/extended_barabasi_albert.so: extended_barabasi_albert.cc extended_barabasi_albert.hh
	@echo "Compiling extended_barabasi_albert.so..."
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< -o $@

clean:
	@echo "Cleaning ExtendedBarabasiAlbertGT C++ extension..."
	rm -f $(BUILD_DIR)/*.so $(BUILD_DIR)/*.o
