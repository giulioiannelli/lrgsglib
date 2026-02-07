# BarabasiAlbertGT C++ extension Makefile
# Builds graph-tool C++ extension for Barabasi-Albert preferential attachment

# Set default target before includes
.DEFAULT_GOAL := all

# Include shared graph-tool C++ build configuration
include ../../../../../../build/gt-cpp.mk

# Target
TARGETS := $(BUILD_DIR)/barabasi_albert.so

.PHONY: all clean

all: $(TARGETS)
	@echo "BarabasiAlbertGT C++ extension built."
	@ls -la $(BUILD_DIR)/*.so 2>/dev/null || true

$(BUILD_DIR)/barabasi_albert.so: barabasi_albert.cc barabasi_albert.hh
	@echo "Compiling barabasi_albert.so..."
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< -o $@

clean:
	@echo "Cleaning BarabasiAlbertGT C++ extension..."
	rm -f $(BUILD_DIR)/*.so $(BUILD_DIR)/*.o
