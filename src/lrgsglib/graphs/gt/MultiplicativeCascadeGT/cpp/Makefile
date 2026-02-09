# MultiplicativeCascadeGT C++ extension Makefile
# Builds graph-tool C++ extension for vectorized cascade edge building

# Set default target before includes
.DEFAULT_GOAL := all

# Include shared graph-tool C++ build configuration
include ../../../../../../build/gt-cpp.mk

# Target
TARGETS := $(BUILD_DIR)/multiplicative_cascade.so

.PHONY: all clean

all: $(TARGETS)
	@echo "MultiplicativeCascadeGT C++ extension built."
	@ls -la $(BUILD_DIR)/*.so 2>/dev/null || true

$(BUILD_DIR)/multiplicative_cascade.so: multiplicative_cascade.cc multiplicative_cascade.hh
	@echo "Compiling multiplicative_cascade.so..."
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< -o $@

clean:
	@echo "Cleaning MultiplicativeCascadeGT C++ extension..."
	rm -f $(BUILD_DIR)/*.so $(BUILD_DIR)/*.o
