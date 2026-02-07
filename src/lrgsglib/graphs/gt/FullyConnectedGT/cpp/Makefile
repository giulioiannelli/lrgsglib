# FullyConnectedGT C++ extension Makefile
# Builds graph-tool C++ extension for complete graph generation

# Set default target before includes
.DEFAULT_GOAL := all

# Include shared graph-tool C++ build configuration
include ../../../../../../build/gt-cpp.mk

# Target
TARGETS := $(BUILD_DIR)/complete_graph.so

.PHONY: all clean

all: $(TARGETS)
	@echo "FullyConnectedGT C++ extension built."
	@ls -la $(BUILD_DIR)/*.so 2>/dev/null || true

$(BUILD_DIR)/complete_graph.so: complete_graph.cc complete_graph.hh
	@echo "Compiling complete_graph.so..."
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< -o $@

clean:
	@echo "Cleaning FullyConnectedGT C++ extension..."
	rm -f $(BUILD_DIR)/*.so $(BUILD_DIR)/*.o
