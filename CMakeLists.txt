cmake_minimum_required(VERSION 3.20)
project(lrgsglib LANGUAGES C CXX)

# ===========================================================================
# Compiler setup
# ===========================================================================
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ===========================================================================
# External dependencies
# ===========================================================================
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)
find_package(pybind11 CONFIG REQUIRED)
find_package(OpenMP QUIET)

# ===========================================================================
# SFMT configuration
# ===========================================================================
set(SFMT_MEXP 19937 CACHE STRING "Mersenne exponent for SFMT")
add_compile_definitions(DSFMT_MEXP=${SFMT_MEXP})

# SIMD detection
include(CheckCCompilerFlag)
check_c_compiler_flag("-msse2" HAS_SSE2)
if(HAS_SSE2)
    add_compile_options(-msse2)
    add_compile_definitions(HAVE_SSE2)
endif()

# ===========================================================================
# Path variables
# ===========================================================================
set(CCORE      ${CMAKE_SOURCE_DIR}/src/lrgsglib/Ccore)
set(STATSYS    ${CCORE}/statsys)
set(SFMT_DIR   ${CCORE}/SFMT)
set(RBIM_SIMC  ${STATSYS}/RBIsingM/simulatorC)
set(RBIM_BASE  ${STATSYS}/RBIsingM/base)
set(RBIM_STORE ${STATSYS}/RBIsingM/storer)
set(SRW_LATT   ${STATSYS}/signedRw/Lattices)
set(VM_DIR     ${STATSYS}/voterM)
set(CP_DIR     ${STATSYS}/contactP)

# ===========================================================================
# Common source files
# ===========================================================================
set(COMMON_SOURCES
    ${CCORE}/LRGSG_utils.c
    ${CCORE}/sfmtrng.c
    ${SFMT_DIR}/SFMT.c
    ${STATSYS}/LRGSG_bindynsys.c
)

set(COMMON_INCLUDES
    ${CCORE}
    ${SFMT_DIR}
    ${STATSYS}
    ${RBIM_SIMC}
)

# Model-specific shared sources
set(SRC_RBIM  ${RBIM_SIMC}/LRGSG_rbim.c)
set(SRC_SA    ${RBIM_SIMC}/LRGSG_sa.c)
set(SRC_PT    ${RBIM_SIMC}/LRGSG_pt.c)
set(SRC_VM    ${VM_DIR}/LRGSG_vm.c)
set(SRC_CP    ${CP_DIR}/LRGSG_cp.c)

# Install destination for compiled binaries inside the Python package
set(BIN_INSTALL_DIR lrgsglib/Ccore/bin)

# ===========================================================================
# Helper functions
# ===========================================================================

# --- Ising simulator (standard: uses LRGSG_rbim) ---
function(add_ising_simulator NAME)
    add_executable(${NAME}
        ${RBIM_SIMC}/${NAME}.c
        ${COMMON_SOURCES}
        ${SRC_RBIM}
    )
    target_include_directories(${NAME} PRIVATE ${COMMON_INCLUDES})
    target_compile_options(${NAME} PRIVATE -O3 -g -Wall -Wextra)
    target_link_libraries(${NAME} PRIVATE m)
    if(OpenMP_FOUND)
        target_link_libraries(${NAME} PRIVATE OpenMP::OpenMP_C)
    endif()
    install(TARGETS ${NAME} DESTINATION ${BIN_INSTALL_DIR})
endfunction()

# --- Ising simulator with Simulated Annealing (extra: LRGSG_sa) ---
function(add_ising_sa_simulator NAME)
    add_executable(${NAME}
        ${RBIM_SIMC}/${NAME}.c
        ${COMMON_SOURCES}
        ${SRC_RBIM}
        ${SRC_SA}
    )
    target_include_directories(${NAME} PRIVATE ${COMMON_INCLUDES})
    target_compile_options(${NAME} PRIVATE -O3 -g -Wall -Wextra)
    target_link_libraries(${NAME} PRIVATE m)
    if(OpenMP_FOUND)
        target_link_libraries(${NAME} PRIVATE OpenMP::OpenMP_C)
    endif()
    install(TARGETS ${NAME} DESTINATION ${BIN_INSTALL_DIR})
endfunction()

# --- Ising simulator with Parallel Tempering (extra: LRGSG_pt) ---
function(add_ising_pt_simulator NAME)
    add_executable(${NAME}
        ${RBIM_SIMC}/${NAME}.c
        ${COMMON_SOURCES}
        ${SRC_RBIM}
        ${SRC_PT}
    )
    target_include_directories(${NAME} PRIVATE ${COMMON_INCLUDES})
    target_compile_options(${NAME} PRIVATE -O3 -g -Wall -Wextra)
    target_link_libraries(${NAME} PRIVATE m)
    if(OpenMP_FOUND)
        target_link_libraries(${NAME} PRIVATE OpenMP::OpenMP_C)
    endif()
    install(TARGETS ${NAME} DESTINATION ${BIN_INSTALL_DIR})
endfunction()

# --- Voter simulator ---
function(add_voter_simulator NAME)
    add_executable(${NAME}
        ${VM_DIR}/${NAME}.c
        ${COMMON_SOURCES}
        ${SRC_VM}
    )
    target_include_directories(${NAME} PRIVATE ${COMMON_INCLUDES})
    target_compile_options(${NAME} PRIVATE -O3 -g -Wall -Wextra)
    target_link_libraries(${NAME} PRIVATE m)
    if(OpenMP_FOUND)
        target_link_libraries(${NAME} PRIVATE OpenMP::OpenMP_C)
    endif()
    install(TARGETS ${NAME} DESTINATION ${BIN_INSTALL_DIR})
endfunction()

# --- Contact Process simulator ---
function(add_contact_simulator NAME)
    add_executable(${NAME}
        ${CP_DIR}/${NAME}.c
        ${COMMON_SOURCES}
        ${SRC_CP}
    )
    target_include_directories(${NAME} PRIVATE ${COMMON_INCLUDES})
    target_compile_options(${NAME} PRIVATE -O3 -g -Wall -Wextra)
    target_link_libraries(${NAME} PRIVATE m)
    if(OpenMP_FOUND)
        target_link_libraries(${NAME} PRIVATE OpenMP::OpenMP_C)
    endif()
    install(TARGETS ${NAME} DESTINATION ${BIN_INSTALL_DIR})
endfunction()

# ===========================================================================
# Ising simulators
# ===========================================================================

# Legacy Ising simulators (Algorithm number only)
add_ising_simulator(IsingSimulator0)
add_ising_simulator(IsingSimulator1)
add_ising_simulator(IsingSimulator2)
add_ising_simulator(IsingSimulator3)
add_ising_simulator(IsingSimulator4)
add_ising_simulator(IsingSimulator5)

# New naming convention simulators
# Algorithm 1: Glauber-Metropolis, b = E/M output
add_ising_simulator(IsingSimulator1b)

# Algorithm 3: Simulated Annealing, b = E/M output
add_ising_sa_simulator(IsingSimulator3b)

# Algorithm 4: Parallel Tempering, b = E/M output
add_ising_pt_simulator(IsingSimulator4b)

# ===========================================================================
# Voter Model simulators
# ===========================================================================
add_voter_simulator(VoterSimulator0)
add_voter_simulator(VoterSimulator1)

# ===========================================================================
# Contact Process simulators
# ===========================================================================

# Unified simulator (consolidates legacy 1, 1a-1g)
add_contact_simulator(ContactSimulator)

# Legacy variants (kept for backwards compatibility)
add_contact_simulator(ContactSimulator0)
add_contact_simulator(ContactSimulator1)
add_contact_simulator(ContactSimulator1a)
add_contact_simulator(ContactSimulator1b)
add_contact_simulator(ContactSimulator1c)
add_contact_simulator(ContactSimulator1d)
add_contact_simulator(ContactSimulator1e)
add_contact_simulator(ContactSimulator1f)
add_contact_simulator(ContactSimulator1g)

# ===========================================================================
# pybind11 Python extensions
# ===========================================================================

# Signed random walk on CSR adjacency matrix (already pybind11)
pybind11_add_module(_random_walk
    ${STATSYS}/random_walk.cpp
)
target_include_directories(_random_walk PRIVATE ${Python_NumPy_INCLUDE_DIRS})
target_compile_options(_random_walk PRIVATE -O3)
install(TARGETS _random_walk DESTINATION lrgsglib)

# Ising model (base) - pybind11 binding
pybind11_add_module(ising_model
    ${RBIM_BASE}/IsingModel.cpp
)
target_compile_options(ising_model PRIVATE -O3)
install(TARGETS ising_model DESTINATION lrgsglib/Ccore/statsys/RBIsingM/base)

# Ising model store - pybind11 binding
pybind11_add_module(ising_model_store
    ${RBIM_STORE}/main.cpp
    ${RBIM_STORE}/ising_model_store.cpp
)
target_include_directories(ising_model_store PRIVATE ${RBIM_STORE})
target_compile_options(ising_model_store PRIVATE -O3)
install(TARGETS ising_model_store DESTINATION lrgsglib/Ccore/statsys/RBIsingM/storer)

# Random walker on lattices - pybind11 binding
pybind11_add_module(random_walker
    ${SRW_LATT}/random_walker.cpp
)
target_compile_options(random_walker PRIVATE -O3)
install(TARGETS random_walker DESTINATION lrgsglib/Ccore/statsys/signedRw/Lattices)

# ===========================================================================
# graph-tool C++ extensions (NOT built by CMake)
# ===========================================================================
# Graph-tool extensions (barabasi_albert.so, triangular_lattice.so, etc.)
# use Boost.Python and link against libgraph_tool_core.so. They CANNOT be
# migrated to pybind11 because graph-tool's internal C++ API uses Boost.Python
# type converters.
#
# These are built separately via:
#   make cpp-make   (from the legacy Makefile)
# or per-graph:
#   cd src/lrgsglib/graphs/gt/<GraphType>/cpp && make
#
# Shared config: build/gt-cpp.mk
# See .agents/guides/gt-cpp-extensions.md for how to add new GT extensions.
