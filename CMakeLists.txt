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
# Path variables  (folder-per-class layout under statsys/)
# ===========================================================================
set(LIB_STATSYS  ${CMAKE_SOURCE_DIR}/src/lrgsglib/statsys)
set(CCORE        ${LIB_STATSYS}/_ccore)
set(SFMT_DIR     ${CCORE}/SFMT)

# Per-model C source directories
set(ISING_CCORE  ${LIB_STATSYS}/IsingDynamics/ccore)
set(RBIM_BASE    ${ISING_CCORE}/base)
set(RBIM_STORE   ${ISING_CCORE}/storer)
set(CP_DIR       ${LIB_STATSYS}/ContactProcess/ccore)
set(VM_DIR       ${LIB_STATSYS}/VoterModel/ccore)
set(SRW_CCORE    ${LIB_STATSYS}/SignedRW/ccore)
set(SRW_LATT     ${SRW_CCORE}/Lattices)
set(KUR_DIR      ${LIB_STATSYS}/KuramotoModel/ccore)
set(RD_DIR       ${LIB_STATSYS}/ReactionDiffusionModel/ccore)
set(CODE_DIR     ${LIB_STATSYS}/CoupledODEModel/ccore)
set(POTTS_DIR    ${LIB_STATSYS}/PottsModel/ccore)
set(XY_DIR       ${LIB_STATSYS}/XYModel/ccore)
set(HBERG_DIR    ${LIB_STATSYS}/HeisenbergModel/ccore)
set(MSPEC_DIR    ${LIB_STATSYS}/MultiSpeciesModel/ccore)

# ===========================================================================
# Common source files
# ===========================================================================
set(COMMON_SOURCES
    ${CCORE}/LRGSG_utils.c
    ${CCORE}/sfmtrng.c
    ${SFMT_DIR}/SFMT.c
    ${CCORE}/LRGSG_bindynsys.c
)

set(COMMON_INCLUDES
    ${CCORE}
    ${SFMT_DIR}
    ${ISING_CCORE}
)

# Model-specific shared sources
set(SRC_RBIM  ${ISING_CCORE}/LRGSG_rbim.c)
set(SRC_SA    ${ISING_CCORE}/LRGSG_sa.c)
set(SRC_PT    ${ISING_CCORE}/LRGSG_pt.c)
set(SRC_VM    ${VM_DIR}/LRGSG_vm.c)
set(SRC_CP    ${CP_DIR}/LRGSG_cp.c)

# Continuous / vector dynamics shared sources
set(SRC_CONTDYNSYS  ${CCORE}/LRGSG_contdynsys.c)
set(SRC_VECDYNSYS   ${CCORE}/LRGSG_vecdynsys.c)
set(SRC_KURAMOTO    ${KUR_DIR}/LRGSG_kuramoto.c)
set(SRC_RD          ${RD_DIR}/LRGSG_rd.c)
set(SRC_CODE        ${CODE_DIR}/LRGSG_code.c)
set(SRC_POTTS       ${POTTS_DIR}/LRGSG_potts.c)
set(SRC_XY          ${XY_DIR}/LRGSG_xy.c)
set(SRC_HBERG       ${HBERG_DIR}/LRGSG_heisenberg.c)
set(SRC_MSPEC       ${MSPEC_DIR}/LRGSG_multispec.c)

# Per-model install destinations
set(ISING_BIN_DIR  lrgsglib/statsys/IsingDynamics/ccore/bin)
set(CP_BIN_DIR     lrgsglib/statsys/ContactProcess/ccore/bin)
set(VM_BIN_DIR     lrgsglib/statsys/VoterModel/ccore/bin)
set(KUR_BIN_DIR    lrgsglib/statsys/KuramotoModel/ccore/bin)
set(RD_BIN_DIR     lrgsglib/statsys/ReactionDiffusionModel/ccore/bin)
set(CODE_BIN_DIR   lrgsglib/statsys/CoupledODEModel/ccore/bin)
set(POTTS_BIN_DIR  lrgsglib/statsys/PottsModel/ccore/bin)
set(XY_BIN_DIR     lrgsglib/statsys/XYModel/ccore/bin)
set(HBERG_BIN_DIR  lrgsglib/statsys/HeisenbergModel/ccore/bin)
set(MSPEC_BIN_DIR  lrgsglib/statsys/MultiSpeciesModel/ccore/bin)

# ===========================================================================
# Helper functions
# ===========================================================================

# --- Ising simulator (standard: uses LRGSG_rbim) ---
function(add_ising_simulator NAME)
    add_executable(${NAME}
        ${ISING_CCORE}/${NAME}.c
        ${COMMON_SOURCES}
        ${SRC_RBIM}
    )
    target_include_directories(${NAME} PRIVATE ${COMMON_INCLUDES})
    target_compile_options(${NAME} PRIVATE -O3 -g -Wall -Wextra)
    target_link_libraries(${NAME} PRIVATE m)
    if(OpenMP_FOUND)
        target_link_libraries(${NAME} PRIVATE OpenMP::OpenMP_C)
    endif()
    install(TARGETS ${NAME} DESTINATION ${ISING_BIN_DIR})
endfunction()

# --- Ising simulator with Simulated Annealing (extra: LRGSG_sa) ---
function(add_ising_sa_simulator NAME)
    add_executable(${NAME}
        ${ISING_CCORE}/${NAME}.c
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
    install(TARGETS ${NAME} DESTINATION ${ISING_BIN_DIR})
endfunction()

# --- Ising simulator with Parallel Tempering (extra: LRGSG_pt) ---
function(add_ising_pt_simulator NAME)
    add_executable(${NAME}
        ${ISING_CCORE}/${NAME}.c
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
    install(TARGETS ${NAME} DESTINATION ${ISING_BIN_DIR})
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
    install(TARGETS ${NAME} DESTINATION ${VM_BIN_DIR})
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
    install(TARGETS ${NAME} DESTINATION ${CP_BIN_DIR})
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
# Kuramoto Model simulators
# ===========================================================================
function(add_kuramoto_simulator NAME)
    add_executable(${NAME}
        ${KUR_DIR}/${NAME}.c
        ${COMMON_SOURCES}
        ${SRC_CONTDYNSYS}
        ${SRC_KURAMOTO}
    )
    target_include_directories(${NAME} PRIVATE ${COMMON_INCLUDES} ${KUR_DIR})
    target_compile_options(${NAME} PRIVATE -O3 -g -Wall -Wextra)
    target_link_libraries(${NAME} PRIVATE m)
    if(OpenMP_FOUND)
        target_link_libraries(${NAME} PRIVATE OpenMP::OpenMP_C)
    endif()
    install(TARGETS ${NAME} DESTINATION ${KUR_BIN_DIR})
endfunction()

add_kuramoto_simulator(KuramotoSimulator0)

# ===========================================================================
# Reaction-Diffusion simulators
# ===========================================================================
function(add_rd_simulator NAME)
    add_executable(${NAME}
        ${RD_DIR}/${NAME}.c
        ${COMMON_SOURCES}
        ${SRC_CONTDYNSYS}
        ${SRC_RD}
    )
    target_include_directories(${NAME} PRIVATE ${COMMON_INCLUDES} ${RD_DIR})
    target_compile_options(${NAME} PRIVATE -O3 -g -Wall -Wextra)
    target_link_libraries(${NAME} PRIVATE m)
    if(OpenMP_FOUND)
        target_link_libraries(${NAME} PRIVATE OpenMP::OpenMP_C)
    endif()
    install(TARGETS ${NAME} DESTINATION ${RD_BIN_DIR})
endfunction()

add_rd_simulator(ReactionDiffusionSimulator0)

# ===========================================================================
# Coupled ODE simulators
# ===========================================================================
function(add_code_simulator NAME)
    add_executable(${NAME}
        ${CODE_DIR}/${NAME}.c
        ${COMMON_SOURCES}
        ${SRC_CONTDYNSYS}
        ${SRC_CODE}
    )
    target_include_directories(${NAME} PRIVATE ${COMMON_INCLUDES} ${CODE_DIR})
    target_compile_options(${NAME} PRIVATE -O3 -g -Wall -Wextra)
    target_link_libraries(${NAME} PRIVATE m)
    if(OpenMP_FOUND)
        target_link_libraries(${NAME} PRIVATE OpenMP::OpenMP_C)
    endif()
    install(TARGETS ${NAME} DESTINATION ${CODE_BIN_DIR})
endfunction()

add_code_simulator(CoupledODESimulator0)

# ===========================================================================
# Potts Model simulators
# ===========================================================================
function(add_potts_simulator NAME)
    add_executable(${NAME}
        ${POTTS_DIR}/${NAME}.c
        ${COMMON_SOURCES}
        ${SRC_VECDYNSYS}
        ${SRC_POTTS}
    )
    target_include_directories(${NAME} PRIVATE ${COMMON_INCLUDES} ${POTTS_DIR})
    target_compile_options(${NAME} PRIVATE -O3 -g -Wall -Wextra)
    target_link_libraries(${NAME} PRIVATE m)
    if(OpenMP_FOUND)
        target_link_libraries(${NAME} PRIVATE OpenMP::OpenMP_C)
    endif()
    install(TARGETS ${NAME} DESTINATION ${POTTS_BIN_DIR})
endfunction()

add_potts_simulator(PottsSimulator0)

# ===========================================================================
# XY Model simulators
# ===========================================================================
function(add_xy_simulator NAME)
    add_executable(${NAME}
        ${XY_DIR}/${NAME}.c
        ${COMMON_SOURCES}
        ${SRC_CONTDYNSYS}
        ${SRC_VECDYNSYS}
        ${SRC_XY}
    )
    target_include_directories(${NAME} PRIVATE ${COMMON_INCLUDES} ${XY_DIR})
    target_compile_options(${NAME} PRIVATE -O3 -g -Wall -Wextra)
    target_link_libraries(${NAME} PRIVATE m)
    if(OpenMP_FOUND)
        target_link_libraries(${NAME} PRIVATE OpenMP::OpenMP_C)
    endif()
    install(TARGETS ${NAME} DESTINATION ${XY_BIN_DIR})
endfunction()

add_xy_simulator(XYSimulator0)

# ===========================================================================
# Heisenberg Model simulators
# ===========================================================================
function(add_heisenberg_simulator NAME)
    add_executable(${NAME}
        ${HBERG_DIR}/${NAME}.c
        ${COMMON_SOURCES}
        ${SRC_CONTDYNSYS}
        ${SRC_VECDYNSYS}
        ${SRC_HBERG}
    )
    target_include_directories(${NAME} PRIVATE ${COMMON_INCLUDES} ${HBERG_DIR})
    target_compile_options(${NAME} PRIVATE -O3 -g -Wall -Wextra)
    target_link_libraries(${NAME} PRIVATE m)
    if(OpenMP_FOUND)
        target_link_libraries(${NAME} PRIVATE OpenMP::OpenMP_C)
    endif()
    install(TARGETS ${NAME} DESTINATION ${HBERG_BIN_DIR})
endfunction()

add_heisenberg_simulator(HeisenbergSimulator0)

# ===========================================================================
# MultiSpecies simulators
# ===========================================================================
function(add_multispec_simulator NAME)
    add_executable(${NAME}
        ${MSPEC_DIR}/${NAME}.c
        ${COMMON_SOURCES}
        ${SRC_VECDYNSYS}
        ${SRC_MSPEC}
    )
    target_include_directories(${NAME} PRIVATE ${COMMON_INCLUDES} ${MSPEC_DIR})
    target_compile_options(${NAME} PRIVATE -O3 -g -Wall -Wextra)
    target_link_libraries(${NAME} PRIVATE m)
    if(OpenMP_FOUND)
        target_link_libraries(${NAME} PRIVATE OpenMP::OpenMP_C)
    endif()
    install(TARGETS ${NAME} DESTINATION ${MSPEC_BIN_DIR})
endfunction()

add_multispec_simulator(MultiSpeciesSimulator0)

# ===========================================================================
# pybind11 Python extensions
# ===========================================================================

# Signed random walk on CSR adjacency matrix (already pybind11)
pybind11_add_module(_random_walk
    ${SRW_CCORE}/random_walk.cpp
)
target_include_directories(_random_walk PRIVATE ${Python_NumPy_INCLUDE_DIRS})
target_compile_options(_random_walk PRIVATE -O3)
install(TARGETS _random_walk DESTINATION lrgsglib)

# Ising model (base) - pybind11 binding
pybind11_add_module(ising_model
    ${RBIM_BASE}/IsingModel.cpp
)
target_compile_options(ising_model PRIVATE -O3)
install(TARGETS ising_model DESTINATION lrgsglib/statsys/IsingDynamics/ccore/base)

# Ising model store - pybind11 binding
pybind11_add_module(ising_model_store
    ${RBIM_STORE}/main.cpp
    ${RBIM_STORE}/ising_model_store.cpp
)
target_include_directories(ising_model_store PRIVATE ${RBIM_STORE})
target_compile_options(ising_model_store PRIVATE -O3)
install(TARGETS ising_model_store DESTINATION lrgsglib/statsys/IsingDynamics/ccore/storer)

# Ising native kernels - pybind11 binding (Metropolis, SA, PT via numpy)
pybind11_add_module(_ising_native
    ${ISING_CCORE}/ising_native.cpp
    ${SRC_RBIM}
    ${SRC_SA}
    ${SRC_PT}
    ${COMMON_SOURCES}
)
target_include_directories(_ising_native PRIVATE
    ${COMMON_INCLUDES}
    ${Python_NumPy_INCLUDE_DIRS}
)
target_compile_options(_ising_native PRIVATE -O3)
install(TARGETS _ising_native DESTINATION lrgsglib/statsys/IsingDynamics/ccore)

# Random walker on lattices - pybind11 binding
pybind11_add_module(random_walker
    ${SRW_LATT}/random_walker.cpp
)
target_compile_options(random_walker PRIVATE -O3)
install(TARGETS random_walker DESTINATION lrgsglib/statsys/SignedRW/ccore/Lattices)

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
