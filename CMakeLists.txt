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

# SIMD detection (platform-aware)
include(CheckCCompilerFlag)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm|aarch64|ARM64")
    # ARM: check for NEON (default on AArch64)
    check_c_compiler_flag("-march=armv8-a" HAS_NEON)
    if(HAS_NEON)
        add_compile_definitions(HAVE_NEON)
    endif()
else()
    # x86/x86_64: check for SSE2
    check_c_compiler_flag("-msse2" HAS_SSE2)
    if(HAS_SSE2)
        add_compile_options(-msse2)
        add_compile_definitions(HAVE_SSE2)
    endif()
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
# Voter shares the gillespie CTMC kernel + cluster recording from _ccore; both
# the VoterSimulator exe and the _voter_native pybind module (which list SRC_VM)
# must link them, else: undefined reference to voter_ctmc_run / clusters_record.
set(SRC_VM
    ${VM_DIR}/LRGSG_vm.c
    ${CCORE}/LRGSG_ctmc.c       # voter_ctmc_run (gillespie: composition-rejection + Fenwick)
    ${CCORE}/LRGSG_clusters.c   # clusters_record (fast cluster-size recompute)
    ${CCORE}/LRGSG_uf.c         # union-find used by cluster recording
)
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

# Unified Ising simulators (canonical names used by Python dispatch)
add_ising_simulator(IsingMetropolis)
add_ising_sa_simulator(IsingSimulatedAnnealing)
add_ising_pt_simulator(IsingParallelTempering)

# NB: the legacy IsingSimulator0..5 / 1b / 3b / 4b targets were removed -- their
# .c sources were dropped in the folder-per-class reorg, so a fresh checkout had
# no files for them and CMake configure failed ("No SOURCES given to target").
# The canonical Metropolis / SA / PT binaries above are the ones Python dispatches to.

# ===========================================================================
# Voter Model simulators
# ===========================================================================
add_voter_simulator(VoterSimulator)

# ===========================================================================
# Contact Process simulators
# ===========================================================================
add_contact_simulator(ContactProcessSIR)
add_contact_simulator(ContactProcessEI)

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

# Ising native kernels - pybind11 binding (Metropolis, SA, PT via numpy).
# LRGSG_cluster.c + LRGSG_uf.c mirror ccore/native/setup.py (the working
# Makefile build) -- without them the module fails to link.
pybind11_add_module(_ising_native
    ${ISING_CCORE}/ising_native.cpp
    ${SRC_RBIM}
    ${SRC_SA}
    ${SRC_PT}
    ${ISING_CCORE}/LRGSG_cluster.c
    ${CCORE}/LRGSG_uf.c
    ${COMMON_SOURCES}
)
target_include_directories(_ising_native PRIVATE
    ${COMMON_INCLUDES}
    ${Python_NumPy_INCLUDE_DIRS}
)
target_compile_options(_ising_native PRIVATE -O3)
install(TARGETS _ising_native DESTINATION lrgsglib/statsys/IsingDynamics/ccore)

# Voter native kernels - pybind11 binding (in-process voter dynamics via numpy)
pybind11_add_module(_voter_native
    ${VM_DIR}/voter_native.cpp
    ${SRC_VM}
    ${COMMON_SOURCES}
)
target_include_directories(_voter_native PRIVATE
    ${COMMON_INCLUDES}
    ${VM_DIR}
    ${Python_NumPy_INCLUDE_DIRS}
)
target_compile_options(_voter_native PRIVATE -O3)
install(TARGETS _voter_native DESTINATION lrgsglib/statsys/VoterModel/ccore)

# Random walker on lattices - pybind11 binding
pybind11_add_module(random_walker
    ${SRW_LATT}/random_walker.cpp
)
target_compile_options(random_walker PRIVATE -O3)
install(TARGETS random_walker DESTINATION lrgsglib/statsys/SignedRW/ccore/Lattices)

# ===========================================================================
# In-process native kernels (one pybind module per model)
# ===========================================================================
# Folded here from the per-model ccore/native/{setup.py,Makefile} builds
# (statsys unification, decision D-B7): ONE build system for every compiled
# artifact. Source lists mirror each setup.py exactly. The per-model
# Makefiles keep working for one-off local rebuilds.
function(add_native_module NAME DEST)
    cmake_parse_arguments(NAT "" "" "SOURCES;INCLUDES" ${ARGN})
    pybind11_add_module(${NAME} ${NAT_SOURCES})
    target_include_directories(${NAME} PRIVATE
        ${COMMON_INCLUDES}
        ${NAT_INCLUDES}
        ${Python_NumPy_INCLUDE_DIRS}
    )
    target_compile_options(${NAME} PRIVATE -O3)
    install(TARGETS ${NAME} DESTINATION ${DEST})
endfunction()

add_native_module(_potts_native lrgsglib/statsys/PottsModel/ccore
    SOURCES ${POTTS_DIR}/potts_native.cpp ${SRC_POTTS} ${SRC_VECDYNSYS}
            ${COMMON_SOURCES}
    INCLUDES ${POTTS_DIR}
)

add_native_module(_xy_native lrgsglib/statsys/XYModel/ccore
    SOURCES ${XY_DIR}/xy_native.cpp ${SRC_XY} ${SRC_CONTDYNSYS}
            ${SRC_VECDYNSYS} ${COMMON_SOURCES}
    INCLUDES ${XY_DIR}
)

add_native_module(_heisenberg_native lrgsglib/statsys/HeisenbergModel/ccore
    SOURCES ${HBERG_DIR}/heisenberg_native.cpp ${SRC_HBERG}
            ${SRC_VECDYNSYS} ${COMMON_SOURCES}
    INCLUDES ${HBERG_DIR}
)

add_native_module(_multispec_native lrgsglib/statsys/MultiSpeciesModel/ccore
    SOURCES ${MSPEC_DIR}/multispec_native.cpp ${SRC_MSPEC}
            ${SRC_VECDYNSYS} ${COMMON_SOURCES}
    INCLUDES ${MSPEC_DIR}
)

add_native_module(_kuramoto_native lrgsglib/statsys/KuramotoModel/ccore
    SOURCES ${KUR_DIR}/kuramoto_native.cpp ${SRC_KURAMOTO}
            ${SRC_CONTDYNSYS} ${COMMON_SOURCES}
    INCLUDES ${KUR_DIR}
)

add_native_module(_rd_native lrgsglib/statsys/ReactionDiffusionModel/ccore
    SOURCES ${RD_DIR}/rd_native.cpp ${SRC_RD} ${SRC_CONTDYNSYS}
            ${COMMON_SOURCES}
    INCLUDES ${RD_DIR}
)

add_native_module(_code_native lrgsglib/statsys/CoupledODEModel/ccore
    SOURCES ${CODE_DIR}/code_native.cpp ${SRC_CODE} ${SRC_CONTDYNSYS}
            ${COMMON_SOURCES}
    INCLUDES ${CODE_DIR}
)

# SignedRW walker kernel: no BinDynSys dependency (its setup.py lists only
# the walker + RNG sources).
add_native_module(_srw_native lrgsglib/statsys/SignedRW/ccore
    SOURCES ${SRW_CCORE}/srw_native.cpp ${SRW_CCORE}/LRGSG_srw.c
            ${CCORE}/LRGSG_utils.c ${CCORE}/sfmtrng.c ${SFMT_DIR}/SFMT.c
    INCLUDES ${SRW_CCORE}
)

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
