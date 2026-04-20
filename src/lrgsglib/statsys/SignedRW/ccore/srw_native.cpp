/**
 * @file srw_native.cpp
 * @brief Pybind11 wrapper for the SignedRW absorb-walker C kernel.
 *
 * Provides a Python-callable function that wraps
 * ``run_absorb_walker`` (see LRGSG_srw.h) for direct invocation with
 * numpy arrays --- no subprocess, no file I/O.
 *
 * Python side:
 *   runlang = "pb_absorb"  ->  absorb_walker_sampling()
 */
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>

extern "C" {
#include "LRGSG_srw.h"
#include "sfmtrng.h"
#include "SFMT/SFMT.h"
}

namespace py = pybind11;


/* ==================================================================
 * Helper: seed the SFMT RNG directly from the user-supplied seed.
 *
 * Note: sfmtrng.c's `__set_seed_SFMT()` hard-codes seed components from
 * time/pid and ignores the external `seed_rand` pointer (the Ising
 * wrapper's `seed_rand = seed_arr` line is therefore a no-op).  We
 * bypass that and call `sfmt_init_by_array` directly so the caller's
 * `seed_val` actually controls the RNG stream.
 * ================================================================== */
static void seed_rng(uint64_t seed_val) {
    uint32_t seed_arr[4];
    seed_arr[0] = static_cast<uint32_t>(seed_val & 0xFFFFFFFF);
    seed_arr[1] = static_cast<uint32_t>((seed_val >> 32) & 0xFFFFFFFF);
    seed_arr[2] = static_cast<uint32_t>(seed_val * 0x9E3779B97F4A7C15ULL);
    seed_arr[3] = static_cast<uint32_t>(seed_val ^ 0xBE11AC1A0ULL);
    sfmt_init_by_array(&sfmt, seed_arr, 4);
}


/* ==================================================================
 * Absorb-walker sampling
 * ================================================================== */
static py::tuple absorb_walker_sampling(
    size_t N,
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<int64_t>& neigh_ptr,
    const py::array_t<uint8_t>& kill_mask,
    const py::array_t<uint8_t>& reflect_mask,
    const py::array_t<int64_t>& start_positions,
    int64_t stop_thresh,
    uint64_t seed_val
) {
    /* ---- validate CSR shapes ---- */
    auto np_buf = neigh_ptr.request();
    if (np_buf.ndim != 1 || static_cast<size_t>(np_buf.size) != N + 1)
        throw std::runtime_error("neigh_ptr must have shape (N+1,).");
    auto ni_buf = neigh_indices.request();
    auto km_buf = kill_mask.request();
    auto rm_buf = reflect_mask.request();
    if (ni_buf.size != km_buf.size || ni_buf.size != rm_buf.size)
        throw std::runtime_error(
            "neigh_indices / kill_mask / reflect_mask must have identical shape."
        );

    auto sp_buf = start_positions.request();
    size_t n_walkers = static_cast<size_t>(sp_buf.size);
    if (n_walkers == 0)
        throw std::runtime_error("start_positions must be non-empty.");

    /* ---- build size_t neigh_ptr (C kernel expects size_t) ---- */
    std::vector<size_t> ptr_buf(N + 1);
    const int64_t *np_ptr = static_cast<const int64_t*>(np_buf.ptr);
    for (size_t i = 0; i <= N; ++i) {
        ptr_buf[i] = static_cast<size_t>(np_ptr[i]);
    }

    const int64_t *ni_ptr = static_cast<const int64_t*>(ni_buf.ptr);
    const uint8_t *km_ptr = static_cast<const uint8_t*>(km_buf.ptr);
    const uint8_t *rm_ptr = static_cast<const uint8_t*>(rm_buf.ptr);
    const int64_t *sp_ptr = static_cast<const int64_t*>(sp_buf.ptr);

    /* ---- allocate output buffers ---- */
    std::vector<int64_t> visits_agg(N, 0);
    std::vector<int64_t> unique_visits(n_walkers, 0);
    std::vector<int64_t> stop_step(n_walkers, -1);
    std::vector<int8_t>  stop_reason(n_walkers, -1);
    std::vector<int64_t> final_position(n_walkers, 0);

    {
        py::gil_scoped_release release;
        seed_rng(seed_val);
        run_absorb_walker(
            N,
            n_walkers,
            ptr_buf.data(),
            ni_ptr,
            km_ptr,
            rm_ptr,
            sp_ptr,
            stop_thresh,
            visits_agg.data(),
            unique_visits.data(),
            stop_step.data(),
            stop_reason.data(),
            final_position.data()
        );
    }

    /* ---- wrap outputs as numpy arrays ---- */
    py::array_t<int64_t> out_visits(N);
    std::memcpy(out_visits.mutable_data(), visits_agg.data(), N * sizeof(int64_t));

    py::array_t<int64_t> out_unique(n_walkers);
    std::memcpy(out_unique.mutable_data(), unique_visits.data(),
                n_walkers * sizeof(int64_t));

    py::array_t<int64_t> out_stop_step(n_walkers);
    std::memcpy(out_stop_step.mutable_data(), stop_step.data(),
                n_walkers * sizeof(int64_t));

    py::array_t<int8_t> out_stop_reason(n_walkers);
    std::memcpy(out_stop_reason.mutable_data(), stop_reason.data(),
                n_walkers * sizeof(int8_t));

    py::array_t<int64_t> out_final_pos(n_walkers);
    std::memcpy(out_final_pos.mutable_data(), final_position.data(),
                n_walkers * sizeof(int64_t));

    return py::make_tuple(
        out_visits, out_unique, out_stop_step, out_stop_reason, out_final_pos);
}


PYBIND11_MODULE(_srw_native, m) {
    m.doc() = R"pbdoc(
        Native SignedRW absorb-walker kernel via pybind11.

        Wraps the C kernel `run_absorb_walker` in
        statsys/SignedRW/ccore/LRGSG_srw.c for in-process invocation
        from Python without subprocess or file I/O.  Kill / sticky
        rules are Python-only in Phase 2; they will be added in a
        follow-up C kernel.
    )pbdoc";

    m.def("absorb_walker_sampling", &absorb_walker_sampling,
        R"pbdoc(
Run an ensemble of absorbing signed random walkers on a CSR graph.

Parameters
----------
N : int
    Number of nodes.
neigh_indices : ndarray[int64]
    Concatenated neighbour indices (length = total degree).
neigh_ptr : ndarray[int64]
    CSR row pointers of length N+1.
kill_mask : ndarray[uint8]
    1 where moving across edge e kills the walker; 0 otherwise.
reflect_mask : ndarray[uint8]
    1 where moving across edge e reflects the walker (stays put); 0 otherwise.
start_positions : ndarray[int64]
    Initial node index for each walker (length n_walkers).
stop_thresh : int
    Walker terminates once its unique-visit count reaches this value.
seed : int
    SFMT RNG seed.

Returns
-------
tuple[ndarray, ndarray, ndarray, ndarray, ndarray]
    (visits_agg, unique_visits, stop_step, stop_reason, final_position)
        )pbdoc",
        py::arg("N"),
        py::arg("neigh_indices"),
        py::arg("neigh_ptr"),
        py::arg("kill_mask"),
        py::arg("reflect_mask"),
        py::arg("start_positions"),
        py::arg("stop_thresh"),
        py::arg("seed"));
}
