#ifdef _MSC_VER
#pragma warning(disable: 26495)
#pragma warning(disable: 6011)
#pragma warning(disable: 4005)
#endif

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <windows.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#ifndef BI_BOARDTYPE
#define BI_BOARDTYPE 0
#endif

#include "cbw.h"

namespace py = pybind11;

// --- Variables Globales ---
static HGLOBAL g_mem_handle = nullptr;
static int g_range = 0;
static long g_rate = 0;
static int g_num_chans_scanned = 0;
static long g_total_buffer_points = 0;
static long long g_total_scans_read = 0;
static long g_last_index = 0; // <<< MODIFICADO: Hacer last_index global

void start_stream(int board_num, int low_chan, int high_chan, int num_samples_per_chan, int rate, int range) {
    if (g_mem_handle) {
        cbStopBackground(board_num, AIFUNCTION);
        cbWinBufFree(g_mem_handle);
        g_mem_handle = nullptr;
    }

    g_num_chans_scanned = high_chan - low_chan + 1;
    if (g_num_chans_scanned <= 0) {
        throw std::runtime_error("High channel must be greater than or equal to low channel.");
    }

    g_total_buffer_points = num_samples_per_chan * g_num_chans_scanned;

    // NO reiniciamos g_total_scans_read aquí para permitir la reanudación
    g_rate = rate;
    g_range = range;
    g_last_index = 0; // <<< MODIFICADO: Reiniciar el índice del búfer en cada nuevo stream

    g_mem_handle = cbWinBufAlloc(g_total_buffer_points);
    if (!g_mem_handle) {
        throw std::runtime_error("Failed to allocate DAQ memory buffer.");
    }

    long actual_rate = rate;
    int err = cbAInScan(board_num, low_chan, high_chan, g_total_buffer_points, &actual_rate, range, g_mem_handle, BACKGROUND | CONTINUOUS);

    g_rate = actual_rate;

    if (err != 0) {
        cbWinBufFree(g_mem_handle);
        g_mem_handle = nullptr;
        char msg[ERRSTRLEN];
        cbGetErrMsg(err, msg);
        throw std::runtime_error("cbAInScan failed: " + std::string(msg));
    }
}

py::array_t<double> read_block(int board_num) {
    if (!g_mem_handle) {
        return py::array_t<double>(std::vector<py::ssize_t>{0, (py::ssize_t)(1 + g_num_chans_scanned)});
    }

    short status = 0;
    long current_count = 0;
    long current_index = 0;

    int err = cbGetIOStatus(board_num, &status, &current_count, &current_index, AIFUNCTION);
    if (err != 0) {
        char msg[ERRSTRLEN];
        cbGetErrMsg(err, msg);
        throw std::runtime_error("cbGetIOStatus failed: " + std::string(msg));
    }

    WORD* buffer = (WORD*)g_mem_handle;
    long points_available = current_index - g_last_index;
    if (points_available < 0) {
        points_available += g_total_buffer_points;
    }

    if (g_num_chans_scanned == 0 || points_available < g_num_chans_scanned) {
        return py::array_t<double>(std::vector<py::ssize_t>{0, (py::ssize_t)(1 + g_num_chans_scanned)});
    }

    long num_scans_available = points_available / g_num_chans_scanned;
    if (num_scans_available == 0) {
        return py::array_t<double>(std::vector<py::ssize_t>{0, (py::ssize_t)(1 + g_num_chans_scanned)});
    }

    long points_to_process = num_scans_available * g_num_chans_scanned;

    py::array_t<double> result(std::vector<py::ssize_t>{ (py::ssize_t)num_scans_available, (py::ssize_t)(1 + g_num_chans_scanned) });
    auto r_acc = result.mutable_unchecked<2>();

    double dt = 1.0 / g_rate;

    // <<< MODIFICADO: Lógica para forzar el inicio en t=0 >>>
    // Si es el primer bloque de datos (g_total_scans_read es 0), guardamos el tiempo del primer punto para usarlo como offset
    double time_offset = 0;
    if (g_total_scans_read == 0) {
        time_offset = (g_total_scans_read)*dt;
    }

    for (long i = 0; i < num_scans_available; ++i) {
        r_acc(i, 0) = (g_total_scans_read + i) * dt - time_offset; // Aplicamos el offset
        for (int j = 0; j < g_num_chans_scanned; ++j) {
            long buffer_idx = (g_last_index + i * g_num_chans_scanned + j) % g_total_buffer_points;
            float voltage;
            cbToEngUnits(board_num, g_range, buffer[buffer_idx], &voltage);
            r_acc(i, j + 1) = (double)voltage;
        }
    }

    g_last_index = (g_last_index + points_to_process) % g_total_buffer_points;
    g_total_scans_read += num_scans_available;

    return result;
}

// <<< MODIFICADO: stop_stream ahora es una "pausa", no reinicia el contador de tiempo >>>
void stop_stream(int board_num) {
    if (g_mem_handle) {
        cbStopBackground(board_num, AIFUNCTION);
        cbWinBufFree(g_mem_handle);
        g_mem_handle = nullptr;
        // YA NO REINICIAMOS g_total_scans_read aquí.
    }
}

// <<< MODIFICADO: Nueva función para el reinicio completo del tiempo >>>
void reset_time_counter() {
    g_total_scans_read = 0;
}

// --- Módulo Pybind11 ---
PYBIND11_MODULE(USB_1608FS_PLUS, m) {
    m.doc() = "Python bindings for OM-USB-1608FS-PLUS";

    m.def("start_stream", &start_stream, "Starts the data acquisition stream",
        py::arg("board_num"), py::arg("low_chan"), py::arg("high_chan"),
        py::arg("num_samples_per_chan"), py::arg("rate"), py::arg("range"));

    m.def("read_block", &read_block, "Reads a block of data from the stream", py::arg("board_num"));
    m.def("stop_stream", &stop_stream, "Stops (pauses) the data acquisition stream", py::arg("board_num"));

    // <<< MODIFICADO: Exponer la nueva función de reseteo a Python >>>
    m.def("reset_time_counter", &reset_time_counter, "Resets the absolute time counter to zero");

    m.attr("BIP10VOLTS") = BIP10VOLTS;
    m.attr("BIP5VOLTS") = BIP5VOLTS;
    m.attr("BIP2VOLTS") = BIP2VOLTS;
    m.attr("BIP1VOLTS") = BIP1VOLTS;
}