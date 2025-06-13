#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <limits>
#include <complex>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#ifdef _WIN32
#include <windows.h>
#include <Commdlg.h>
#pragma comment(lib, "Comdlg32.lib")
#endif

const float PI_DEVICE_VAL = 3.1415926535f;
const int NUM_SUBWOOFERS_FIXED = 4; // Max number of subs fixed at 4
const int OPTIMAL_THREADS_PER_BLOCK = 512;

__constant__ int c_num_freq_points_k;
__constant__ float c_max_spread_s_param_k;
__constant__ float c_PI_DEVICE_K;
__constant__ int c_num_active_subwoofers_k;

inline void gpuAssert(cudaError_t code, const char* file, int line, bool abort = true) {
    if (code != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) {
#ifdef _WIN32
            std::cerr << "Press Enter to exit..." << std::endl;
            std::cin.clear(); std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); std::cin.get();
#endif
            exit(code);
        }
    }
}
#define CUDA_CHECK(ans) { gpuAssert((ans), __FILE__, __LINE__); }

struct DeviceOptResultPayload4Subs {
    float metric_to_maximize;
    float delay_s[NUM_SUBWOOFERS_FIXED];
    float pol[NUM_SUBWOOFERS_FIXED];
};

struct OptResultPayload4Subs {
    double metric_to_maximize;
    float delay_s[NUM_SUBWOOFERS_FIXED];
    float pol[NUM_SUBWOOFERS_FIXED];
};

struct ChannelDataSourceSub {
    std::vector<float> frequencies_hz;
    std::vector<float> magnitudes_linear;
    std::vector<float> phases_radian;
    size_t num_points = 0;
    std::string filepath_source;
    bool load_and_process_file(const std::string& filepath, float min_freq = 15.0f, float max_freq = 250.0f) {
        this->filepath_source = filepath;
        std::ifstream file(filepath);
        if (!file.is_open()) { std::cerr << "Error opening " << filepath << std::endl; return false; }
        frequencies_hz.clear(); magnitudes_linear.clear(); phases_radian.clear(); num_points = 0;
        std::string line;
        char bom_check[3]; file.read(bom_check, 3);
        if (!(bom_check[0] == (char)0xEF && bom_check[1] == (char)0xBB && bom_check[2] == (char)0xBF)) { file.seekg(0); }
        while (std::getline(file, line)) {
            line.erase(0, line.find_first_not_of(" \t\n\r\f\v"));
            line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
            if (line.empty() || line[0] == '*' || line[0] == '#' || line.rfind("//", 0) == 0) continue;
            std::string p_line = line; std::replace(p_line.begin(), p_line.end(), '\t', ',');
            std::stringstream ss(p_line); std::string seg; std::vector<float> parts;
            try {
                while (std::getline(ss, seg, ',')) {
                    seg.erase(0, seg.find_first_not_of(" \t\n\r\f\v"));
                    seg.erase(seg.find_last_not_of(" \t\n\r\f\v") + 1);
                    if (!seg.empty()) parts.push_back(std::stof(seg));
                }
            }
            catch (const std::exception&) { continue; }
            if (parts.size() == 3) {
                float fq = parts[0];
                if (fq >= min_freq && fq <= max_freq && std::isfinite(parts[1]) && std::isfinite(parts[2])) {
                    frequencies_hz.push_back(fq);
                    magnitudes_linear.push_back(powf(10.0f, parts[1] / 20.0f));
                    phases_radian.push_back(parts[2] * (PI_DEVICE_VAL / 180.0f));
                }
            }
        }
        file.close();
        if (frequencies_hz.empty()) { std::cerr << "No valid data in " << min_freq << "-" << max_freq << "Hz from " << filepath << std::endl; return false; }
        num_points = frequencies_hz.size();
        return true;
    }
};

__global__ void optimize_4subs_only_kernel(
    const float* d_frequencies_f,
    const float* d_sub_mags_flat,
    const float* d_sub_phases_flat,
    const float* delay_steps_s_f_d,
    int num_discrete_delay_steps,
    float* d_block_best_results_output
) {
    extern __shared__ float s_data[];
    float* s_frequencies = s_data;
    DeviceOptResultPayload4Subs* s_payload_cache = (DeviceOptResultPayload4Subs*)&s_frequencies[c_num_freq_points_k];
    int tid = threadIdx.x;
    int block_size = blockDim.x;
    int block_id = blockIdx.x;
    int global_tid = block_id * block_size + tid;
    int total_grid_threads = gridDim.x * block_size;
    for (int i = tid; i < c_num_freq_points_k; i += block_size) {
        s_frequencies[i] = d_frequencies_f[i];
    }
    __syncthreads();
    float thread_best_metric = -1.0f;
    float thread_best_d_s[NUM_SUBWOOFERS_FIXED] = { 0.0f };
    float thread_best_p[NUM_SUBWOOFERS_FIXED] = { 0.0f };
    unsigned long long n_delays_ull = num_discrete_delay_steps;
    unsigned long long n_delays_power_n_subs = 1;

    if (n_delays_ull > 1) {
        for (int s_iter = 0; s_iter < c_num_active_subwoofers_k; ++s_iter) {
            if (n_delays_power_n_subs > ULLONG_MAX / n_delays_ull) {
                n_delays_power_n_subs = ULLONG_MAX; break;
            }
            n_delays_power_n_subs *= n_delays_ull;
        }
    }
    unsigned long long num_pol_comb = (1ULL << c_num_active_subwoofers_k);
    unsigned long long total_combinations;
    total_combinations = (n_delays_power_n_subs > 0 && num_pol_comb > ULLONG_MAX / n_delays_power_n_subs) ? ULLONG_MAX : n_delays_power_n_subs * num_pol_comb;
    if (total_combinations == 0) total_combinations = 1;
    float current_d_s_f_arr[NUM_SUBWOOFERS_FIXED];
    float pol_val_s_arr[NUM_SUBWOOFERS_FIXED];
    for (unsigned long long comb_idx = global_tid; comb_idx < total_combinations; comb_idx += total_grid_threads) {
        unsigned long long temp_idx = comb_idx;
        for (int s_idx = 0; s_idx < NUM_SUBWOOFERS_FIXED; ++s_idx) {
            pol_val_s_arr[s_idx] = 0.0f;
            current_d_s_f_arr[s_idx] = (num_discrete_delay_steps > 0) ? delay_steps_s_f_d[0] : 0.0f;
        }
        for (int s_idx = 0; s_idx < c_num_active_subwoofers_k; ++s_idx) {
            pol_val_s_arr[s_idx] = (float)(temp_idx & 1);
            temp_idx >>= 1;
        }
        if (n_delays_ull > 1) {
            for (int s_idx = 0; s_idx < c_num_active_subwoofers_k; ++s_idx) {
                current_d_s_f_arr[s_idx] = delay_steps_s_f_d[temp_idx % n_delays_ull];
                temp_idx /= n_delays_ull;
            }
        }
        bool spread_constraint_violated = false;
        for (int i = 0; i < c_num_active_subwoofers_k; ++i) {
            for (int j = i + 1; j < c_num_active_subwoofers_k; ++j) {
                if (fabsf(current_d_s_f_arr[i] - current_d_s_f_arr[j]) > (c_max_spread_s_param_k + 1e-7f)) {
                    spread_constraint_violated = true;
                }
            }
        }
        if (spread_constraint_violated) continue;
        float current_sum_linear_mags = 0.0f;
        for (int f_idx = 0; f_idx < c_num_freq_points_k; ++f_idx) {
            float total_r_freq = 0.0f;
            float total_i_freq = 0.0f;
            float omega_f = 2.0f * c_PI_DEVICE_K * s_frequencies[f_idx];
            size_t sub_base_idx = (size_t)f_idx * NUM_SUBWOOFERS_FIXED;
#pragma unroll
            for (int s_idx = 0; s_idx < NUM_SUBWOOFERS_FIXED; ++s_idx) {
                size_t sub_data_idx = sub_base_idx + s_idx;
                float mag = d_sub_mags_flat[sub_data_idx];
                float phase_orig = d_sub_phases_flat[sub_data_idx];
                float phase_shifted = phase_orig + omega_f * current_d_s_f_arr[s_idx] + pol_val_s_arr[s_idx] * c_PI_DEVICE_K;
                float s, c;
                sincosf(phase_shifted, &s, &c);
                total_r_freq += mag * c;
                total_i_freq += mag * s;
            }
            current_sum_linear_mags += sqrtf(total_r_freq * total_r_freq + total_i_freq * total_i_freq);
        }
        if (current_sum_linear_mags > thread_best_metric) {
            thread_best_metric = current_sum_linear_mags;
#pragma unroll
            for (int s_idx = 0; s_idx < NUM_SUBWOOFERS_FIXED; ++s_idx) {
                thread_best_d_s[s_idx] = current_d_s_f_arr[s_idx];
                thread_best_p[s_idx] = pol_val_s_arr[s_idx];
            }
        }
    }
    s_payload_cache[tid].metric_to_maximize = thread_best_metric;
#pragma unroll
    for (int s_idx = 0; s_idx < NUM_SUBWOOFERS_FIXED; ++s_idx) {
        s_payload_cache[tid].delay_s[s_idx] = thread_best_d_s[s_idx];
        s_payload_cache[tid].pol[s_idx] = thread_best_p[s_idx];
    }
    __syncthreads();
    for (int s_val = block_size / 2; s_val > 0; s_val >>= 1) {
        if (tid < s_val) {
            if (s_payload_cache[tid + s_val].metric_to_maximize > s_payload_cache[tid].metric_to_maximize) {
                s_payload_cache[tid] = s_payload_cache[tid + s_val];
            }
        }
        __syncthreads();
    }
    if (tid == 0 && s_payload_cache[0].metric_to_maximize >= 0.0f) {
        int out_idx = block_id * (1 + NUM_SUBWOOFERS_FIXED * 2);
        d_block_best_results_output[out_idx + 0] = s_payload_cache[0].metric_to_maximize;
#pragma unroll
        for (int s_idx = 0; s_idx < NUM_SUBWOOFERS_FIXED; ++s_idx) {
            d_block_best_results_output[out_idx + 1 + s_idx * 2 + 0] = s_payload_cache[0].delay_s[s_idx];
            d_block_best_results_output[out_idx + 1 + s_idx * 2 + 1] = s_payload_cache[0].pol[s_idx];
        }
    }
}

#ifdef _WIN32
std::string get_open_filepath_subs(const char* dialog_title_param) {
    OPENFILENAMEA ofn; ZeroMemory(&ofn, sizeof(ofn)); char szFile[MAX_PATH] = { 0 };
    ofn.lStructSize = sizeof(ofn); ofn.hwndOwner = GetConsoleWindow(); ofn.lpstrFile = szFile; ofn.nMaxFile = sizeof(szFile);
    ofn.lpstrFilter = "Text Files (*.txt)\0*.txt\0All Files (*.*)\0*.*\0"; ofn.nFilterIndex = 1;
    ofn.lpstrTitle = dialog_title_param; ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST | OFN_NOCHANGEDIR;
    if (GetOpenFileNameA(&ofn) == TRUE) return std::string(ofn.lpstrFile);
    return "";
}
#endif

std::string format_with_commas_subs(unsigned long long n) {
    if (n == ULLONG_MAX) return "overflow";
    std::string s = std::to_string(n); if (s.empty()) return "0"; int len = s.length();
    int num_commas = (len - 1) / 3; if (num_commas == 0) return s; std::string formatted_s;
    formatted_s.resize(len + num_commas); int read_idx = len - 1; int write_idx = formatted_s.length() - 1;
    int count = 0; while (read_idx >= 0) {
        formatted_s[write_idx--] = s[read_idx--]; count++;
        if (count == 3 && read_idx >= 0) { formatted_s[write_idx--] = ','; count = 0; }
    } return formatted_s;
}

float get_user_float(const std::string& prompt, float default_val) {
    std::cout << prompt << " (default: " << default_val << "): ";
    std::string line;
    std::getline(std::cin, line);
    if (line.empty()) {
        return default_val;
    }
    try {
        return std::stof(line);
    }
    catch (const std::exception&) {
        std::cout << "Invalid input. Using default value." << std::endl;
        return default_val;
    }
}

std::string get_filepath_from_args(int argc, char** argv, int index) {
    if (index + 1 >= argc) {
        return "";
    }
    return argv[index + 1];
}

int main(int argc, char** argv) {
#ifndef _WIN32
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " sub1.txt [sub2.txt ... sub4.txt]\n";
        return 1;
    }
#endif
    std::cout << "--- Multiple Subwoofer (max. 4) Delay & Inversion Optimizer (GPU) ---" << std::endl;
    std::cout << "\n--- Configuration ---" << std::endl;
    std::cout << "Enter desired parameters or press Enter to use defaults." << std::endl;
    float min_freq = get_user_float("Minimum frequency to analyze (Hz)", 15.0f);
    float max_freq = get_user_float("Maximum frequency to analyze (Hz)", 250.0f);
    std::vector<ChannelDataSourceSub> sub_data_sources(NUM_SUBWOOFERS_FIXED);
    int num_active_subwoofers = 0;
    for (int i = 0; i < NUM_SUBWOOFERS_FIXED; ++i) {
#ifdef _WIN32
        std::string title = "Select Measurement File for Subwoofer " + std::to_string(i + 1) + " (or Cancel to finish)";
        std::cout << "\n" << title << std::endl;
        std::string sub_p = get_open_filepath_subs(title.c_str());
#else
        std::string sub_p = get_filepath_from_args(argc, argv, i);
#endif
        if (sub_p.empty()) {
            if (i == 0) {
                std::cerr << "Must load at least one subwoofer file. Aborting." << std::endl;
                std::cin.get(); return 1;
            }
#ifdef _WIN32
            std::cout << "No more files selected. Proceeding with " << i << " subwoofer(s)." << std::endl;
#endif
            break;
        }
        if (!sub_data_sources[i].load_and_process_file(sub_p, min_freq, max_freq)) {
            std::cerr << "Failed to load or process file. Aborting." << std::endl;
            std::cin.get(); return 1;
        }
        std::cout << "Loaded SW" << (i + 1) << ": " << sub_data_sources[i].filepath_source
            << " (" << sub_data_sources[i].num_points << "pts)" << std::endl;
        num_active_subwoofers++;
    }
    int num_freq_points_host = static_cast<int>(sub_data_sources[0].num_points);
    for (int i = 1; i < num_active_subwoofers; ++i) {
        if (sub_data_sources[i].num_points != num_freq_points_host) {
            std::cerr << "Subwoofer " << (i + 1) << " has different number of points. Aborting." << std::endl;
            std::cin.get(); return 1;
        }
        for (int k = 0; k < num_freq_points_host; ++k) {
            if (std::fabs(sub_data_sources[i].frequencies_hz[k] - sub_data_sources[0].frequencies_hz[k]) > 1e-3f) {
                std::cerr << "Subwoofer " << (i + 1) << " has different frequency axis. Aborting." << std::endl;
                std::cin.get(); return 1;
            }
        }
    }
    std::cout << "\nAll subwoofers use a common frequency axis." << std::endl;
    std::cout << "  Frequency Range: " << std::fixed << std::setprecision(1) << sub_data_sources[0].frequencies_hz.front()
        << " Hz to " << sub_data_sources[0].frequencies_hz.back() << " Hz" << std::endl;
    std::cout << "  Frequency Points: " << num_freq_points_host << std::endl;
    float min_delay_ms = get_user_float("\nMinimum delay (ms)", -20.0f);
    float max_delay_ms = get_user_float("Maximum delay (ms)", 20.0f);
    float step_delay_ms = get_user_float("Delay step (ms)", 0.1f);
    float max_spread_ms = get_user_float("Maximum delay spread between any two subs (ms)", 20.0f);
    const float min_delay_s = min_delay_ms / 1000.0f;
    const float max_delay_s = max_delay_ms / 1000.0f;
    const float step_delay_s = step_delay_ms / 1000.0f;
    const float max_spread_s = max_spread_ms / 1000.0f;
    std::vector<float> h_delay_steps_s_f;
    for (float d_s = min_delay_s; d_s <= max_delay_s + (step_delay_s / 2.0f); d_s += step_delay_s) {
        h_delay_steps_s_f.push_back(d_s);
    }
    int num_discrete_delay_steps_host = static_cast<int>(h_delay_steps_s_f.size());
    if (num_discrete_delay_steps_host == 0) {
        std::cerr << "Error: 0 delay steps generated. Check delay parameters." << std::endl;
        std::cin.get(); return 1;
    }
    std::cout << std::defaultfloat;
    std::cout << "\nUsing " << num_discrete_delay_steps_host << " delay steps from " << min_delay_ms << "ms to " << max_delay_ms << "ms with " << step_delay_ms << "ms step." << std::endl;
    std::vector<float> h_sub_mags_flat((size_t)NUM_SUBWOOFERS_FIXED * num_freq_points_host);
    std::vector<float> h_sub_phases_flat((size_t)NUM_SUBWOOFERS_FIXED * num_freq_points_host);
    for (int f_idx = 0; f_idx < num_freq_points_host; ++f_idx) {
        for (int s_idx = 0; s_idx < NUM_SUBWOOFERS_FIXED; ++s_idx) {
            size_t flat_idx = (size_t)f_idx * NUM_SUBWOOFERS_FIXED + s_idx;
            if (s_idx < num_active_subwoofers) {
                h_sub_mags_flat[flat_idx] = sub_data_sources[s_idx].magnitudes_linear[f_idx];
                h_sub_phases_flat[flat_idx] = sub_data_sources[s_idx].phases_radian[f_idx];
            }
            else {
                h_sub_mags_flat[flat_idx] = 0.0f;
                h_sub_phases_flat[flat_idx] = 0.0f;
            }
        }
    }
    cudaDeviceProp dev_prop; int dev_id = 0;
    CUDA_CHECK(cudaGetDeviceProperties(&dev_prop, dev_id));
    CUDA_CHECK(cudaSetDevice(dev_id));
    std::cout << "Using GPU: " << dev_prop.name << std::endl;
    CUDA_CHECK(cudaMemcpyToSymbol(c_PI_DEVICE_K, &PI_DEVICE_VAL, sizeof(float)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_num_freq_points_k, &num_freq_points_host, sizeof(int)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_max_spread_s_param_k, &max_spread_s, sizeof(float)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_num_active_subwoofers_k, &num_active_subwoofers, sizeof(int)));
    float* d_frequencies_gpu = nullptr;
    float* d_sub_mags_flat_gpu = nullptr;
    float* d_sub_phases_flat_gpu = nullptr;
    float* d_delay_steps_gpu = nullptr;
    CUDA_CHECK(cudaMalloc((void**)&d_frequencies_gpu, sub_data_sources[0].frequencies_hz.size() * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d_frequencies_gpu, sub_data_sources[0].frequencies_hz.data(), sub_data_sources[0].frequencies_hz.size() * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMalloc((void**)&d_sub_mags_flat_gpu, h_sub_mags_flat.size() * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d_sub_mags_flat_gpu, h_sub_mags_flat.data(), h_sub_mags_flat.size() * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMalloc((void**)&d_sub_phases_flat_gpu, h_sub_phases_flat.size() * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d_sub_phases_flat_gpu, h_sub_phases_flat.data(), h_sub_phases_flat.size() * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMalloc((void**)&d_delay_steps_gpu, h_delay_steps_s_f.size() * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d_delay_steps_gpu, h_delay_steps_s_f.data(), h_delay_steps_s_f.size() * sizeof(float), cudaMemcpyHostToDevice));
    int threads_per_block = std::min(OPTIMAL_THREADS_PER_BLOCK, dev_prop.maxThreadsPerBlock);
    size_t shared_mem_bytes = (size_t)num_freq_points_host * sizeof(float) + (size_t)threads_per_block * sizeof(DeviceOptResultPayload4Subs);
    if (shared_mem_bytes > dev_prop.sharedMemPerBlock) {
        std::cerr << "Error: Required shared memory (" << shared_mem_bytes << " bytes) exceeds device limit (" << dev_prop.sharedMemPerBlock << " bytes)." << std::endl;
        std::cin.get(); return 1;
    }
    unsigned long long n_delays_ul_host = num_discrete_delay_steps_host;
    unsigned long long n_delays_power_n_subs_host = 1;
    if (n_delays_ul_host > 1) {
        for (int i = 0; i < num_active_subwoofers; ++i) {
            if (n_delays_power_n_subs_host > ULLONG_MAX / n_delays_ul_host) {
                n_delays_power_n_subs_host = ULLONG_MAX; break;
            }
            n_delays_power_n_subs_host *= n_delays_ul_host;
        }
    }
    unsigned long long num_pol_comb_host = (1ULL << num_active_subwoofers);
    unsigned long long total_combinations_host_calc;
    total_combinations_host_calc = (n_delays_power_n_subs_host > 0 && num_pol_comb_host > ULLONG_MAX / n_delays_power_n_subs_host) ? ULLONG_MAX : n_delays_power_n_subs_host * num_pol_comb_host;
    if (total_combinations_host_calc == 0) total_combinations_host_calc = 1;
    std::cout << "Total combinations to check: " << format_with_commas_subs(total_combinations_host_calc) << std::endl;
    int blocks_per_grid = dev_prop.multiProcessorCount * 32;
    if (total_combinations_host_calc != ULLONG_MAX && total_combinations_host_calc < (unsigned long long)blocks_per_grid * threads_per_block) {
        blocks_per_grid = (int)((total_combinations_host_calc + threads_per_block - 1) / threads_per_block);
    }
    blocks_per_grid = std::max(1, blocks_per_grid);
    blocks_per_grid = std::min({ blocks_per_grid, dev_prop.maxGridSize[0] });
    std::cout << "Launching kernel: " << blocks_per_grid << " blocks, " << threads_per_block << " threads/block." << std::endl;
    int results_per_block_elements = 1 + NUM_SUBWOOFERS_FIXED * 2;
    std::vector<float> h_block_results((size_t)blocks_per_grid * results_per_block_elements, -1.0f);
    float* d_block_results_gpu = nullptr;
    CUDA_CHECK(cudaMalloc((void**)&d_block_results_gpu, h_block_results.size() * sizeof(float)));
    CUDA_CHECK(cudaMemset(d_block_results_gpu, 0, h_block_results.size() * sizeof(float)));
    cudaEvent_t start_ev, stop_ev;
    CUDA_CHECK(cudaEventCreate(&start_ev)); CUDA_CHECK(cudaEventCreate(&stop_ev));
    auto host_wall_start_tm = std::chrono::high_resolution_clock::now();
    CUDA_CHECK(cudaEventRecord(start_ev, 0));
    optimize_4subs_only_kernel << <blocks_per_grid, threads_per_block, shared_mem_bytes >> > (
        d_frequencies_gpu,
        d_sub_mags_flat_gpu,
        d_sub_phases_flat_gpu,
        d_delay_steps_gpu,
        num_discrete_delay_steps_host,
        d_block_results_gpu
        );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaEventRecord(stop_ev, 0));
    CUDA_CHECK(cudaEventSynchronize(stop_ev));
    float gpu_ms;
    CUDA_CHECK(cudaEventElapsedTime(&gpu_ms, start_ev, stop_ev));
    auto host_wall_end_tm = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> host_wall_s = host_wall_end_tm - host_wall_start_tm;
    CUDA_CHECK(cudaMemcpy(h_block_results.data(), d_block_results_gpu, h_block_results.size() * sizeof(float), cudaMemcpyDeviceToHost));
    double overall_best_metric = -std::numeric_limits<double>::infinity();
    OptResultPayload4Subs best_config_final;
    for (int i = 0; i < blocks_per_grid; ++i) {
        int res_idx = i * results_per_block_elements;
        float block_metric_f = h_block_results[res_idx];
        if (static_cast<double>(block_metric_f) > overall_best_metric) {
            overall_best_metric = static_cast<double>(block_metric_f);
            best_config_final.metric_to_maximize = overall_best_metric;
            for (int s_idx = 0; s_idx < NUM_SUBWOOFERS_FIXED; ++s_idx) {
                best_config_final.delay_s[s_idx] = h_block_results[res_idx + 1 + s_idx * 2 + 0];
                best_config_final.pol[s_idx] = h_block_results[res_idx + 1 + s_idx * 2 + 1];
            }
        }
    }
    std::cout << "\n--- " << num_active_subwoofers << "-Subwoofer Optimization Complete ---" << std::endl;
    std::cout << "Host Wall Time: " << std::fixed << std::setprecision(3) << host_wall_s.count() << " s" << std::endl;
    std::cout << "GPU Kernel Time: " << std::fixed << std::setprecision(3) << gpu_ms / 1000.0 << " s" << std::endl;
    if (overall_best_metric >= 0.0) {
        std::cout << "\nBest Configuration Found:" << std::endl;
        float min_ref_delay_s = std::numeric_limits<float>::max();
        if (num_active_subwoofers > 0) {
            for (int s_idx = 0; s_idx < num_active_subwoofers; ++s_idx) {
                min_ref_delay_s = std::min(min_ref_delay_s, best_config_final.delay_s[s_idx]);
            }
        }
        else {
            min_ref_delay_s = 0.0f; // Should not happen, but safe
        }
        std::cout << "(Delays are shown relative to the earliest subwoofer)" << std::endl;
        std::cout << std::fixed;
        for (int s_idx = 0; s_idx < num_active_subwoofers; ++s_idx) {
            // Subtract the minimum delay from each sub's delay to get the relative value.
            float relative_delay_ms = (best_config_final.delay_s[s_idx] - min_ref_delay_s) * 1000.0f;
            std::cout << "  SW" << (s_idx + 1) << " Delay: " << std::setprecision(3) << relative_delay_ms << " ms, Polarity: "
                << (best_config_final.pol[s_idx] > 0.5f ? "Inverted" : "Normal") << std::endl;
        }
        //std::cout << "  Sum Lin Mags: " << std::scientific << std::setprecision(6) << overall_best_metric << std::defaultfloat << std::endl;
        if (num_freq_points_host > 0) {
            double avg_lin_mag = overall_best_metric / num_freq_points_host;
            if (avg_lin_mag > 1e-9) {
                std::cout << "  Avg SPL: " << std::fixed << std::setprecision(2) << 20.0 * log10(avg_lin_mag) << " dB" << std::endl;
            }
        }
    }
    else {
        std::cout << "\nNo valid optimal solution found." << std::endl;
    }
    std::cout << std::defaultfloat;
    CUDA_CHECK(cudaEventDestroy(start_ev)); CUDA_CHECK(cudaEventDestroy(stop_ev));
    if (d_frequencies_gpu) CUDA_CHECK(cudaFree(d_frequencies_gpu));
    if (d_sub_mags_flat_gpu) CUDA_CHECK(cudaFree(d_sub_mags_flat_gpu));
    if (d_sub_phases_flat_gpu) CUDA_CHECK(cudaFree(d_sub_phases_flat_gpu));
    if (d_delay_steps_gpu) CUDA_CHECK(cudaFree(d_delay_steps_gpu));
    if (d_block_results_gpu) CUDA_CHECK(cudaFree(d_block_results_gpu));
    CUDA_CHECK(cudaDeviceReset());
#ifdef _WIN32
    std::cout << "\nPress Enter to exit." << std::endl;
    std::cin.get();
#endif
    return 0;
}