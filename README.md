Multi-Subwoofer Time Delay & Polarity Optimizer (CUDA | Windows | NVIDIA GPU)
A GPU-accelerated brute-force optimizer for aligning multiple subwoofers via relative time delay and polarity (inversion). This tool ensures optimal results by exhaustively evaluating all possible combinations within the specified delay range.

🚀 Key Features
Brute-force precision: Exhaustive search guarantees the best alignment—no approximations or shortcuts.

GPU-accelerated: Optimized for memory coalescing and parallel execution on NVIDIA GPUs using CUDA.

Blazing fast: Near-instant results for up to 3 subs. With 4 subs, it can power through up to 1 trillion iterations in just a few hours on an RTX 3090 (Compute Capability 8.6).

Flexible resolution: Use fine delay steps (e.g. 0.1 ms) or coarser ones for faster results.

Human-hearing-aware optimization: Uses linearly weighted magnitude sums that approximate human hearing when REW exports use the PPO scale.

📦 Requirements
Windows PC

NVIDIA GPU (Compute Capability 5.2 / Maxwell or newer — RTX 3000+ recommended)

CUDA 12.9

Visual Studio 2022 (for building from source) or use the included .exe

Subwoofer measurements exported from REW as .txt files

Tab- or comma-separated formats supported

Use 96 PPO for best results

All subwoofers must be measured at the same mic position with a timing reference

🛠️ Build Instructions
To build from source:

Open the project in Visual Studio 2022

Set CUDA to version 12.9

Enable "Use Fast Math" in CUDA build settings

Build and run!

Or, just use the precompiled binary:
_win64_RTX3000+seriesGPU.exe

📁 Input Files
Export each subwoofer's frequency and phase response from REW (.txt format). Measurements should be aligned to a timing reference and taken from the same mic position.

👨‍💻 About
This tool uses a deliberately simple (and "dumb") brute-force method, exploiting modern GPU speed to guarantee optimal alignment. While not elegant, it's effective and fast—especially on high-end NVIDIA GPUs.

Sample screen output:

--- Multiple Subwoofer (max. 4) Delay & Inversion Optimizer (GPU) ---

--- Configuration ---
Enter desired parameters or press Enter to use defaults.
Minimum frequency to analyze (Hz) (default: 15): 20
Maximum frequency to analyze (Hz) (default: 250): 250

Select Measurement File for Subwoofer 1 (or Cancel to finish)
Loaded SW1: C:\Users\Ronin\Desktop\SW1o.txt (350pts)

Select Measurement File for Subwoofer 2 (or Cancel to finish)
Loaded SW2: C:\Users\Ronin\Desktop\SW2o.txt (350pts)

Select Measurement File for Subwoofer 3 (or Cancel to finish)
Loaded SW3: C:\Users\Ronin\Desktop\SW3o.txt (350pts)

Select Measurement File for Subwoofer 4 (or Cancel to finish)
No more files selected. Proceeding with 3 subwoofer(s).

All subwoofers use a common frequency axis.
  Frequency Range: 20.0 Hz to 248.5 Hz
  Frequency Points: 350

Minimum delay (ms) (default: -20.0):
Maximum delay (ms) (default: 20.0):
Delay step (ms) (default: 0.1): 0.05
Maximum delay spread between any two subs (ms) (default: 20.0):

Using 801 delay steps from -2e+01ms to 2e+01ms with 0.05ms step.
Using GPU: NVIDIA GeForce RTX 3090
Total combinations to check: 4,111,379,208
Launching kernel: 2624 blocks, 512 threads/block.

--- 3-Subwoofer Optimization Complete ---
Host Wall Time: 3.304 s
GPU Kernel Time: 3.304 s

Best Configuration Found:
  SW1 Delay: 15.250 ms, Polarity: Normal
  SW2 Delay: 15.650 ms, Polarity: Normal
  SW3 Delay: 19.000 ms, Polarity: Normal
  Avg SPL: 81.14 dB

Press Enter to exit.



