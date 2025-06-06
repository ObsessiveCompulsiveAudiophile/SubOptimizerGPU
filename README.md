# SubOptimizerGPU
Multiple subwoofer relative time delay and inversion optimizer (Windows only, requires NVDIA GPU)

Ultrafast (memory coalescing optimized GPU parallelization) brute force (results are optimal) muiltiple subwoofer time delay and inversion optimizer. It's almost instant for any time delay step for up to 3 subs. With 4 subs it can chew through a trillion iterations in a couple of hours (RTX3090 compute 86) but you can use 0.1ms steps or larger for much quicker results. It needs frequency and phase responses of all subs as .txt files. Comma or tab separated REW .txt exports (prefer 96 ppo) will work. Subwoofers should have been measured with timing reference at the same mic position.

Compile with "Use Fast Math" (or run the .exe file directly - Nvdia RTX3000 series or higher)

Export subwoofer measurements as .txt from REW (use 96 ppo fro optimal results)

Written in VS Studio 2022 with Cuda 12.9

OCA


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



