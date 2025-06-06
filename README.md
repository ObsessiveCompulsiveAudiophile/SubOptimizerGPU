# SubOptimizerGPU
Multiple subwoofer relative time delay and inversion optimizer (Windows only, requires NVDIA GPU)

Ultrafast (memory coalescing optimized GPU parallelization) brute force (results are optimal) muiltiple subwoofer time delay and inversion optimizer. It's almost instant for any time delay step for up to 3 subs. With 4 subs it can chew through a trillion iterations in a couple of hours (RTX3090 compute 86) but you can use 0.1ms steps or larger for much quicker results. It needs frequency and phase responses of all subs as .txt files. Comma or tab separated REW .txt exports (prefer 96 ppo) will work. Subwoofers should have been measured with timing reference at the same mic position.

Compile with "Use Fast Math" (or run the .exe file directly - RTXC3000 series or higher)

Export subwoofer measurements as .txt from REW (use 96 ppo fro optimal results)

Written in VS Studio 2022 with Cuda 12.9

OCA




