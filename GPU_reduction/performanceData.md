# Performance comparison of various reductions and timestepping implementations on different hardware and systems

## Summary Tables

All run statistics are the result of 1,000 runs with 5 runs for warmup
beforehand. The "atomic" vs. non-atomic reductions differ only in how the final
grid reduction is done. In the "atomic" case it's done with atomic operations
and the non-atomic case uses a second kernel launch

### 512^3 Grid Test

| System: C-3PO/V100  | Average Time | Standard Deviation | Fastest Run | Slowest Run | Speedup |
|---------------------|--------------|--------------------|-------------|-------------|---------|
| Atomic Max          |  |  |  |  |  |
| Non-Atomic Max      |  |  |  |  |  |
| Original DTI Kernel |  |  |  |  |  |
| New DTI Kernel      |  |  |  |  |  |

## Results from running on a V100 on C-3PO

### 256^3 Run - C-3PO

nvcc -g -O3 -std=c++11 -arch sm_70 -fmad=false -I/usr/local/cuda-11.4/include -DCUDA_BUILD gpu_reductions.cu -o gpu_reductions.exe
prop.maxThreadsPerMultiProcessor = 2048
prop.multiProcessorCount         = 80
prop.maxThreadsPerBlock          = 1024

numThreads = 163840
numBlocks  = 160

number of trials = 1000
Grid size        = 256

Timer name: AtomicMax Reduction Timer
  Number of trials: 1000, Total time: 173.855ms, Average Time: 173.855µs, Standard Deviation: 2.54315µs, Fastest Run: 168.593µs, Slowest Run: 182.547µs

Timer name: Max Reduction Timer
  Number of trials: 1000, Total time: 177.751ms, Average Time: 177.751µs, Standard Deviation: 3.01136µs, Fastest Run: 171.542µs, Slowest Run: 186.331µs

Timer name: DTI Original Timer
  Number of trials: 1000, Total time: 1.22304s, Average Time: 1.22304ms, Standard Deviation: 31.1317µs, Fastest Run: 1.19895ms, Slowest Run: 1.38855ms

Timer name: DTI New Timer
  Number of trials: 1000, Total time: 1.41235s, Average Time: 1.41235ms, Standard Deviation: 421.304µs, Fastest Run: 1.27623ms, Slowest Run: 3.35978ms

maxReducedAtomic = 4
maxReduced       = 4
oldDTI           = 1.14468
newDTI           = 1.14468

### 512^3 Run - C-3PO

prop.maxThreadsPerMultiProcessor = 2048
prop.multiProcessorCount         = 80
prop.maxThreadsPerBlock          = 1024

numThreads = 163840
numBlocks  = 160

number of trials = 1000
Grid size        = 512

Timer name: AtomicMax Reduction Timer
  Number of trials: 1000, Total time: 1.28774s, Average Time: 1.28774ms, Standard Deviation: 6.69806µs, Fastest Run: 1.26784ms, Slowest Run: 1.31017ms

Timer name: Max Reduction Timer
  Number of trials: 1000, Total time: 1.29181s, Average Time: 1.29181ms, Standard Deviation: 8.73028µs, Fastest Run: 1.27522ms, Slowest Run: 1.44953ms

Timer name: DTI Original Timer
  Number of trials: 1000, Total time: 9.35504s, Average Time: 9.35504ms, Standard Deviation: 49.2086µs, Fastest Run: 9.22509ms, Slowest Run: 9.59188ms

Timer name: DTI New Timer
  Number of trials: 1000, Total time: 7.69866s, Average Time: 7.69866ms, Standard Deviation: 463.887µs, Fastest Run: 7.53224ms, Slowest Run: 10.6045ms

maxReducedAtomic = 4
maxReduced       = 4
oldDTI           = 1.15467
newDTI           = 1.15467

## Results from running on a V100 on CRC H2P PPC-n0

### 256^3 Run - CRC

prop.maxThreadsPerMultiProcessor = 2048
prop.multiProcessorCount         = 80
prop.maxThreadsPerBlock          = 1024

numThreads = 163840
numBlocks  = 160

number of trials = 1000
Grid size        = 256

Timer name: AtomicMax Reduction Timer
  Number of trials: 1000, Total time: 187.349ms, Average Time: 187.349µs, Standard Deviation: 3.361µs, Fastest Run: 180.557µs, Slowest Run: 195.215µs

Timer name: Max Reduction Timer
  Number of trials: 1000, Total time: 190.728ms, Average Time: 190.728µs, Standard Deviation: 2.56455µs, Fastest Run: 184.612µs, Slowest Run: 199.431µs

Timer name: DTI Original Timer
  Number of trials: 1000, Total time: 1.2716s, Average Time: 1.2716ms, Standard Deviation: 13.4786µs, Fastest Run: 1.24724ms, Slowest Run: 1.63597ms

Timer name: DTI New Timer
  Number of trials: 1000, Total time: 1.43593s, Average Time: 1.43593ms, Standard Deviation: 143.489µs, Fastest Run: 1.41477ms, Slowest Run: 4.95555ms

maxReducedAtomic = 4
maxReduced       = 4
oldDTI           = 1.14468
newDTI           = 1.14468

### 512^3 Run - CRC

prop.maxThreadsPerMultiProcessor = 2048
prop.multiProcessorCount         = 80
prop.maxThreadsPerBlock          = 1024

numThreads = 163840
numBlocks  = 160

number of trials = 1000
Grid size        = 512

Timer name: AtomicMax Reduction Timer
  Number of trials: 1000, Total time: 1.30471s, Average Time: 1.30471ms, Standard Deviation: 9.93422µs, Fastest Run: 1.2827ms, Slowest Run: 1.42029ms

Timer name: Max Reduction Timer
  Number of trials: 1000, Total time: 1.30728s, Average Time: 1.30728ms, Standard Deviation: 11.5186µs, Fastest Run: 1.2899ms, Slowest Run: 1.5356ms

Timer name: DTI Original Timer
  Number of trials: 1000, Total time: 9.58485s, Average Time: 9.58485ms, Standard Deviation: 24.7759µs, Fastest Run: 9.52506ms, Slowest Run: 9.83532ms

Timer name: DTI New Timer
  Number of trials: 1000, Total time: 7.70318s, Average Time: 7.70318ms, Standard Deviation: 124.398µs, Fastest Run: 7.67089ms, Slowest Run: 10.2437ms

maxReducedAtomic = 4
maxReduced       = 4
oldDTI           = 1.15467
newDTI           = 1.15467

## Results from running on a single GCD on Crusher

### 256^3 Run - Crusher

prop.maxThreadsPerMultiProcessor = 2048
prop.multiProcessorCount         = 110
prop.maxThreadsPerBlock          = 1024

numThreads = 225280
numBlocks  = 220

number of trials = 1000
Grid size        = 256

Timer name: AtomicMax Reduction Timer
  Number of trials: 1000, Total time: 143.947ms, Average Time: 143.947µs, Standard Deviation: 525.102ns, Fastest Run: 142.484µs, Slowest Run: 150.629µs

Timer name: Max Reduction Timer
  Number of trials: 1000, Total time: 147.132ms, Average Time: 147.132µs, Standard Deviation: 515.63ns, Fastest Run: 145.881µs, Slowest Run: 151.1µs

Timer name: DTI Original Timer
  Number of trials: 1000, Total time: 791.041ms, Average Time: 791.041µs, Standard Deviation: 1.44385µs, Fastest Run: 786.921µs, Slowest Run: 812.711µs

Timer name: DTI New Timer
  Number of trials: 1000, Total time: 724.83ms, Average Time: 724.83µs, Standard Deviation: 19.7034µs, Fastest Run: 684.234µs, Slowest Run: 820.155µs

maxReducedAtomic = 4
maxReduced       = 4
oldDTI           = 1.14468
newDTI           = 1.14468

### 512^3 Run - Crusher

prop.maxThreadsPerMultiProcessor = 2048
prop.multiProcessorCount         = 110
prop.maxThreadsPerBlock          = 1024

numThreads = 225280
numBlocks  = 220

number of trials = 1000
Grid size        = 512

Timer name: AtomicMax Reduction Timer
  Number of trials: 1000, Total time: 804.097ms, Average Time: 804.097µs, Standard Deviation: 2.19387µs, Fastest Run: 799.716µs, Slowest Run: 816.929µs

Timer name: Max Reduction Timer
  Number of trials: 1000, Total time: 807.347ms, Average Time: 807.347µs, Standard Deviation: 2.20808µs, Fastest Run: 802.992µs, Slowest Run: 816.167µs

Timer name: DTI Original Timer
  Number of trials: 1000, Total time: 5.9293s, Average Time: 5.9293ms, Standard Deviation: 11.1097µs, Fastest Run: 5.91888ms, Slowest Run: 6.00364ms

Timer name: DTI New Timer
  Number of trials: 1000, Total time: 6.80555s, Average Time: 6.80555ms, Standard Deviation: 203.338µs, Fastest Run: 6.27068ms, Slowest Run: 7.57096ms

maxReducedAtomic = 4
maxReduced       = 4
oldDTI           = 1.15467
newDTI           = 1.15467
