# Performance comparison of various reductions and timestepping implementations on different hardware and systems

## Results from running on a V100 on C-3PO

### 512^3 Run - C-3PO

prop.maxThreadsPerMultiProcessor = 2048
prop.multiProcessorCount         = 80
prop.maxThreadsPerBlock          = 1024

numThreads = 163840
numBlocks  = 160

Timer name: AtomicMax Reduction Timer
  Number of trials: 1000, Total time: 1.28748s, Average Time: 1.28748ms, Standard Deviation: 7.29302µs, Fastest Run: 1.27067ms, Slowest Run: 1.31923ms

Timer name: Max Reduction Timer
  Number of trials: 1000, Total time: 1.29094s, Average Time: 1.29094ms, Standard Deviation: 6.84559µs, Fastest Run: 1.27254ms, Slowest Run: 1.31207ms

Timer name: DTI Original Timer
  Number of trials: 1000, Total time: 7.9751s, Average Time: 7.9751ms, Standard Deviation: 72.4559µs, Fastest Run: 7.80672ms, Slowest Run: 8.50795ms

Timer name: DTI New Timer
  Number of trials: 1000, Total time: 382.29ms, Average Time: 382.29µs, Standard Deviation: 186.041µs, Fastest Run: 357.272µs, Slowest Run: 2.38161ms

### 256^3 Run - C-3PO

prop.maxThreadsPerMultiProcessor = 2048
prop.multiProcessorCount         = 80
prop.maxThreadsPerBlock          = 1024

numThreads = 163840
numBlocks  = 160

Timer name: AtomicMax Reduction Timer
  Number of trials: 1000, Total time: 177.18ms, Average Time: 177.18µs, Standard Deviation: 8.63067µs, Fastest Run: 169.306µs, Slowest Run: 239.635µs

Timer name: Max Reduction Timer
  Number of trials: 1000, Total time: 178.13ms, Average Time: 178.13µs, Standard Deviation: 3.03273µs, Fastest Run: 172.074µs, Slowest Run: 186.373µs

Timer name: DTI Original Timer
  Number of trials: 1000, Total time: 1.00592s, Average Time: 1.00592ms, Standard Deviation: 17.7278µs, Fastest Run: 950.118µs, Slowest Run: 1.0389ms

Timer name: DTI New Timer
  Number of trials: 1000, Total time: 388.19ms, Average Time: 388.19µs, Standard Deviation: 223.459µs, Fastest Run: 356.062µs, Slowest Run: 2.40895ms

## Results from running on a V100 on CRC H2P PPC-n0

### 512^3 Run - CRC

prop.maxThreadsPerMultiProcessor = 2048
prop.multiProcessorCount         = 80
prop.maxThreadsPerBlock          = 1024

numThreads = 163840
numBlocks  = 160

Timer name: AtomicMax Reduction Timer
  Number of trials: 1000, Total time: 1.30623s, Average Time: 1.30623ms, Standard Deviation: 10.4087µs, Fastest Run: 1.28777ms, Slowest Run: 1.42357ms

Timer name: Max Reduction Timer
  Number of trials: 1000, Total time: 1.30989s, Average Time: 1.30989ms, Standard Deviation: 11.4767µs, Fastest Run: 1.28761ms, Slowest Run: 1.54134ms

Timer name: DTI Original Timer
  Number of trials: 1000, Total time: 7.59033s, Average Time: 7.59033ms, Standard Deviation: 46.452µs, Fastest Run: 7.55999ms, Slowest Run: 7.85941ms

Timer name: DTI New Timer
  Number of trials: 1000, Total time: 476.758ms, Average Time: 476.758µs, Standard Deviation: 112.033µs, Fastest Run: 468.648µs, Slowest Run: 4.00577ms

### 256^3 Run - CRC

prop.maxThreadsPerMultiProcessor = 2048
prop.multiProcessorCount         = 80
prop.maxThreadsPerBlock          = 1024

numThreads = 163840
numBlocks  = 160

Timer name: AtomicMax Reduction Timer
  Number of trials: 1000, Total time: 188.855ms, Average Time: 188.855µs, Standard Deviation: 3.38397µs, Fastest Run: 182.645µs, Slowest Run: 195.791µs

Timer name: Max Reduction Timer
  Number of trials: 1000, Total time: 191.099ms, Average Time: 191.099µs, Standard Deviation: 2.87885µs, Fastest Run: 185.142µs, Slowest Run: 200.409µs

Timer name: DTI Original Timer
  Number of trials: 1000, Total time: 958.653ms, Average Time: 958.653µs, Standard Deviation: 4.61462µs, Fastest Run: 949.03µs, Slowest Run: 1.0472ms

Timer name: DTI New Timer
  Number of trials: 1000, Total time: 477.213ms, Average Time: 477.213µs, Standard Deviation: 115.592µs, Fastest Run: 467.828µs, Slowest Run: 3.97577ms

## Results from running on a single GCD on Crusher

### 512^3 Run - Crusher

prop.maxThreadsPerMultiProcessor = 2048
prop.multiProcessorCount         = 110
prop.maxThreadsPerBlock          = 1024

numThreads = 225280
numBlocks  = 220

Timer name: AtomicMax Reduction Timer
  Number of trials: 1000, Total time: 807.75ms, Average Time: 807.75µs, Standard Deviation: 2.42729µs, Fastest Run: 802.158µs, Slowest Run: 817.657µs

Timer name: Max Reduction Timer
  Number of trials: 1000, Total time: 811.217ms, Average Time: 811.217µs, Standard Deviation: 2.26458µs, Fastest Run: 805.855µs, Slowest Run: 820.062µs

Timer name: DTI Original Timer
  Number of trials: 1000, Total time: 5.77784s, Average Time: 5.77784ms, Standard Deviation: 9.1389µs, Fastest Run: 5.76928ms, Slowest Run: 6.03113ms

Timer name: DTI New Timer
  Number of trials: 1000, Total time: 55.7436ms, Average Time: 55.7436µs, Standard Deviation: 296.163ns, Fastest Run: 55.297µs, Slowest Run: 59.344µs

### 256^3 Run - Crusher

prop.maxThreadsPerMultiProcessor = 2048
prop.multiProcessorCount         = 110
prop.maxThreadsPerBlock          = 1024

numThreads = 225280
numBlocks  = 220

Timer name: AtomicMax Reduction Timer
  Number of trials: 1000, Total time: 143.293ms, Average Time: 143.293µs, Standard Deviation: 1.84799µs, Fastest Run: 142.064µs, Slowest Run: 192.089µs

Timer name: Max Reduction Timer
  Number of trials: 1000, Total time: 146.281ms, Average Time: 146.281µs, Standard Deviation: 397.09ns, Fastest Run: 145.349µs, Slowest Run: 148.817µs

Timer name: DTI Original Timer
  Number of trials: 1000, Total time: 742.353ms, Average Time: 742.353µs, Standard Deviation: 1.7847µs, Fastest Run: 738.871µs, Slowest Run: 780.781µs

Timer name: DTI New Timer
  Number of trials: 1000, Total time: 55.7243ms, Average Time: 55.7243µs, Standard Deviation: 3.38105µs, Fastest Run: 55.116µs, Slowest Run: 130.1µs
