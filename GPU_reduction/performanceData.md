# Performance comparison of various reductions and timestepping implementations on different hardware and systems

## Results from running on a V100 on C-3PO

### Run 1

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

### Run 2

prop.maxThreadsPerMultiProcessor = 2048
prop.multiProcessorCount         = 80
prop.maxThreadsPerBlock          = 1024

numThreads = 163840
numBlocks  = 160

Timer name: AtomicMax Reduction Timer
  Number of trials: 1000, Total time: 1.28675s, Average Time: 1.28675ms, Standard Deviation: 6.72352µs, Fastest Run: 1.26728ms, Slowest Run: 1.31333ms

Timer name: Max Reduction Timer
  Number of trials: 1000, Total time: 1.29077s, Average Time: 1.29077ms, Standard Deviation: 8.80753µs, Fastest Run: 1.27375ms, Slowest Run: 1.45683ms

Timer name: DTI Original Timer
  Number of trials: 1000, Total time: 7.97112s, Average Time: 7.97112ms, Standard Deviation: 55.5687µs, Fastest Run: 7.69729ms, Slowest Run: 8.5001ms

Timer name: DTI New Timer
  Number of trials: 1000, Total time: 386.332ms, Average Time: 386.332µs, Standard Deviation: 209.528µs, Fastest Run: 357.445µs, Slowest Run: 2.42069ms

## Results from running on a V100 on CRC H2P PPC-n0

### Run 1

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

### Run 2

prop.maxThreadsPerMultiProcessor = 2048
prop.multiProcessorCount         = 80
prop.maxThreadsPerBlock          = 1024

numThreads = 163840
numBlocks  = 160

Timer name: AtomicMax Reduction Timer
  Number of trials: 1000, Total time: 1.30014s, Average Time: 1.30014ms, Standard Deviation: 8.39229µs, Fastest Run: 1.28702ms, Slowest Run: 1.48747ms

Timer name: Max Reduction Timer
  Number of trials: 1000, Total time: 1.3128s, Average Time: 1.3128ms, Standard Deviation: 8.5439µs, Fastest Run: 1.28892ms, Slowest Run: 1.33816ms

Timer name: DTI Original Timer
  Number of trials: 1000, Total time: 7.66384s, Average Time: 7.66384ms, Standard Deviation: 31.8569µs, Fastest Run: 7.44163ms, Slowest Run: 7.73476ms

Timer name: DTI New Timer
  Number of trials: 1000, Total time: 485.294ms, Average Time: 485.294µs, Standard Deviation: 157.341µs, Fastest Run: 472.819µs, Slowest Run: 4.00286ms

## Results from running on a single GCD on Crusher

### Run 1

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

### Run 2

prop.maxThreadsPerMultiProcessor = 2048
prop.multiProcessorCount         = 110
prop.maxThreadsPerBlock          = 1024

numThreads = 225280
numBlocks  = 220

Timer name: AtomicMax Reduction Timer
  Number of trials: 1000, Total time: 806.584ms, Average Time: 806.584µs, Standard Deviation: 2.11728µs, Fastest Run: 802.007µs, Slowest Run: 813.42µs

Timer name: Max Reduction Timer
  Number of trials: 1000, Total time: 810.639ms, Average Time: 810.639µs, Standard Deviation: 2.12174µs, Fastest Run: 806.025µs, Slowest Run: 818.368µs

Timer name: DTI Original Timer
  Number of trials: 1000, Total time: 5.77731s, Average Time: 5.77731ms, Standard Deviation: 2.84564µs, Fastest Run: 5.7694ms, Slowest Run: 5.78748ms

Timer name: DTI New Timer
  Number of trials: 1000, Total time: 55.8442ms, Average Time: 55.8442µs, Standard Deviation: 1.26596µs, Fastest Run: 55.286µs, Slowest Run: 84.934µs
