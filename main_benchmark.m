out = experimentBenchmark([4 8 16 32 64 96], 50);
plotBenchmark(out, 'Save', true, 'File', 'Benchmark', 'Format', 'png');
