results.CMAS.algorithm_time  = [3.204 3.342 3.545 3.662 3.774];
results.CMAS.convergence_time= [49.2 69.9 93.0 107.8 133.8];
results.CMAS.overhead        = [5168 17783 46934 108202 268168];
results.CMAS.energy          = [1824 4460 8792 17538 35124];
results.DMAS.algorithm_time  = [2.910 2.813 2.795 2.758 2.844];
results.DMAS.convergence_time= [48.2 66.9 87.7 106.8 132.2];
results.DMAS.overhead        = [63550 458133 2391743 9711007 42365727];
results.DMAS.energy          = [1827 4473 8839 17602 35257];

% map into the struct plotScalability expects
out.N     = [10 25 50 100 200];
out.tCMAS = results.CMAS.algorithm_time;    out.tDMAS = results.DMAS.algorithm_time;
out.TcC   = results.CMAS.convergence_time;  out.TcD   = results.DMAS.convergence_time;
out.ovC   = results.CMAS.overhead;          out.ovD   = results.DMAS.overhead;
out.eC    = results.CMAS.energy;            out.eD    = results.DMAS.energy;

plotScalability(out, 'Save', true, 'File', 'ScalabilityPlot', 'Format', 'png');
