out = experimentGapSweep([12 16 20 24 28 32], 20, true);   % adds the mission-energy gap column
% plotGapSweep(out, 'energy');                               % energy gap vs L and vs diameter
% plotGapSweep(out, 'energy', 'Layout', 'stacked', 'Save', true, 'File', 'GapSweepEnergy');
plotGapSweep(out, 'energy', 'Layout', 'stacked');
