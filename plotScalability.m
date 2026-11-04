function h = plotScalability(out, opts)
%PLOTSCALABILITY  2x2 scalability panel: C-MAS vs D-MAS over team size N.
%
%   h = PLOTSCALABILITY(out) draws the four panels of Fig.~\ref{fig:ScalabilityPlot}
%   from EXPERIMENTSCALABILITY output: (a) per-cycle runtime, (b) convergence time,
%   (c) communication overhead, (d) mission energy, each versus N with C-MAS (blue)
%   and D-MAS (red).
%
%   Name-Value options (as in PLOTINSTANCE): 'Save' (false), 'File'
%   ('ScalabilityPlot1'), 'Format' ('png'), 'Resolution' (300).
%
%   See also EXPERIMENTSCALABILITY.

arguments
    out              (1,1) struct
    opts.Save        (1,1) logical = false
    opts.File        (1,:) char    = 'ScalabilityPlot1'
    opts.Format      (1,:) char    = 'png'
    opts.Resolution  (1,1) double  = 300
end

cC = [0 0.447 0.741];      % C-MAS blue
cD = [0.85 0.325 0.098];   % D-MAS red
N  = out.N;

h = figure('Name','Scalability: C-MAS vs D-MAS','Color','w','Position',[60 60 980 720]);
t = tiledlayout(h, 2, 2, 'TileSpacing','compact', 'Padding','compact');

% (a) per-cycle runtime (log-y: spans small values, D-MAS << C-MAS at scale)
ax = nexttile(t); hold(ax,'on'); box(ax,'on'); grid(ax,'on'); set(ax,'YScale','log');
plot(ax, N, out.tCMAS, '-o', 'Color',cC, 'MarkerFaceColor',cC, 'LineWidth',1.7);
plot(ax, N, out.tDMAS, '-s', 'Color',cD, 'MarkerFaceColor',cD, 'LineWidth',1.7);
xlabel(ax,'Number of agents $N$','Interpreter','latex');
ylabel(ax,'Per-cycle runtime (ms)','Interpreter','latex');
title(ax,'(a) Algorithm runtime','Interpreter','latex');
legend(ax, {'\texttt{C-MAS}','\texttt{D-MAS}'}, 'Interpreter','latex','Location','northwest');

% (b) convergence time
ax = nexttile(t); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, N, out.TcC, '-o', 'Color',cC, 'MarkerFaceColor',cC, 'LineWidth',1.7);
plot(ax, N, out.TcD, '-s', 'Color',cD, 'MarkerFaceColor',cD, 'LineWidth',1.7);
xlabel(ax,'Number of agents $N$','Interpreter','latex');
ylabel(ax,'Convergence time $T_c$ (s)','Interpreter','latex');
title(ax,'(b) Convergence time','Interpreter','latex');
legend(ax, {'\texttt{C-MAS}','\texttt{D-MAS}'}, 'Interpreter','latex','Location','northwest');

% (c) communication overhead (log-y: D-MAS grows O(N^2) vs O(N))
ax = nexttile(t); hold(ax,'on'); box(ax,'on'); grid(ax,'on'); set(ax,'YScale','log');
plot(ax, N, out.ovC, '-o', 'Color',cC, 'MarkerFaceColor',cC, 'LineWidth',1.7);
plot(ax, N, out.ovD, '-s', 'Color',cD, 'MarkerFaceColor',cD, 'LineWidth',1.7);
xlabel(ax,'Number of agents $N$','Interpreter','latex');
ylabel(ax,'Communication overhead (Kbits)','Interpreter','latex');
title(ax,'(c) Communication overhead','Interpreter','latex');
legend(ax, {'\texttt{C-MAS}','\texttt{D-MAS}'}, 'Interpreter','latex','Location','northwest');

% (d) energy
ax = nexttile(t); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, N, out.eC, '-o', 'Color',cC, 'MarkerFaceColor',cC, 'LineWidth',1.7);
plot(ax, N, out.eD, '-s', 'Color',cD, 'MarkerFaceColor',cD, 'LineWidth',1.7);
xlabel(ax,'Number of agents $N$','Interpreter','latex');
ylabel(ax,'Mission energy (mAh)','Interpreter','latex');
title(ax,'(d) Energy consumption','Interpreter','latex');
legend(ax, {'\texttt{C-MAS}','\texttt{D-MAS}'}, 'Interpreter','latex','Location','northwest');

title(t, 'Scalability versus team size $N$', 'Interpreter','latex', 'FontWeight','bold');

if opts.Save
    fname = [opts.File '.' opts.Format];
    if ismember(lower(opts.Format), {'pdf','eps'})
        exportgraphics(h, fname, 'ContentType','vector');
    else
        exportgraphics(h, fname, 'Resolution', opts.Resolution);
    end
    fprintf('Scalability figure saved as %s\n', fname);
end
end
