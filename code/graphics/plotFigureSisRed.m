%{
% Scatter plots of time from division to red off for
% siblings for the 4 conditions.
% (supplementary figure)
%
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function plotFigureSisRed(mCycleMeasA, mCycleMeasB, par)
s = RandStream('mt19937ar', 'Seed', 100 + 4 + par.iDataset);
RandStream.setGlobalStream(s);
idxCycle = initTableCycleMeas();

% variables
if (1)
	vRedA = mCycleMeasA(:,idxCycle.timeRedOff);
	vRedB = mCycleMeasB(:,idxCycle.timeRedOff);
end
par.yMax = 25;
% xOffs = 5;

axes(par.hAxes);
scatter(vRedA, vRedB, '.');
xlim([0, par.yMax]);
ylim([0, par.yMax]);
[r, ~, rlo, rup] = corrcoef(vRedA, vRedB, 'rows', 'complete');
r = r(1,2); rlo = rlo(1,2); rup = rup(1,2);
fprintf('RED: r = %0.2f (%0.2f - %0.2f)\n', r, rlo, rup);

n = sum(~isnan(vRedA) & ~isnan(vRedB));
%text(1.5, 13, ...
%	[par.sTitle, ' (N = ', num2str(n), ')']);	
text(1.5, 15, sprintf('r = %.02f', r));
text(1.5, 12, sprintf('(%.02f—%.02f)', rlo, rup));
% text(-2, 16.5, sprintf('\\bf{%c}', 'A' + par.iDataset - 1), ...
% 	'FontSize', par.letterSize);
end