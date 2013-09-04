%{
(A-D) Scatter plots of sibling division times for the 4
conditions. (E-H) Scatter plots of time green for the siblings.
%
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function plotFigureSis( ...
	mCellDataA, mCellDataB, mCycleMeasA, mCycleMeasB, par)
s = RandStream('mt19937ar', 'Seed', 4 + par.iDataset);
RandStream.setGlobalStream(s);

idxData = initTableCellData();
idxCycle = initTableCycleMeas();

% variables
if (1)
	vTimeToDivA = mCellDataA(:,idxData.stop) - mCellDataA(:,idxData.start);
	vTimeToDivB = mCellDataB(:,idxData.stop) - mCellDataB(:,idxData.start);	
	vGrnOnEndA = vTimeToDivA - mCycleMeasA(:,idxCycle.timeGrnOn);
	vGrnOnEndB = vTimeToDivB - mCycleMeasB(:,idxCycle.timeGrnOn);
end

par.xMax = 25;
xOffs = 0.65;
axes(par.hAxesTop);
scatter(vTimeToDivA, vTimeToDivB, '.');
xlim([xOffs, par.xMax]);
ylim([xOffs, par.xMax]);
[r, ~, rlo, rup] = corrcoef(vTimeToDivA, vTimeToDivB);
r = r(1,2); rlo = rlo(1,2); rup = rup(1,2);
fprintf('TOT: r = %0.2f (%0.2f - %0.2f)\n', r, rlo, rup);
text(par.text.x + xOffs, 27, ...
	[par.sTitle, ' (N = ', num2str(numel(vTimeToDivA)), ')']);
text(par.text.x + xOffs, 25, sprintf('r = %.02f', r));
text(-4, 27, sprintf('\\bf{%c}', 'A' + par.iDataset - 1), ...
	'FontSize', par.letterSize);

xOffs = -2;
axes(par.hAxesBot);
scatter(vGrnOnEndA, vGrnOnEndB, '.');
xlim([0, par.xMax]);
ylim([0, par.xMax]);
[r, ~, rlo, rup] = corrcoef(vGrnOnEndA, vGrnOnEndB, 'rows', 'complete');
r = r(1,2); rlo = rlo(1,2); rup = rup(1,2);
fprintf('GRN: r = %0.2f (%0.2f - %0.2f)\n', r, rlo, rup);
text(par.text.x + xOffs, 27, sprintf('r = %.02f', r));
text(-3, 27, sprintf('\\bf{%c}', 'E' + par.iDataset - 1), ...
	'FontSize', par.letterSize);
end