%{
% Scatter plots of [timeToDiv] VS [grnOnEnd] for
% the 4 conditions, with a straight line fit through the origin
% with the equation on the graph, an a pie chart in the other corner.
%
% Andrey Kan, et al.
% akan@wehi.edu.au
%
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function plotFigureFits( ...
	mCellDataA, mCellDataB, mCellDataX, ...
	mCycleMeasA, mCycleMeasB, mCycleMeasX, par)

s = RandStream('mt19937ar', 'Seed', 2 + par.iDataset);
RandStream.setGlobalStream(s);

idxData = initTableCellData();
idxCycle = initTableCycleMeas();

% variables
if (1)
	vTimeToDivA = mCellDataA(:,idxData.stop) - mCellDataA(:,idxData.start);
	vTimeToDivB = mCellDataB(:,idxData.stop) - mCellDataB(:,idxData.start);
	vTimeToDivX = mCellDataX(:,idxData.stop) - mCellDataX(:,idxData.start);
	
	vGrnOnEndA = vTimeToDivA - mCycleMeasA(:,idxCycle.timeGrnOn);
	vGrnOnEndB = vTimeToDivB - mCycleMeasB(:,idxCycle.timeGrnOn);
	vGrnOnEndX = vTimeToDivX - mCycleMeasX(:,idxCycle.timeGrnOn);
	
	vNaN = nan(size(vTimeToDivX));
end

fprintf('[timeToDiv] VS [grnOnEnd], through origin\n');
mVar1 = [[vTimeToDivA; vTimeToDivX], [vTimeToDivB; vNaN]];
mVar2 = [[vGrnOnEndA; vGrnOnEndX], [vGrnOnEndB; vNaN]];

stats = fitLineOrigin(mVar1, mVar2, par);

v1 = mVar1(:);
v2 = mVar2(:);
vValid = (~isnan(v1) & ~isnan(v2));
v1 = v1(vValid);
v2 = v2(vValid);
fprintf('mean time to div = %.02f\n', mean(v1));
fprintf('std time to div = %.02f\n', std(v1));

[r, ~, rlo, rup] = corrcoef(v1(:), v2(:), 'rows', 'complete');
r = r(1,2); rlo = rlo(1,2); rup = rup(1,2);
fprintf('model: r = %0.2f (%0.2f—%0.2f)\n', r, rlo, rup);

axes(par.hAxes);
scatter(v1, v2, '.');
xlim([0, par.xMax]);
ylim([0, par.yMax]);

line([0, par.xMax], [0, stats.beta*par.xMax], 'Color', 'b'); hold on;
line([0, par.xMax], [0, stats.vCI(1)*par.xMax], ...
	'Color', 'r', 'LineStyle', '--', 'LineWidth', 0.8); hold on;
line([0, par.xMax], [0, stats.vCI(2)*par.xMax], ...
	'Color', 'r', 'LineStyle', '--', 'LineWidth', 0.8); hold on;

sText2 = sprintf('y = %.02f * x', stats.beta);
sText3 = sprintf('(%.02f; %.02f)', stats.vCI(1), stats.vCI(2));
sText4 = sprintf('r = %0.2f', r);

text(par.text.x, par.text.y1, ...
	[par.sTitle, ' (N = ', num2str(numel(v1)), ')']);
text(par.text.x, par.text.y2, sText2);
text(par.text.x, par.text.y3, sText3);
text(par.text.x, par.text.y4, sText4);

text(25.6, 9.5, 'G1');
text(17, 3.6, 'S/G2/M');
text(-2.8, 21.5, sprintf('\\bf{%c}', 'A' + par.iDataset - 1), ...
	'FontSize', par.letterSize);

drawCircles(stats.beta, par);
end

function stats = fitLineOrigin(mVar1, mVar2, par)
if (all(isnan(mVar1(:))) || all(isnan(mVar2(:))))
	fprintf('no data for one of the variables\n');
	return
end

% parameter estimates
beta = fitLineOriginSample(mVar1, mVar2);
disp(beta);

% 95% confidence intervals
vCI = bootci(par.nSamples, @fitLineOriginSample, mVar1, mVar2);
disp(vCI);

stats.beta = beta;
stats.vCI = vCI;
end

function beta = fitLineOriginSample(mVar1, mVar2)
v1 = mVar1(:);
v2 = mVar2(:);
vValid = (~isnan(v1) & ~isnan(v2));
v1 = v1(vValid);
v2 = v2(vValid);

v0 = unique(v1);
if (numel(v0) < 3)
	beta = NaN;
	return
end

beta = sum(v1.*v2)/sum(v1.*v1);
end