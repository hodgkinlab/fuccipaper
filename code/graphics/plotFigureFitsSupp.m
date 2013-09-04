%{
% Scatter plots of [timeToDiv] VS [grnOnEnd] for
% the 4 conditions, with a straight line fit through
% the origin and allowing a non-zero intercept.
% (supplementary figure)
%
% Andrey Kan + Hodgkin Lab
% akan@wehi.edu.au
%
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function plotFigureFitsSupp( ...
	mCellDataA, mCellDataB, mCellDataX, ...
	mCycleMeasA, mCycleMeasB, mCycleMeasX, par)

s = RandStream('mt19937ar', 'Seed', 100 + 2 + par.iDataset);
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

fprintf('[timeToDiv] VS [grnOnEnd], allowing intercept\n');
mVar1 = [[vTimeToDivA; vTimeToDivX], [vTimeToDivB; vNaN]];
mVar2 = [[vGrnOnEndA; vGrnOnEndX], [vGrnOnEndB; vNaN]];
stats = fitLinearModel(mVar1, mVar2, par);

v1 = mVar1(:);
v2 = mVar2(:);
vValid = (~isnan(v1) & ~isnan(v2));
v1 = v1(vValid);
v2 = v2(vValid);

[r, ~, rlo, rup] = corrcoef(v1(:), v2(:), 'rows', 'complete');
r = r(1,2); rlo = rlo(1,2); rup = rup(1,2);
fprintf('r = %0.2f (%0.2f—%0.2f)\n', r, rlo, rup);

axes(par.hAxes);
scatter(v1, v2, '.');
xlim([0, par.xMax]);
ylim([0, par.yMax]);

vLims = [0, par.xMax];
xMean = nanmean(v1);
yMean = nanmean(v2);
vCoeffsLow = stats.mCI(1,:);
vCoeffsHi = stats.mCI(2,:);

hold on;
line(vLims, (stats.vCoeffs(1).*vLims + stats.vCoeffs(2)), 'Color', 'b');

hold on;
line(vLims, (vCoeffsLow(1).*vLims - vCoeffsLow(1)*xMean + yMean), ...
	'color', 'r', 'linestyle', ':', 'LineWidth', 0.8);
hold on;
line(vLims, (vCoeffsHi(1).*vLims - vCoeffsHi(1)*xMean + yMean), ...
	'color', 'r', 'linestyle', ':', 'LineWidth', 0.8);

if (stats.vCoeffs(2) >= 0)
	cSgn = '+';
else
	cSgn = '-';
end
sText2 = sprintf('y = %.02f * x %s %.02f', ...
	stats.vCoeffs(1), cSgn, abs(stats.vCoeffs(2)));
sText3 = sprintf('(%.02f; %.02f)', vCoeffsLow(1), vCoeffsHi(1));
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

drawCircles(stats.vCoeffs(1), par);
end

function stats = fitLinearModel(mVar1, mVar2, par)
if (all(isnan(mVar1(:))) || all(isnan(mVar2(:))))
	fprintf('no data for one of the variables\n');
	return
end

% parameter estimates
vCoeffs = fitLinearModelSample(mVar1, mVar2);
disp(vCoeffs);

% 95% confidence intervals
mCI = bootci(par.nSamples, @fitLinearModelSample, mVar1, mVar2);
% disp(mCI);

stats.vCoeffs = vCoeffs;
stats.mCI = mCI;
end

function vCoeffs = fitLinearModelSample(mVar1, mVar2)
%{
	% sanity check
	vValid = ~isnan(mVar1(:,2));
	v1 = mVar1(vValid,1);
	v2 = mVar1(vValid,2);
	r = corr(v1, v2);
	fprintf('%.3f\n', r);
	v1 = mVar2(vValid,1);
	v2 = mVar2(vValid,2);
	r = corr(v1, v2);
	fprintf('%.3f\n', r);
%}
v1 = mVar1(:);
v2 = mVar2(:);
vValid = (~isnan(v1) & ~isnan(v2));
v1 = v1(vValid);
v2 = v2(vValid);

v0 = unique(v1);
if (numel(v0) < 3)
	vCoeffs = [NaN NaN];
	return
end

vCoeffs = lscov([v1, ones(size(v1))], v2, 1./(v1.^2));
% vCoeffs = flipud(deming(v1, v2, var(v2)/var(v1)));
end