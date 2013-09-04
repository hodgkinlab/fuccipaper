%{
Christmas tree plots of siblings for the four conditions. (A-D)
unscaled, (E-H) scaled with bars and histograms underneath to show the
average time of events and the distribution,
%
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function plotFigureBars( ...
	mCellDataA, mCellDataB, mCycleMeasA, mCycleMeasB, par)
s = RandStream('mt19937ar', 'Seed', 3 + par.iDataset);
RandStream.setGlobalStream(s);

par.xMax = 25;

idxData = initTableCellData();
vTimeToDivA = mCellDataA(:,idxData.stop) - mCellDataA(:,idxData.start);
vTimeToDivB = mCellDataB(:,idxData.stop) - mCellDataB(:,idxData.start);
[vTimeToDivA, vOrder] = sort(vTimeToDivA, 'descend');
vTimeToDivB = vTimeToDivB(vOrder);
mCycleMeasA = mCycleMeasA(vOrder,:);
mCycleMeasB = mCycleMeasB(vOrder,:);

lw = 0.9;
cl = [0.99, 1, 1];
nPairs = numel(vTimeToDivA);

% bars, normalized
axes(par.hAxes(2));
par.xSep = 0.005;
mCycleMeasNormA = mCycleMeasA ...
	./ repmat(vTimeToDivA, 1, size(mCycleMeasA, 2));
mCycleMeasNormB = mCycleMeasB ...
	./ repmat(vTimeToDivB, 1, size(mCycleMeasB, 2));
plotBars(ones(size(vTimeToDivA)), mCycleMeasNormA, par, -1);
plotBars(ones(size(vTimeToDivB)), mCycleMeasNormB, par, 1);
hold on;
line([0, 0], [0, numel(vTimeToDivA)], 'linewidth', lw, 'color', cl);
xlim([-1.05, 1.05]);
ylim([0, nPairs + 1]);

% bars, not normalized
axes(par.hAxes(1));
plotBars(vTimeToDivA, mCycleMeasA, par, -1);
plotBars(vTimeToDivB, mCycleMeasB, par, 1);

hold on;
line([0, 0], [0, nPairs], 'linewidth', lw, 'color', cl);
xlim([-par.xMax, par.xMax]);

yText = 0.9*nPairs;

hold on;
text(-33, yText, sprintf('\\bf{%c}', 'A' + par.iDataset - 1), ...
	'FontSize', par.letterSize);
hold on;
text(22, yText, sprintf('\\bf{%c}', 'E' + par.iDataset - 1), ...
	'FontSize', par.letterSize);

text(-29, 0.9*nPairs, par.sTitle);
text(-29, 0.74*nPairs, ['N = ', num2str(nPairs), ' pairs']);
ylim([0, nPairs + 1]);


yMax = 0.75;

% histograms, normalized
axes(par.hAxes(4));
par.xSep = 0.005;
par.nBins = 20;
mCycleMeasNormA = mCycleMeasA ...
	./ repmat(vTimeToDivA, 1, size(mCycleMeasA, 2));
mCycleMeasNormB = mCycleMeasB ...
	./ repmat(vTimeToDivB, 1, size(mCycleMeasB, 2));
plotHistograms(ones(size(vTimeToDivA)), mCycleMeasNormA, par, -1);
plotHistograms(ones(size(vTimeToDivB)), mCycleMeasNormB, par, 1);
hold on;
line([0, 0], [0, numel(vTimeToDivA)], 'LineWidth', lw, 'Color', cl);
xlim([-1.05, 1.05]);
ylim([-yMax, yMax]);

% histograms, not normalized
axes(par.hAxes(3));
par.nBins = 30;
plotHistograms(vTimeToDivA, mCycleMeasA, par, -1);
plotHistograms(vTimeToDivB, mCycleMeasB, par, 1);
hold on;
line([0, 0], [0, numel(vTimeToDivA)], 'linewidth', lw, 'color', cl);
xlim([-par.xMax, par.xMax]);
ylim([-yMax, yMax]);

end

function plotBars(vTimeToDiv, mCycleMeas, par, flipCoeff)
idxCycle = initTableCycleMeas();
vTimeRedOn	= mCycleMeas(:,idxCycle.timeRedOn);
vTimeRedOff = mCycleMeas(:,idxCycle.timeRedOff);
vTimeGrnOn	= mCycleMeas(:,idxCycle.timeGrnOn);
vTimeGrnOff = mCycleMeas(:,idxCycle.timeGrnOff);

vTimeTwoOn	= max(vTimeRedOn, vTimeGrnOn);
vTimeTwoOff = min(vTimeRedOff, vTimeGrnOff);
vTwoValid = (vTimeTwoOn < vTimeTwoOff) & (~isnan(vTimeRedOn));

nCells = numel(vTimeToDiv);
for iCell = 1:nCells
	vY = [0, 0, 1, 1] + (iCell - 1);
	
	vX = flipCoeff*[0, vTimeToDiv(iCell), vTimeToDiv(iCell), 0];
	hold on;
	fill(vX, vY, 'k', 'EdgeColor', 'k'); hold on;
	
	vX = flipCoeff* ...
		[vTimeRedOn(iCell), vTimeRedOff(iCell), ...
		vTimeRedOff(iCell), vTimeRedOn(iCell)];
	hold on;
	fill(vX, vY, 'r', 'EdgeColor', 'r'); hold on;
	
	vX = flipCoeff* ...
		[vTimeGrnOn(iCell), vTimeGrnOff(iCell), ...
		vTimeGrnOff(iCell), vTimeGrnOn(iCell)];
	hold on;
	fill(vX, vY, 'g', 'EdgeColor', 'g'); hold on;
	
	if (vTwoValid(iCell))
		vX = flipCoeff* ...
			[vTimeTwoOn(iCell), vTimeTwoOff(iCell), ...
			vTimeTwoOff(iCell), vTimeTwoOn(iCell)];
		hold on;
		fill(vX, vY, par.clr.darkGrn, 'EdgeColor', par.clr.darkGrn);
	end
end
end

function plotHistograms(vTimeToDiv, mCycleMeas, par, flipCoeff)
idxCycle = initTableCycleMeas();
vTimeRedOn	= mCycleMeas(:,idxCycle.timeRedOn);
vTimeRedOff = mCycleMeas(:,idxCycle.timeRedOff);
vTimeGrnOn	= mCycleMeas(:,idxCycle.timeGrnOn);
vTimeGrnOff = mCycleMeas(:,idxCycle.timeGrnOff);

uTimeDiv = mean(vTimeToDiv);
uRedOn = nanmean(vTimeRedOn);
uRedOff = nanmean(vTimeRedOff);
uGrnOn = nanmean(vTimeGrnOn);
uGrnOff = nanmean(vTimeGrnOff);

par.dyBar = 0;

% histograms
% par.nBins = 10;
par.lw = 1.4;
clrLightRed = [1, 0.3, 0.3];
clrDarkRed = [0.8, 0, 0];
clrLightGrn = [0.3, 1, 0.3];
clrDarkGrn = [0, 0.8, 0];

%{
par.nBins = sshist(vTimeRedOff);
[vCounts, vBarPos] = hist(vTimeRedOff, par.nBins);
dPos = vBarPos(2) - vBarPos(1);
vBarPos = [vBarPos(1) - dPos, vBarPos, vBarPos(end) + dPos];
vCounts = [0, vCounts, 0] ./ sum(~isnan(vTimeRedOff));
hold on;
stairs(flipCoeff*vBarPos, par.dyBar + vCounts, ...
	'color', 'r', 'linewidth', par.lw, 'linestyle', ':');
yMax = max([yMax, vCounts]);

par.nBins = sshist(vTimeRedOn);
[vCounts, vBarPos] = hist(vTimeRedOn, par.nBins);
dPos = vBarPos(2) - vBarPos(1);
vBarPos = [vBarPos(1) - dPos, vBarPos, vBarPos(end) + dPos];
vCounts = [0, vCounts, 0] ./ sum(~isnan(vTimeRedOn));
hold on;
stairs(flipCoeff*vBarPos, par.dyBar + vCounts, ...
	'color', 'r', 'linewidth', par.lw);
yMax = max([yMax, vCounts]);

par.nBins = sshist(vTimeGrnOff);
[vCounts, vBarPos] = hist(vTimeGrnOff, par.nBins);
dPos = vBarPos(2) - vBarPos(1);
vBarPos = [vBarPos(1) - dPos, vBarPos, vBarPos(end) + dPos];
vCounts = [0, vCounts, 0] ./ sum(~isnan(vTimeGrnOn));
hold on;
stairs(flipCoeff*vBarPos, -vCounts, ...
	'color', 'g', 'linewidth', par.lw, 'linestyle', ':');
yMax = max([yMax, vCounts]);

par.nBins = sshist(vTimeGrnOn);
[vCounts, vBarPos] = hist(vTimeGrnOn, par.nBins);
dPos = vBarPos(2) - vBarPos(1);
vBarPos = [vBarPos(1) - dPos, vBarPos, vBarPos(end) + dPos];
vCounts = [0, vCounts, 0] ./ sum(~isnan(vTimeGrnOn));
hold on;
stairs(flipCoeff*vBarPos, -vCounts, ...
	'color', 'g', 'linewidth', par.lw);
yMax = max([yMax, vCounts]);
%}

binStep = max(vTimeToDiv)/par.nBins;
par.vEdges = 0:binStep:max(vTimeToDiv);
par.vEdges2 = par.vEdges + binStep/2;
par.vEdges2 = [par.vEdges(1), par.vEdges2(1:(end - 1)), par.vEdges(end)];

vCounts = histc(vTimeRedOn, par.vEdges);
vCounts = [0; vCounts(1:(end - 1)); 0] ./ sum(~isnan(vTimeRedOn));
hold on;
% stairs(flipCoeff*par.vEdges, par.dyBar + vCounts, ...
% 	'color', 'r', 'linewidth', par.lw);
fill(flipCoeff*par.vEdges2, par.dyBar + vCounts, ...
	clrLightRed, 'edgecolor', clrLightRed);

vCounts = histc(vTimeRedOff, par.vEdges);
vCounts = [0; vCounts(1:(end - 1)); 0] ./ sum(~isnan(vTimeRedOff));
hold on;
% stairs(flipCoeff*par.vEdges, par.dyBar + vCounts, ...
% 	'color', 'r', 'linewidth', par.lw);
line(flipCoeff*par.vEdges2, par.dyBar + vCounts, ...
	'color', clrDarkRed, 'linewidth', par.lw);

vCounts = histc(vTimeGrnOn, par.vEdges);
vCounts = [0; vCounts(1:(end - 1)); 0] ./ sum(~isnan(vTimeGrnOn));
hold on;
% stairs(flipCoeff*par.vEdges, -vCounts, ...
% 	'color', 'g', 'linewidth', par.lw);
fill(flipCoeff*par.vEdges2, - vCounts, ...
	clrLightGrn, 'edgecolor', clrLightGrn, 'linewidth', par.lw);

vCounts = histc(vTimeGrnOff, par.vEdges);
% vCounts = vCounts ./ sum(~isnan(vTimeGrnOn));
vCounts = [0; vCounts(1:(end - 1)); 0] ./ sum(~isnan(vTimeGrnOff));
hold on;
% stairs(flipCoeff*par.vEdges, -vCounts, ...
% 	'color', 'g', 'linewidth', par.lw);
line(flipCoeff*par.vEdges2, - vCounts, ...
	'color', clrDarkGrn, 'linewidth', par.lw);

hold on;
line(flipCoeff*[0, max(vTimeToDiv)], [0, 0], ...
	'color', 'k', 'linewidth', par.lw)
%{
% horizontal bar
par.lw = 3.5;

vX = flipCoeff*[0, uTimeDiv, uTimeDiv, 0];
vY = [0, 0, par.dyBar, par.dyBar];
hold on;
fill(vX, vY, 'k', 'edgecolor', 'k');

hold on;
line(flipCoeff*[uRedOn, uRedOn], [0, par.dyBar], ...
	'linewidth', par.lw, 'color', 'r');
hold on; line(flipCoeff*[uRedOff, uRedOff], [0, par.dyBar], ...
	'color', 'r', 'linewidth', 2);
hold on; line(flipCoeff*[uRedOff, uRedOff], [0, par.dyBar], ...
	'color', 'b', 'linewidth', 2, 'linestyle', ':');
hold on;
line(flipCoeff*[uGrnOn, uGrnOn], [0, par.dyBar], ...
	'linewidth', par.lw, 'color', 'g');
hold on;
line(flipCoeff*[uGrnOff, uGrnOff], [0, par.dyBar], ...
	'color', 'g', 'linewidth', par.lw);
hold on;
line(flipCoeff*[uGrnOff, uGrnOff], [0, par.dyBar], ...
	'color', 'b', 'linewidth', par.lw, 'linestyle', ':');
	%}
	% hold on;
	% line(flipCoeff*[uTimeDiv, uTimeDiv], [0, par.dyBar], ...
	% 	'linewidth', par.lw, 'color', 'c');
	
end