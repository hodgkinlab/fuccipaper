%{
Histogram versions of alpha-plots
%
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function plotAlphaHist(vAxesAHist, mCellData, mCycleMeas, models, par)
if (isempty(vAxesAHist))
	return
end

idxData = initTableCellData();
idxCycle = initTableCycleMeas();
vc = mCellData(:,idxData.stop) - mCellData(:,idxData.start);
va = mCycleMeas(:,idxCycle.timeGrnOn);
vb = vc - mCycleMeas(:,idxCycle.timeGrnOn);

v = vAxesAHist{par.iDataset};
h1 = v(1);
h2 = v(2);
h3 = v(3);

iModLExp = 3;
iModExpG = 4;
iModLogn = 5;

xmax = 30;
vt = 0.025:0.25:xmax;
mu = models{iModLogn}.vParam(1); % stretched lognormal
sigma = models{iModLogn}.vParam(2);
r1 = models{iModLogn}.vParam(3);
r2 = models{iModLogn}.vParam(4);
vpdf = lognpdf(vt, mu, sigma);
vpdf = vpdf ./ max(vpdf);

t0 = models{iModLExp}.vParam(1);
lambda = models{iModLExp}.vParam(2);
vLExp = lambda.*exp(-lambda.*(vt - t0));
vLExp(vt < t0) = 0;
vLExp = vLExp ./ max(vLExp);

uHat = models{iModExpG}.vParam(1);
sHat = models{iModExpG}.vParam(2);
lHat = models{iModExpG}.vParam(3);

lw = 1;

axes(h1);
[n, xout] = hist(va, 1:xmax);
bar(xout, n./max(n));
hold on;
plot(vt.*r1, vpdf, 'g--', 'LineWidth', lw);

vy = lambda.*exp(-lambda.*(vt));
hold on;
plot(vt, vy ./ max(vy), 'r', 'LineWidth', lw);

vy = lHat.*exp(-lHat.*(vt));
hold on;
plot(vt, vy ./ max(vy), 'c', 'LineWidth', lw);
xlim([0, xmax]);
set(gca, 'YTick', []);
text(20, 0.85, sprintf('m = %0.1f h', nanmean(va)));
text(20, 0.7, sprintf('s = %0.1f h', nanstd(va)));

axes(h2);
[n, xout] = hist(vb, 1:xmax);
bar(xout, n./max(n));
hold on;
plot(vt.*r2, vpdf, 'g--', 'LineWidth', lw);
hold on;
line([0 xmax], [0 0], 'color', 'r', 'LineWidth', lw);
hold on;
line([t0 t0], [0, 1], 'color', 'r', 'LineWidth', lw);

vy = normpdf(vt, uHat, sHat);
hold on;
plot(vt, vy ./ max(vy), 'c', 'LineWidth', lw);

xlim([0, xmax]);
set(gca, 'YTick', []);
text(20, 0.85, sprintf('m = %0.1f h', nanmean(vb)));
text(20, 0.7, sprintf('s = %0.1f h', nanstd(vb)));

axes(h3);
[n, xout] = hist(vc, 1:xmax);
bar(xout, n./max(n));
hold on;
plot(vt, vpdf, 'g--', 'LineWidth', lw);
hold on;
plot(vt, vLExp, 'r', 'LineWidth', lw);

% EMG
vy = 0.5*lHat ...
	.*exp(0.5.*lHat.*(2*uHat + lHat*sHat*sHat - 2.*vt)) ...
	.*erfc((uHat + lHat*sHat*sHat - vt)./sqrt(2)./sHat);
hold on;
plot(vt, vy ./ max(vy), 'c', 'LineWidth', lw);

xlim([0, xmax]);
set(gca, 'YTick', []);
text(20, 0.85, sprintf('m = %0.1f h', nanmean(vc)));
text(20, 0.7, sprintf('s = %0.1f h', nanstd(vc)));