%{
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function rss = evaluateModel(modelFunc, vPar, ctx, mObserved, vDispTimes)
vPredicted = 100.*modelFunc(ctx, vPar);
mPredicted = repmat(vPredicted, 1, 3);
mResiduals = (mPredicted - mObserved).^2;
%mResiduals = log((mPredicted - mObserved).^2);
rss = sum(mResiduals(:));

if (nargin >= 5)
	fprintf('fitted parameters:\n');
	ctx.printFcn(vPar);
	fprintf('RSS = %.3f\n', rss);
	
	mTimes = repmat(ctx.vTimes, 1, 3);
	ctx.vTimes = vDispTimes;
	
	figure;
	plot(mTimes(:), mObserved(:), 'marker', '*', 'linestyle', 'none');
	hold on;
	plot(vDispTimes, 100.*modelFunc(ctx, vPar), 'color', 'g');
	xlim([0, vDispTimes(end)]);
 	ylim([0, 12]);
end
end