%{
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function vPar = fitDataBrdU(ctx, mData, modelFunc, bExhaustive)
if (bExhaustive)
	% exhaustive search optimization
	nPar = numel(ctx.listPar);
	vSizes = cellfun(@numel, ctx.listPar);
	nCmb = prod(vSizes);
	mPar = NaN(nCmb, nPar);
	for i = 1:nPar
		x = ctx.listPar{i}(:);
		s = vSizes; s(i) = [];
		x = reshape(x(:, ones(1, prod(s))), [length(x) s]);
		m = permute(x, [2:i 1 i+1:nPar]);
		mPar(:,i) = m(:);
	end
	vRSS = NaN(1, nCmb);
	fprintf('combinations total: %d; current:\n%05d', nCmb, 1);
	for i = 1:nCmb
		if (~bitand(i, 127))
			fprintf('\b\b\b\b\b');
			fprintf('%05d', i);
		end
		vCurPar = mPar(i,:);
		vRSS(i)	= evaluateModel(modelFunc, vCurPar, ctx, mData);
	end
	fprintf('\n');
	[~, i] = min(vRSS);
	vPar = mPar(i,:);
else
	% optimization using fmincon
	vPar0	= cellfun(@mean, ctx.listPar);
	vParLB	= cellfun(@min, ctx.listPar);
	vParUB	= cellfun(@max, ctx.listPar);
	
	opts = optimset( ...
		'Display', 'final', 'Algorithm', ctx.sOptimAlg, ...
		'Maxiter', ctx.nMaxIters, 'MaxFunEvals', ctx.nFunEvals, ...
		'TolFun', ctx.tolerance, 'TolX', ctx.tolerance);
	
	objfun = @(vPar) evaluateModel(modelFunc, vPar, ctx, mData);
	vPar = fmincon( ...
		objfun, vPar0, [], [], [], [], vParLB, vParUB, [], opts);
end

evaluateModel(modelFunc, vPar, ctx, mData, ctx.vDispTimes);
end