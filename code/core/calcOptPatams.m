%{
% Andrey Kan, et al.
% akan@wehi.edu.au
%
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function vParamOpt = calcOptPatams(vc, iExp)
if (1)
% pre-computed values
m = [
7.30000, 0.17000 
8.40000, 0.28500
8.50000, 0.21500
7.50000, 0.30500];
vParamOpt = m(iExp,:);
return
end

vParam1 = 1:0.1:30;
vParam2 = 0:0.005:30;
n1 = numel(vParam1);
n2 = numel(vParam2);
m = NaN(n1, n2);
matlabpool open
parfor i1 = 1:n1
	for i2 = 1:n2
		vParam = [vParam1(i1), vParam2(i2)];
		[~, m(i1,i2)] = kstest(vc, ...
			[vc, 1 - exp(-vParam(2)*(vc - vParam(1)))]);
	end
end
matlabpool close
[~, i] = max(m(:));
[i1, i2] = ind2sub(size(m), i(1));

vParamOpt = [vParam1(i1), vParam2(i2)];
fprintf('%0.5f, %0.5f\n', vParam1(i1), vParam2(i2));

end