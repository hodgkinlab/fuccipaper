%{
% Simple stretched G2M model is based on
% the steady state model for a growing culture
% (see Dowling et al., J. Royal Society, 2005)
%
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function y = constantG2M(ctx, vPar)
u = vPar(1); % mean of the source lognormal distribution
s = vPar(2); % st.dev. of the source lognormal distribution
a = vPar(3); % t_G2M = r * t_total

y = NaN(size(ctx.vTimes));
for i = 1:numel(ctx.vTimes)
	y(i) = calcConstProp(u, s, a, ...
		ctx.numIntStep, (u + ctx.limSdevMax*s), ctx.vTimes(i));
end
end

function f = calcConstProp(mu, sigma, tG2M, dt, tEnd, pulseLen)
if (pulseLen >= tG2M)
	f = 0;
	return
end

% PDF over the space of total and remaining times
[m, s] = lognparams(mu, sigma);
k = fzero('normalise', 0, [], m, s, 0:dt:tEnd);
fun = @(tot, rem) 2.*k.*exp(-k.*(tot - rem)).*lognpdf(tot, m, s);

% integrate over the triangular region bounded by
% 0 <= tot <= (tG2M - pulseLen)  and
% 0 <= rem <= (tot)
maxRemTime = @(tot) (tot);
f1 = integral2(fun, 0, (tG2M - pulseLen), 0, maxRemTime);

% integrate over the region bounded by
% (tG2M - pulseLen) <= tot <= tEnd  and
% 0 <= rem <= (tG2M - pulseLen)
f2 = integral2(fun, (tG2M - pulseLen), tEnd, 0, (tG2M - pulseLen));

f = min(f1 + f2, 1);
end