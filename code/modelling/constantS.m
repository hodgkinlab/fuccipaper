%{
% Simple stretched SG2M + constant S model is based on
% the steady state model for a growing culture
% (see Dowling et al., J. Royal Society, 2005)
%
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function y = constantS(ctx, vPar)
u = vPar(1); % mean of the source lognormal distribution
s = vPar(2); % st.dev. of the source lognormal distribution
c = vPar(3); % duration of the S phase
r = vPar(4); % T_SG2M = r*T_div

y = NaN(size(ctx.vTimes));
for i = 1:numel(ctx.vTimes)
	y(i) = calcPropG2M(u, s, c, r, ...
		ctx.numIntStep, (u + ctx.limSdevMax*s), ctx.vTimes(i));
end
end

%{
% mu, sigma are mean and standard deviation of division time;
% dt is the time increment, dt=0.01 for good results;
% tEnd s the maximum time in the simulation, want tEnd>mu+4*sigma;
% (internal functions were implemented by Mark Dowling)
%}
function f = calcPropG2M(mu, sigma, c, r, dt, tEnd, pulseLen)
% PDF over the space of total and remaining times
[m, s] = lognparams(mu, sigma);
k = fzero('normalise', 0, [], m, s, 0:dt:tEnd);
fun = @(tot, rem) 2.*k.*exp(-k.*(tot - rem)).*lognpdf(tot, m, s);

% integrate over the triangular region bounded by
% ((pulseLen + c)/r) <= tot <= tEnd and
% 0 <= rem <= (r*tot - pulseLen - c)
maxRemTime = @(tot) (r*tot - pulseLen - c);
f = integral2(fun, (pulseLen + c)/r, tEnd, 0, maxRemTime);
f = min(f, 1);
end