%{
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function out = normalise(k,m,s,tau)

Y = lognequildivfun(tau,m,s,k);

out = trapz(tau,Y)-0.5;
