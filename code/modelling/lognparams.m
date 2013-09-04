%{
% input: lognormal mean and standard deviation
% output: parameters of the associated normal distribution
%
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function [m,s] = lognparams(mu,sigma)

m = log(mu.^2./sqrt(sigma.^2+mu.^2));
s = sqrt(log((sigma./mu).^2+1));