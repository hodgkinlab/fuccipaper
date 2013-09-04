%{
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function printConstantS(vPar)
fprintf('mean = %0.2f\n', vPar(1)/60);
fprintf('sdev = %0.2f\n', vPar(2)/60);
fprintf('S-phase = %0.2f\n', vPar(3)/60);
fprintf('r(SG2M) = %d%%\n', round(100*vPar(4)));
end