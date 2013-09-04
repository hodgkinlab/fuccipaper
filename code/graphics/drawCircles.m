%{
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function drawCircles(prop, par)
step = 0.001;
width = 1;

% circles
vX = par.circle.r1*sin(0:step:2*pi) + par.circle.x0;
vY = par.circle.r1*cos(0:step:2*pi) + par.circle.y0;
hold on; plot(vX, vY, 'color', 'k', 'linewidth', width);
vX = par.circle.r2*sin(0:step:2*pi) + par.circle.x0;
vY = par.circle.r2*cos(0:step:2*pi) + par.circle.y0;
hold on; plot(vX, vY, 'color', 'k', 'linewidth', width);

% arrow
%{
hold on;
line([par.circle.x0, par.circle.x0], ...
	[par.circle.r2, par.circle.r2 + 1] + par.circle.y0, ...
	'color', 'k', 'linewidth', width);
vX = (par.circle.r2 + 1)*sin(0.17:step:0.5) + par.circle.x0;
vY = (par.circle.r2 + 1)*cos(0.17:step:0.5) + par.circle.y0;
hold on; plot(vX, vY, 'color', 'k', 'linewidth', width);
hold on;
line([vX(end), vX(end) - 0.5], ...
	[vY(end), vY(end) + 0.7], ...
	'color', 'k', 'linewidth', width);
%}

angleOn = (1 - prop)*2*pi;
vAng2 = angleOn:step:2*pi;
vAng1 = fliplr(vAng2);
vX = [par.circle.r2*sin(vAng2), par.circle.r1*sin(vAng1)];
vY = [par.circle.r2*cos(vAng2), par.circle.r1*cos(vAng1)];
vX = vX + par.circle.x0;
vY = vY + par.circle.y0;

hold on;
fill(vX, vY, par.clr.lightGrn, ...
	'edgecolor', par.clr.darkGrn, 'linewidth', width);
end
