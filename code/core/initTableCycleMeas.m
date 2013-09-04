%{
% Define the structure of the table with cell cycle measurements
%
% Andrey Kan
% akan@wehi.edu.au
%
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function idx = initTableCycleMeas()

idx.timeRedMax	= 0;
idx.levRedMax	= 0;
idx.timeRedOn	= 0;
idx.timeRedOff	= 0;

idx.timeGrnMax	= 0;
idx.levGrnMax	= 0;
idx.timeGrnOn	= 0;
idx.timeGrnOff	= 0;

fields = fieldnames(idx);
for i = 1:numel(fields)
	idx.(fields{i}) = i;
end
idx.last = numel(fields);
end