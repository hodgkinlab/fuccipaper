%{
% Defines the structure of the table that lists cells
%
% Andrey Kan, et al.
% akan@wehi.edu.au
%
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function idx = initTableCellData()
idx.cell	= 0;
idx.start	= 0;
idx.stop	= 0;
idx.dvsn	= 0;
idx.fate	= 0; % 0 - died; 1 - divided; 2 - lost
idx.pos		= 0;
idx.well	= 0;
idx.mom		= 0;
idx.sis		= 0;
idx.child1	= 0;
idx.child2	= 0;

fields = fieldnames(idx);
for i = 1:numel(fields)
	idx.(fields{i}) = i;
end
idx.last = numel(fields);
end