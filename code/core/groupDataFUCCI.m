%{
% Group data for processing for the FUCCI paper. This is needed,
% because an experiment can involve different types of cells
% or conditions in different chambers.
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
function groupDataFUCCI(par)
idxData = initTableCellData();

sExpName = '20111118'; nMinPerFrame = 3;
par.sInfoFile = fullfile(par.sDataPath, [sExpName, '-info.mat']);
par.sDataFile = fullfile(par.sDataPath, [sExpName, '-cycle-meas.mat']);
load(par.sDataFile, 'mCellData');

vGroups = ones(size(mCellData, 1), 1);
groupNames = {'B cells, CpG'};
save(par.sInfoFile, 'nMinPerFrame', 'vGroups', 'groupNames');
%{
par.sExpName = '20111118-Jie';
par.sInfoFile = fullfile(par.sDataPath, [par.sExpName, '-info.mat']);
par.sDataFile = fullfile(par.sDataPath, [par.sExpName, '-cycle-meas.mat']);

load(par.sInfoFile, 'nMinPerFrame');
load(par.sDataFile, 'mCellData');
vDvsn = mCellData(:,idxData.dvsn); %#ok<*SUSENS>

vGroups = zeros(size(mCellData, 1), 1);
vGroups(vDvsn == 1) = 1;
%vGroups(vDvsn == 2) = 2;
%vGroups(vDvsn == 3) = 3;
% groupNames = {'B-cells, CpG gen.1', 'B-cells, CpG gen.2', 'B-cells, CpG gen.3'};
vGroups(vDvsn == 2) = 1;
vGroups(vDvsn == 3) = 1;
groupNames = {'B-cells, CpG (Jie)'};
save(par.sInfoFile, 'nMinPerFrame', 'vGroups', 'groupNames');
%}
%{
par.sExpName = '20111118-NICTA';
par.sInfoFile = fullfile(par.sDataPath, [par.sExpName, '-info.mat']);
par.sDataFile = fullfile(par.sDataPath, [par.sExpName, '-cycle-meas.mat']);

load(par.sInfoFile, 'nMinPerFrame');
load(par.sDataFile, 'mCellData');

vGroups = ones(size(mCellData, 1), 1); %#ok<*NASGU>
groupNames = {'B-cells, CpG (NICTA)'};
save(par.sInfoFile, 'nMinPerFrame', 'vGroups', 'groupNames');
%}


par.sExpName = '20120309'; nMinPerFrame = 4;
par.sInfoFile = fullfile(par.sDataPath, [par.sExpName, '-info.mat']);
par.sDataFile = fullfile(par.sDataPath, [par.sExpName, '-cycle-meas.mat']);
load(par.sDataFile, 'mCellData');

vPos = mCellData(:,idxData.pos);
vGroups = zeros(size(mCellData, 1), 1);
vGroups((vPos >= 49) & (vPos <= 72)) = 1;
groupNames = {'CD8+ T cells'}; % +IL2
save(par.sInfoFile, 'nMinPerFrame', 'vGroups', 'groupNames');


par.sExpName = '20120316'; nMinPerFrame = 10;
par.sInfoFile = fullfile(par.sDataPath, [par.sExpName, '-info.mat']);
par.sDataFile = fullfile(par.sDataPath, [par.sExpName, '-cycle-meas.mat']);
load(par.sDataFile, 'mCellData');

vPos = mCellData(:,idxData.pos);
vGroups = zeros(size(mCellData, 1), 1);
vGroups((vPos >= 1) & (vPos <= 30)) = 1;
groupNames = {'OT-I T cells'}; % high affinity N, +IL2
save(par.sInfoFile, 'nMinPerFrame', 'vGroups', 'groupNames');


par.sExpName = '20120413'; nMinPerFrame = 10;
par.sInfoFile = fullfile(par.sDataPath, [par.sExpName, '-info.mat']);
par.sDataFile = fullfile(par.sDataPath, [par.sExpName, '-cycle-meas.mat']);
load(par.sDataFile, 'mCellData');

vPos = mCellData(:,idxData.pos);
vGroups = zeros(size(mCellData, 1), 1);
vGroups((vPos >= 1) & (vPos <= 36)) = 1;
groupNames = {'B cells, \alphaCD40'}; % 40 ug/mL aCD40
save(par.sInfoFile, 'nMinPerFrame', 'vGroups', 'groupNames');
end