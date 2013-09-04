%{
% Fit stretched G2M model to BrdU experiment results
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
clc, clearvars, close all; addpath(genpath('.'));
reset(RandStream.getGlobalStream, 2013);

% load data
sDataPath = fullfile('..', 'data');
if (exist(fullfile(sDataPath, 'BrdU-data.mat'), 'file') == 2)
	load(fullfile(sDataPath, 'BrdU-data.mat'), 'mCpG', 'mCD8', 'vTimes');
else
	sDataFile = fullfile(sDataPath, 'SH1.93 G2M summary.xlsx');
	mCpG = xlsread(sDataFile, 4, 'E2:G7');
	mCD8 = xlsread(sDataFile, 4, 'Q2:S7');
	vTimes = xlsread(sDataFile, 4, 'A2:A7');
	save(fullfile(sDataPath, 'BrdU-data.mat'), 'mCpG', 'mCD8', 'vTimes');
end

% define parameters
if (1)
	% times are measured in minutes
	ctx.vTimes		= vTimes;
	ctx.vDispTimes	= 1:5:275;
	ctx.numIntStep	= 2;
	ctx.limSdevMax	= 6;
	ctx.nMaxIters	= Inf;
	ctx.nFunEvals	= Inf;
	ctx.tolerance	= 1.0e-9;
	ctx.sOptimAlg	= 'interior-point';
end

% CpG B cells, MLE lognormal uLog = 12.34; sLog = 3.48 hours
ctx.printFcn = @printStretchedG2M;
ctx.listPar = {(12.34)*60, (3.48)*60, 0.01:0.001:1};
fitDataBrdU(ctx, mCpG, @stretchedG2M, 1);

ctx.printFcn = @printConstantS;
ctx.listPar = {(12.34)*60, (3.48)*60, 0:0.5:30*60, 0.73};
fitDataBrdU(ctx, mCpG, @constantS, 1);

fprintf('---------------\n');

% CD8+ T cells, MLE lognormal uLog = 12.95; sLog = 3.46 hours
ctx.printFcn = @printStretchedG2M;
ctx.listPar = {(12.95)*60, (3.46)*60, 0.01:0.001:1};
fitDataBrdU(ctx, mCD8, @stretchedG2M, 1);

ctx.printFcn = @printConstantS;
ctx.listPar = {(12.95)*60, (3.46)*60, 0:0.5:30*60, 0.65};
fitDataBrdU(ctx, mCD8, @constantS, 1);