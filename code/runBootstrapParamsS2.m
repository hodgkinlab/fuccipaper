%{
% bootstrapped CI for stretched and constant G2M models
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
function runBootstrapParamsS2()
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

nSamples = 1000;
parStretCpG = NaN(1, nSamples);
parConstCpG = NaN(1, nSamples);
parStretCD8 = NaN(1, nSamples);
parConstCD8 = NaN(1, nSamples);

n = numel(mCpG);
rowIdx = repmat(1:size(mCpG, 1), size(mCpG, 2), 1)';

matlabpool open;
parfor i = 1:nSamples
	randIndex = randsample(3, n, true);
	colIdx = reshape(randIndex, size(mCpG));
	currCpG = mCpG(sub2ind(size(mCpG), rowIdx, colIdx));
	currCD8 = mCD8(sub2ind(size(mCD8), rowIdx, colIdx));
	parSample = fitOneSample(ctx, currCpG, currCD8);
	parStretCpG(i) = parSample(1);
	parConstCpG(i) = parSample(2);
	parStretCD8(i) = parSample(3);
	parConstCD8(i) = parSample(4);
end
matlabpool close;

prctile(parStretCpG, [2.5, 97.5])
prctile(parConstCpG, [2.5, 97.5])./60
prctile(parStretCD8, [2.5, 97.5])
prctile(parConstCD8, [2.5, 97.5])./60

end

function parSample = fitOneSample(ctx, currCpG, currCD8)

% CpG B cells, MLE lognormal uLog = 12.34; sLog = 3.48 hours
ctx.printFcn = @printStretchedG2M;
ctx.listPar = {(12.34)*60, (3.48)*60, 0.01:0.001:1};
par = fitDataSample(ctx, currCpG, @stretchedG2M);
parSample(1) = par(3);

ctx.printFcn = @printConstantS;
ctx.listPar = {(12.34)*60, (3.48)*60, 0:0.5:30*60, 0.73};
par = fitDataSample(ctx, currCpG, @constantS);
parSample(2) = par(3);

fprintf('---------------\n');

% CD8+ T cells, MLE lognormal uLog = 12.95; sLog = 3.46 hours
ctx.printFcn = @printStretchedG2M;
ctx.listPar = {(12.95)*60, (3.46)*60, 0.01:0.001:1};
par = fitDataSample(ctx, currCD8, @stretchedG2M);
parSample(3) = par(3);

ctx.printFcn = @printConstantS;
ctx.listPar = {(12.95)*60, (3.46)*60, 0:0.5:30*60, 0.65};
par = fitDataSample(ctx, currCD8, @constantS);
parSample(4) = par(3);

end

function vPar = fitDataSample(ctx, mData, modelFunc)
% optimization using fmincon
vPar0	= cellfun(@mean, ctx.listPar);
vParLB	= cellfun(@min, ctx.listPar);
vParUB	= cellfun(@max, ctx.listPar);

opts = optimset( ...
	'Display', 'final', 'Algorithm', ctx.sOptimAlg, ...
	'Maxiter', ctx.nMaxIters, 'MaxFunEvals', ctx.nFunEvals, ...
	'TolFun', ctx.tolerance, 'TolX', ctx.tolerance);

objfun = @(vPar) evaluateModel(modelFunc, vPar, ctx, mData);
vPar = fmincon( ...
	objfun, vPar0, [], [], [], [], vParLB, vParUB, [], opts);
end