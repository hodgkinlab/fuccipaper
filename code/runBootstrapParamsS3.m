%{
% Bootstrapped confidence intervals for supplementary table S3
%
% Data is drawn from experiments:
% 20111118 - B cells, CpG
% 20120413 - B cells, aCD40
% 20120309 - CD8 T cells
% 20120316 - OT1 cells
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
function runBootstrapParamsS3()
clc, clearvars, close all, addpath(genpath('.'));
reset(RandStream.getGlobalStream, 2013);
par.sDataPath = '..\data';
par.sExpNames = {'20111118', '20120413', '20120309', '20120316'};

par.nSamples = 1000;% used for computing bootstrap confidence intervals
par.nParams = 6;

groupDataFUCCI(par);
par.vPropSG2M = [0.73 0.78 0.65 0.72];
par.vPropG1 = 1 - par.vPropSG2M;

idxData		= initTableCellData();
idxCycle	= initTableCycleMeas();
par.idxData = idxData;
par.idxCycle = idxCycle;

paramVal = NaN(numel(par.sExpNames), par.nParams);
boundsLo = NaN(numel(par.sExpNames), par.nParams);
boundsHi = NaN(numel(par.sExpNames), par.nParams);

iDataset = 1;
matlabpool open;
for iExp = 1:numel(par.sExpNames)
	fprintf('%d\n', iExp);
	
	% load current experiment
	if (1)
		par.sExpName = par.sExpNames{iExp};
		par.sInfoFile = fullfile( ...
			par.sDataPath, [par.sExpName, '-info.mat']);
		par.sDataFile = fullfile( ...
			par.sDataPath, [par.sExpName, '-cycle-meas.mat']);
		% recall that times are measured in hours
		load(par.sInfoFile, 'nMinPerFrame', 'vGroups', 'groupNames');
		load(par.sDataFile, 'mCellData', 'mCycleMeas');
	end
	% filter positions of interest (see groupDataFUCCI)
	if (1)
		mCellData = mCellData(vGroups > 0,:);
		mCycleMeas = mCycleMeas(vGroups > 0,:);
		vGroups = vGroups(vGroups > 0);
		if (isempty(vGroups))
			error('empty selection');
		end
	end
	
	mCellDataAll = mCellData;
	mCycleMeasAll = mCycleMeas;
	nGroups = max(vGroups);
	if (nGroups ~= 1)
		error('in our case we should have exactly one group');
	end
	for iGroup = 1:nGroups
		par.iDataset = iDataset;
		% select a group of cells
		if (1)
			vSelect = (vGroups == iGroup);
			nCells = sum(vSelect);
			if (nCells == 0)
				error('empty selection');
			else
				%fprintf('(%d cells) ===\n', nCells);
			end
			mCellData = mCellDataAll(vSelect,:);
			mCycleMeas = mCycleMeasAll(vSelect,:);
		end
		
		par.mCellData = mCellData;
		par.mCycleMeas = mCycleMeas;
		
		currVal = estimateParams(par, false);
		paramVal(iExp,:) = round(100*currVal)/100;
		
		mParams = NaN(par.nSamples, par.nParams);
		parfor i = 1:par.nSamples
			mParams(i,:) = estimateParams(par, true);
		end
		
		mLims = prctile(mParams, [2.5, 97.5]);
		mLims = round(100*mLims)/100;
		boundsLo(iExp,:) = mLims(1,:);
		boundsHi(iExp,:) = mLims(2,:);
		
		a = mParams(:,5);
		fprintf('exp-gauss failed in %d%% of cases\n', ...
			round(100*sum(isnan(a))/numel(a)));
		
		
		iDataset = iDataset + 1;
	end
end
matlabpool close;
disp(num2str(boundsLo));
disp(' ');
disp(num2str(paramVal));
disp(' ');
disp(num2str(boundsHi));
disp(' ');
disp(num2str(boundsHi - paramVal >= 0));
disp(' ');
disp(num2str(boundsLo - paramVal <= 0));
disp(' ');
disp(num2str(boundsHi - boundsLo >= 0));

save(fullfile(par.sDataPath, 'boundsS3.mat'), ...
	'paramVal', 'boundsLo', 'boundsHi');
end

function mParamsSample = estimateParams(par, doSampling)
mParamsSample = NaN(1, par.nParams);

idxData = par.idxData;
idxCycle = par.idxCycle;
mCellData = par.mCellData;
mCycleMeas = par.mCycleMeas;

T = mCellData(:,idxData.stop) - mCellData(:,idxData.start);
b2 = T - mCycleMeas(:,idxCycle.timeGrnOn);
ab1 = mCycleMeas(:,idxCycle.timeGrnOn);

vSelect = ~isnan(b2);
T = T(vSelect);
b2 = b2(vSelect);
ab1 = ab1(vSelect);

if (doSampling)
	randMap = randsample(1:numel(T), numel(T), true);
	T = T(randMap);
	b2 = b2(randMap);
	ab1 = ab1(randMap);
end

s = std(T);
sZ = std(b2);
mParamsSample(2) = sZ;
mParamsSample(3) = s;
mParamsSample(4) = sZ/s;

if (s < sZ)
	uaUpBnd = NaN;
else
	uaUpBnd = sqrt(s*s - sZ*sZ);
end
mParamsSample(5) = uaUpBnd;

% m = mean(T);
y = skewness(T);
% uHat = m - s*nthroot(y/2, 3);
% sHat = s*s*(1 - nthroot(y/2, 3).^2);
kHat = s*nthroot(y/2, 3);
mParamsSample(6) = kHat;

cv = cov(ab1(:), b2(:));
mParamsSample(1) = cv(1,2);
end