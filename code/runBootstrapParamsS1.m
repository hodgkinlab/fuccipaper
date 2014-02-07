%{
% Bootstrapped confidence intervals for supplementary table S1
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
function runBootstrapParamsS1()
clc, clearvars, close all, addpath(genpath('.'));
reset(RandStream.getGlobalStream, 2013);
par.sDataPath = '..\data';
par.sExpNames = {'20111118', '20120413', '20120309', '20120316'};

par.nSamples = 1000;% used for computing bootstrap confidence intervals
par.nParams = 14;

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
		% split siblings
		if (1)
			% shuffle stuff
			vPerm = randperm(nCells);
			mCellData = mCellData(vPerm,:);
			mCycleMeas = mCycleMeas(vPerm,:);
			
			% this code is for multiple cell wells;
			% for single cell wells, we can simply look at well IDs
			vMoms = mCellData(:,idxData.mom);
			vWells = mCellData(:,idxData.well);
			if (min(vMoms) < 1)
				error('mother cell index should start with 1');
			end
			maxMom = max(vMoms);
			vUniqueMoms = vWells .* maxMom + vMoms;
			
			% put sister cells together
			[vUniqueMoms, vOrder] = sort(vUniqueMoms);
			mCellData = mCellData(vOrder,:);
			mCycleMeas = mCycleMeas(vOrder,:);
			
			% mark pairs and singles
			vSingle = false(1, nCells);
			iCell = 1;
			while (iCell < nCells)
				if (vUniqueMoms(iCell) == vUniqueMoms(iCell + 1))
					iCell = iCell + 2;
				else
					vSingle(iCell) = true;
					iCell = iCell + 1;
				end
			end
			if (iCell == nCells)
				vSingle(end) = true;
			end
			% unpaired cells
			% 			mCellDataX	= mCellData(vSingle,:);
			% 			mCycleMeasX = mCycleMeas(vSingle,:);
			
			% split siblings
			mCellDataY	= mCellData(~vSingle,:);
			mCycleMeasY = mCycleMeas(~vSingle,:);
			mCellDataA	= mCellDataY(1:2:end,:);
			mCycleMeasA = mCycleMeasY(1:2:end,:);
			% 			mCellDataB	= mCellDataY(2:2:end,:);
			% 			mCycleMeasB = mCycleMeasY(2:2:end,:);
			
			% 			nUnpaired = size(mCellDataX, 1);
			% 			nPairs = size(mCellDataA, 1);
			%fprintf('#pairs = %d; #unpaired = %d; ratio = %0.2f\n\n', ...
			%	nPairs, nUnpaired, (double(nUnpaired)/double(2*nPairs)));
		end
		par.mCellDataA = mCellDataA;
		par.mCycleMeasA = mCycleMeasA;
		
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
		
		a = mParams(:,3);
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

for i = 1:4
	fprintf('(%0.2f;\n%0.2f)\n', boundsLo(1,2), boundsHi(1,2));
	fprintf('(%0.2f;\n%0.2f)\n\n', boundsLo(1,1), boundsHi(1,1));
	
	fprintf('(%0.2f;\n%0.2f)\n', boundsLo(1,3), boundsHi(1,3));
	fprintf('(%0.2f;\n%0.2f)\n', boundsLo(1,4), boundsHi(1,4));
	fprintf('(%0.2f;\n%0.2f)\n\n', boundsLo(1,5), boundsHi(1,5));
	
	fprintf('(%0.2f;\n%0.2f)\n', boundsLo(1,6), boundsHi(1,6));
	fprintf('(%0.2f;\n%0.2f)\n\n', boundsLo(1,7), boundsHi(1,7));
	
	fprintf('(%0.2f;\n%0.2f)\n', boundsLo(1,8), boundsHi(1,8));
	fprintf('(%0.2f;\n%0.2f)\n\n', boundsLo(1,9), boundsHi(1,9));
	
	fprintf('(%0.2f;\n%0.2f)\n', boundsLo(1,11), boundsHi(1,11));
	fprintf('(%0.2f;\n%0.2f)\n\n', boundsLo(1,10), boundsHi(1,10));
	
	fprintf('(%0.2f;\n%0.2f)\n', boundsLo(1,13), boundsHi(1,13));
	fprintf('(%0.2f;\n%0.2f)\n', boundsLo(1,14), boundsHi(1,14));
	fprintf('(%0.2f;\n%0.2f)\n\n', boundsLo(1,12), boundsHi(1,12));
	
	fprintf(' ---------------- \n\n');
end

save(fullfile(par.sDataPath, 'boundsS1.mat'), ...
	'paramVal', 'boundsLo', 'boundsHi');
end

function mParamsSample = estimateParams(par, doSampling)
mParamsSample = NaN(1, par.nParams);

idxData = par.idxData;
idxCycle = par.idxCycle;
mCellDataA = par.mCellDataA;
mCycleMeasA = par.mCycleMeasA;

vc = mCellDataA(:,idxData.stop) - mCellDataA(:,idxData.start);
va = mCycleMeasA(:,idxCycle.timeGrnOn);

if (doSampling)
	randMap = randsample(1:numel(vc), numel(vc), true);
	vc = vc(randMap);
	va = va(randMap);
end

vb = vc - va;
va = va(:);
vb = vb(:);

% 'lag-exponential';
[f, x] = ecdf(vc);
x = x(f < 1); f = f(f < 1); f = log(1 - f);
p = polyfit(x, f, 1);
vParam(1) = -p(2)/p(1); vParam(2) = -p(1);
mParamsSample(1) = vParam(1);
mParamsSample(2) = vParam(2);

% 'exp+gauss';
m = mean(vc);
s = std(vc);
y = skewness(vc);

if (y <= 0)
	uHat = NaN;
	sHat = NaN;
	lHat = NaN;
else
	uHat = m - s*nthroot(y/2, 3);
	sHat = s*sqrt(1 - nthroot(y/2, 3).^2);
	lHat = 1 / (s*nthroot(y/2, 3));
end
mParamsSample(3) = lHat;
mParamsSample(4) = uHat;
mParamsSample(5) = sHat;

% 'stretched lognormal';
vParam = lognfit(vc);
[a, b] = lognstat(vParam(1), vParam(2));
mParamsSample(6) = a;
mParamsSample(7) = sqrt(b);

% 'stretched inverse Gaussian';
vParam = mle(vc, 'distribution', 'InverseGaussian');
mParamsSample(8) = vParam(1);
mParamsSample(9) = vParam(2);

% 'lag-exponential 2';
vParam(1) = min(vc);
vParam(2) = 1 / (mean(vc) - vParam(1));
mParamsSample(10) = vParam(1);
mParamsSample(11) = vParam(2);

% 'exp+gauss 2';
lHat = 1 / mean(va);
uHat = mean(vb);
sHat = std(vb, 1);
mParamsSample(12) = uHat;
mParamsSample(13) = sHat;
mParamsSample(14) = lHat;

end