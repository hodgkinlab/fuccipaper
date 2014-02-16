%{
% Fit different distributions to the total division time
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
function runTotDivDist()
clc, clearvars, close all, addpath(genpath('.'));
reset(RandStream.getGlobalStream, 2013);
par.sDataPath = '..\data';
par.sExpNames = {'20111118', '20120413', '20120309', '20120316'};

par.nSamples = 1000;% used for computing bootstrap confidence intervals
par.nParams = 14;

groupDataFUCCI(par);

idxData		= initTableCellData();
idxCycle	= initTableCycleMeas();
par.idxData = idxData;
par.idxCycle = idxCycle;

paramVal = NaN(numel(par.sExpNames), par.nParams);

vPropSG2M = [0.73 0.78 0.65 0.72];
vPropG1 = 1 - vPropSG2M;

iDataset = 1;
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
		
		estimateParams(par, vPropG1(iExp));
				
		iDataset = iDataset + 1;
	end
end
end

function mParamsSample = estimateParams(par, kG1)
mParamsSample = NaN(1, par.nParams);

idxData = par.idxData;
idxCycle = par.idxCycle;
mCellDataA = par.mCellDataA;
mCycleMeasA = par.mCycleMeasA;

vc = mCellDataA(:,idxData.stop) - mCellDataA(:,idxData.start);
va = mCycleMeasA(:,idxCycle.timeGrnOn);

vb = vc - va;
va = va(:);
vb = vb(:);
t = linspace(0, max(vc), 50);

figure;
[f, x] = ecdf(vc);
plot(x, f, 'LineStyle', 'none', 'Marker', '*');


% stretched lognormal
par = lognfit(vc);
y = logncdf(t, par(1), par(2));
hold on; plot(t, y, 'color', 'k');

[~, pa] = kstest(va, [va, logncdf(va./kG1, par(1), par(2))]);
[~, pb] = kstest(vb, [vb, logncdf(vb./(1 - kG1), par(1), par(2))]);
[~, pc] = kstest(vc, [vc, logncdf(vc, par(1), par(2))]);

[m, v] = lognstat(par(1), par(2));
fprintf('lognormal, mean = %0.3f, std = %0.3f\n', m, sqrt(v));
fprintf(' pa = %0.5f\n', pa);
fprintf(' pb = %0.5f\n', pb);
fprintf(' pc = %0.5f\n', pc);


% stretched gamma
par = gamfit(vc);
y = gamcdf(t, par(1), par(2));
hold on; plot(t, y, 'color', 'r');

[~, pa] = kstest(va, [va, gamcdf(va./kG1, par(1), par(2))]);
[~, pb] = kstest(vb, [vb, gamcdf(vb./(1 - kG1), par(1), par(2))]);
[~, pc] = kstest(vc, [vc, gamcdf(vc, par(1), par(2))]);

[m, v] = gamstat(par(1), par(2));
fprintf('gamma, mean = %0.3f, std = %0.3f\n', m, sqrt(v));
fprintf(' pa = %0.5f\n', pa);
fprintf(' pb = %0.5f\n', pb);
fprintf(' pc = %0.5f\n', pc);


% stretched weibull
par = wblfit(vc);
y = wblcdf(t, par(1), par(2));
hold on; plot(t, y, 'color', 'g');

[~, pa] = kstest(va, [va, wblcdf(va./kG1, par(1), par(2))]);
[~, pb] = kstest(vb, [vb, wblcdf(vb./(1 - kG1), par(1), par(2))]);
[~, pc] = kstest(vc, [vc, wblcdf(vc, par(1), par(2))]);

[m, v] = wblstat(par(1), par(2));
fprintf('weib, mean = %0.3f, std = %0.3f\n', m, sqrt(v));
fprintf(' pa = %0.5f\n', pa);
fprintf(' pb = %0.5f\n', pb);
fprintf(' pc = %0.5f\n', pc);


% stretched normal
[par1, par2] = normfit(vc);
par = [par1, par2];
y = normcdf(t, par(1), par(2));
hold on; plot(t, y, 'color', 'c');

[~, pa] = kstest(va, [va, normcdf(va./kG1, par(1), par(2))]);
[~, pb] = kstest(vb, [vb, normcdf(vb./(1 - kG1), par(1), par(2))]);
[~, pc] = kstest(vc, [vc, normcdf(vc, par(1), par(2))]);

fprintf('normal, mean = %0.3f, std = %0.3f\n', par1, par2);
fprintf(' pa = %0.5f\n', pa);
fprintf(' pb = %0.5f\n', pb);
fprintf(' pc = %0.5f\n', pc);


% stretched inverse Gaussian
par = mle(vc, 'distribution', 'InverseGaussian');
invgcdf = @(x, u, l) ...
	normcdf(sqrt(l./x).*(x/u - 1), 0, 1) ...
	+ exp(2*l/u)*normcdf(-sqrt(l./x).*(x/u + 1), 0, 1);

y = invgcdf(t, par(1), par(2));
hold on; plot(t, y, 'color', 'y');

f = invgcdf(va./kG1, par(1), par(2));
f = round(10000*f)/10000;
[~, pa] = kstest(va, [va, f]);
f = invgcdf(vb./(1 - kG1), par(1), par(2));
f = round(10000*f)/10000;
[~, pb] = kstest(vb, [vb, f]);
f = invgcdf(vc, par(1), par(2));
f = round(10000*f)/10000;
[~, pc] = kstest(vc, [vc, f]);

m = par(1);
s = sqrt(par(1)^3/par(2));
fprintf('inv-gauss, mean = %0.3f, std = %0.3f\n', m, s);
fprintf(' pa = %0.5f\n', pa);
fprintf(' pb = %0.5f\n', pb);
fprintf(' pc = %0.5f\n', pc);

end