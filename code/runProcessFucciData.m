%{
% Process data/generate figures for the FUCCI paper
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
clc, clearvars, close all, addpath(genpath('.'));
reset(RandStream.getGlobalStream, 2013);

par.sDataPath = '..\data';
par.sFigPath = '..\data';
par.sExpNames = {'20111118', '20120413', '20120309', '20120316'};
par.sFigFormat = '-dpdf';

groupDataFUCCI(par);
idxData		= initTableCellData();
idxCycle	= initTableCycleMeas();

% layout of figures
if (1)
	vAxesFits		= [];
	vAxesFitsSupp	= [];
	vAxesFitsG1		= [];
	vAxesSisTop		= [];
	vAxesSisBot		= [];
	vAxesSisRed		= [];
	vAxesBars		= [];
	vAxesCorrRG		= [];
	vAxesAHist		= [];
	
	% other plotting parameters
	par.xMax = 30;
	par.yMax = 20;
	par.circle.x0 = 24; par.circle.y0 = 6;
	par.circle.r1 = 1.5; par.circle.r2 = 3;
	par.clr.lightGrn = 'g';% [0.7, 1, 0.7];
	par.clr.darkGrn = [0, 0.4, 0];
	par.text.x = 3.5;
	par.text.y1 = 18;
	par.text.y2 = 16.5;
	par.text.y3 = 15.3;
	par.text.y4 = 13.8;
	
	par.textSize = 6;
	par.letterSize = 8;
	
	set(0, 'defaultaxesfontsize', par.textSize);
	set(0, 'defaulttextfontsize', par.textSize);
end
% fits for S/G2/M: 1.5 column, width = 4.49 inches
if (0)
	width = 4.49; ratio = 0.77;
	fFits = figure('Name', 'fits', 'Position', [5 50 900 900*ratio]);
	set(gcf, ...
		'PaperUnits', 'inches', 'PaperPosition', [0 0 width width*ratio]);
	
	x0 = 0.07; y0 = 0.095;
	x1 = 0.04; y1 = 0.10;
	dx = 0.44; dy = 0.38;
	
	h1 = subplot('Position', [x0, y0 + dy + y1, dx, dy]);
	h2 = subplot('Position', [x0 + dx + x1, y0 + dy + y1, dx, dy]);
	h3 = subplot('Position', [x0, y0, dx, dy]);
	h4 = subplot('Position', [x0 + dx + x1, y0, dx, dy]);
	vAxesFits = [h1, h2, h3, h4];
end
% fits with intercept: 1.5 column, width = 4.49 inches
if (0)
	width = 4.49; ratio = 0.77;
	fSupp = figure('Name', 'fits supp', 'Position', [5 50 900 900*ratio]);
	set(gcf, ...
		'PaperUnits', 'inches', 'PaperPosition', [0 0 width width*ratio])
	
	x0 = 0.07; y0 = 0.095;
	x1 = 0.04; y1 = 0.10;
	dx = 0.44; dy = 0.38;
	
	h1 = subplot('Position', [x0, y0 + dy + y1, dx, dy]);
	h2 = subplot('Position', [x0 + dx + x1, y0 + dy + y1, dx, dy]);
	h3 = subplot('Position', [x0, y0, dx, dy]);
	h4 = subplot('Position', [x0 + dx + x1, y0, dx, dy]);
	vAxesFitsSupp = [h1, h2, h3, h4];
end
% fits for times in red: 1.5 column, width = 4.49 inches
if (0)
	width = 4.49; ratio = 0.75;
	fG1 = figure('Name', 'fits G1', 'Position', [5 50 900 900*ratio]);
	set(gcf, ...
		'PaperUnits', 'inches', 'PaperPosition', [0 0 width width*ratio])
	
	x0 = 0.07; y0 = 0.095;
	x1 = 0.04; y1 = 0.11;
	dx = 0.44; dy = 0.37;
	
	h1 = subplot('Position', [x0, y0 + dy + y1, dx, dy]);
	h2 = subplot('Position', [x0 + dx + x1, y0 + dy + y1, dx, dy]);
	h3 = subplot('Position', [x0, y0, dx, dy]);
	h4 = subplot('Position', [x0 + dx + x1, y0, dx, dy]);
	vAxesFitsG1 = [h1, h2, h3, h4];
end
% Tred VS Tgrn: 1.5 column, width = 4.49 inches
if (1)
	width = 4.49; ratio = 0.75;
	fRG = figure('Name', 'red-grn', 'Position', [5 50 900 900*ratio]);
	set(gcf, ...
		'PaperUnits', 'inches', 'PaperPosition', [0 0 width width*ratio])
	
	x0 = 0.07; y0 = 0.095;
	x1 = 0.04; y1 = 0.07;
	dx = 0.44; dy = 0.41;
	
	h1 = subplot('Position', [x0, y0 + dy + y1, dx, dy]);
	h2 = subplot('Position', [x0 + dx + x1, y0 + dy + y1, dx, dy]);
	h3 = subplot('Position', [x0, y0, dx, dy]);
	h4 = subplot('Position', [x0 + dx + x1, y0, dx, dy]);
	vAxesCorrRG = [h1, h2, h3, h4];
end
% sister correlations: two columns, width = 7.01 inches
if (0)
	width = 7.01; ratio = 0.55;
	fSis = figure('Name', 'sis', 'Position', [5 50 900 900*ratio]);
	set(gcf, ...
		'PaperUnits', 'inches', 'PaperPosition',[0 0 width width*ratio]);
	
	x0 = 0.06;
	y0 = 0.12;
	x1 = 0.04;
	dx = 0.20;
	dy = 0.35;
	y1 = y0 + dy + 0.13;
	
	h1 = subplot('Position', [x0 + 0*(x1 + dx), y1, dx, dy]);
	h2 = subplot('Position', [x0 + 1*(x1 + dx), y1, dx, dy]);
	h3 = subplot('Position', [x0 + 2*(x1 + dx), y1, dx, dy]);
	h4 = subplot('Position', [x0 + 3*(x1 + dx), y1, dx, dy]);
	vAxesSisTop = [h1, h2, h3, h4];
	
	h5 = subplot('Position', [x0 + 0*(x1 + dx), y0, dx, dy]);
	h6 = subplot('Position', [x0 + 1*(x1 + dx), y0, dx, dy]);
	h7 = subplot('Position', [x0 + 2*(x1 + dx), y0, dx, dy]);
	h8 = subplot('Position', [x0 + 3*(x1 + dx), y0, dx, dy]);
	vAxesSisBot = [h5, h6, h7, h8];
end
% sister corr. red: two columns, width = 7.01 inches
if (0)
	width = 7.01; ratio = 0.28;
	fSRed = figure('Name', 'sis red', 'Position', [5 50 900 900*ratio]);
	set(gcf, ...
		'PaperUnits', 'inches', 'PaperPosition', [0 0 width width*ratio])
	x0 = 0.06;
	y0 = 0.15;
	x1 = 0.04;
	dx = 0.20;
	dy = 0.70;
	
	h1 = subplot('Position', [x0 + 0*(x1 + dx), y0, dx, dy]);
	h2 = subplot('Position', [x0 + 1*(x1 + dx), y0, dx, dy]);
	h3 = subplot('Position', [x0 + 2*(x1 + dx), y0, dx, dy]);
	h4 = subplot('Position', [x0 + 3*(x1 + dx), y0, dx, dy]);
	vAxesSisRed = [h1, h2, h3, h4];
end
% bars: single column, width = 3.42 inches
if (0)
	width = 3.42; ratio = 1.53;
	fBars = figure('Name', 'bars', 'Position', [5 50 600 600*ratio]);
	set(gcf, ...
		'PaperUnits', 'inches', 'PaperPosition', [0 0 width width*ratio])
	dx1 = 0.60; dx2 = 0.30;
	dy1 = 0.11; dy2 = 0.138;
	x0 = 0.10; y0 = 0.01; dy = dy1 + dy2;
	
	h1 = subplot('Position', [x0, y0 + 3*dy + dy2, dx1, dy1]);
	h2 = subplot('Position', [x0 + dx1, y0 + 3*dy + dy2, dx2, dy1]);
	h3 = subplot('Position', [x0, y0 + 3*dy, dx1, dy2]);
	h4 = subplot('Position', [x0 + dx1, y0 + 3*dy, dx2, dy2]);
	v1 = [h1, h2, h3, h4];
	
	h1 = subplot('Position', [x0, y0 + 2*(dy1 + dy2) + dy2, dx1, dy1]);
	h2 = subplot('Position', [x0 + dx1, y0 + 2*dy + dy2, dx2, dy1]);
	h3 = subplot('Position', [x0, y0 + 2*(dy1 + dy2), dx1, dy2]);
	h4 = subplot('Position', [x0 + dx1, y0 + 2*(dy1+dy2), dx2, dy2]);
	v2 = [h1, h2, h3, h4];
	
	h1 = subplot('Position', [x0, y0 + dy + dy2, dx1, dy1]);
	h2 = subplot('Position', [x0 + dx1, y0 + dy + dy2, dx2, dy1]);
	h3 = subplot('Position', [x0, y0 + dy, dx1, dy2]);
	h4 = subplot('Position', [x0 + dx1, y0 + dy, dx2, dy2]);
	v3 = [h1, h2, h3, h4];
	
	h1 = subplot('Position', [x0, y0 + 0*(dy1 + dy2) + dy2, dx1, dy1]);
	h2 = subplot('Position', [x0 + dx1, y0 + dy2, dx2, dy1]);
	h3 = subplot('Position', [x0, y0 + 0*(dy1 + dy2), dx1, dy2]);
	h4 = subplot('Position', [x0 + dx1, y0, dx2, dy2]);
	v4 = [h1, h2, h3, h4];
	vAxesBars = {v1, v2, v3, v4};
end
% alpha-plot histograms
if (0)
	width = 4.5;
	ratio = 1;
	fAHist = figure('Position', [50 50 1000 1000*ratio]);
	set(gcf, 'PaperUnits', 'inches', ...
		'PaperPosition', [0 0 width width*ratio]);
	h1 = subplot(4, 3, 1);
	h2 = subplot(4, 3, 2);
	h3 = subplot(4, 3, 3);
	v1 = [h1, h2, h3];
	h1 = subplot(4, 3, 4);
	h2 = subplot(4, 3, 5);
	h3 = subplot(4, 3, 6);
	v2 = [h1, h2, h3];
	h1 = subplot(4, 3, 7);
	h2 = subplot(4, 3, 8);
	h3 = subplot(4, 3, 9);
	v3 = [h1, h2, h3];
	h1 = subplot(4, 3, 10);
	h2 = subplot(4, 3, 11);
	h3 = subplot(4, 3, 12);
	v4 = [h1, h2, h3];
	vAxesAHist = {v1, v2, v3, v4};
end

% used for computing bootstrap confidence intervals
par.nSamples = 1000;

vPropSG2M = [0.73 0.78 0.65 0.72];
vPropG1 = 1 - vPropSG2M;
sFigNames = {'CpG', 'aCD40', 'CD8', 'OT-I'};
iDataset = 1;
for iExp = 1:numel(par.sExpNames)
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
	for iGroup = 1:nGroups
		par.iDataset = iDataset;
		par.sTitle = groupNames{iGroup}; %#ok<SUSENS>
		fprintf('=== %s ', par.sTitle);
		% select a group of cells
		if (1)
			vSelect = (vGroups == iGroup);
			nCells = sum(vSelect);
			if (nCells == 0)
				error('empty selection');
			else
				fprintf('(%d cells) ===\n', nCells);
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
			mCellDataX	= mCellData(vSingle,:);
			mCycleMeasX = mCycleMeas(vSingle,:);
			
			% split siblings
			mCellDataY	= mCellData(~vSingle,:);
			mCycleMeasY = mCycleMeas(~vSingle,:);
			mCellDataA	= mCellDataY(1:2:end,:);
			mCycleMeasA = mCycleMeasY(1:2:end,:);
			mCellDataB	= mCellDataY(2:2:end,:);
			mCycleMeasB = mCycleMeasY(2:2:end,:);
			
			nUnpaired = size(mCellDataX, 1);
			nPairs = size(mCellDataA, 1);
			fprintf('#pairs = %d; #unpaired = %d; ratio = %0.2f\n\n', ...
				nPairs, nUnpaired, (double(nUnpaired)/double(2*nPairs)));
		end
		% alpha plots and comparison with SM models
		if (1)
			vc = mCellDataA(:,idxData.stop) - mCellDataA(:,idxData.start);
			va = mCycleMeasA(:,idxCycle.timeGrnOn);
			vb = vc - va;
			va = va(:);
			vb = vb(:);
			
			vProbStep = 1/double(numel(vc));
			vt = 0:0.01:35;
			xLimA = 10;
			xLimB = 20;
			%xLimC = 30;
			
			models = {};
			
			[~, x, flo, fup] = ecdf(vc);
			%[~, xa, floa, fupa] = ecdf(va); % this affects random stream
			%[~, xb, flob, fupb] = ecdf(vb);
			xa = []; xb = []; floa = []; fupa = []; flob = []; fupb = [];
			if (1)
				model = []; model.name = 'conf. interval';
				model.clr = 'b';
				model.style = '--';
				
				model.vt = x;
				model.vPred = 1 - flo;
				
				model.vtA = xa;
				model.vtB = xb;
				model.vPredA = 1 - floa;
				model.vPredB = 1 - flob;
				models = [models, model];
			end
			if (1)
				model = []; model.name = 'conf. interval';
				model.clr = 'b';
				model.style = '--';
				
				model.vt = x;
				model.vPred = 1 - fup;
				
				model.vtA = xa;
				model.vtB = xb;
				model.vPredA = 1 - fupa;
				model.vPredB = 1 - fupb;
				models = [models, model];
			end
			
			model = []; model.name = 'lag-exponential';
			if (1)
				model.clr = 'r';
				model.style = '--';
				
				%vParam(1) = min(vc);
				%vParam(2) = 1 / (mean(vc) - vParam(1));
				
 				[f, x] = ecdf(vc);
 				x = x(f < 1); f = f(f < 1); f = log(1 - f);
 				p = polyfit(x, f, 1);
 				vParam(1) = -p(2)/p(1); vParam(2) = -p(1);
				
				%vParam = calcOptPatams(vc, iExp);
				model.vParam = vParam;
				
				model.vt = vt;
				model.vPred = exp(-vParam(2)*(vt - vParam(1)));
				
				model.vtA = vt;
				model.vPredA = exp(-vParam(2)*vt);
				
				model.vtB = vt;
				model.vPredB = (vt <= vParam(1)) + realmin;
				
				%lambda = vParam(2);
				%t0 = vParam(1);
				%like = NaN; % prod(lambda*exp(-lambda*(vc - t0)));
				%model.aic = 2*2 - 2*log(like);
				[~, model.pa] = kstest(va, [va, 1 - exp(-vParam(2)*va)]);
				[~, model.pb] = kstest(vb, [vb, (vb >= vParam(1))]);
				[~, model.pc] = kstest(vc, ...
					[vc, 1 - exp(-vParam(2)*(vc - vParam(1)))]);
				
				%round(100.*vParam)./100
				models = [models, model];
			end
			
			% NOTE: this re-uses vParam
			model = []; model.name = 'stretched lag-exponential';
			if (0)
				model.clr = 'm';
				model.style = '-';
				
				vt = 0:1:40;
				model.vt = vt;
				model.vPred = exp(-vParam(2)*(vt - vParam(1)));
				
				model.vtA = vt.*vPropG1(iExp);
				model.vPredA = model.vPred;
				
				model.vtB = vt.*vPropSG2M(iExp);
				model.vPredB = model.vPred;
				
				%lambda = vParam(2);
				%t0 = vParam(1);
				%like = prod(lambda*exp(-lambda*(vc - t0)));
				%model.aic = 2*3 - 2*log(like);
				[~, model.pa] = kstest(va, ...
					[va, 1 - exp(-vParam(2)*(va./vPropG1(iExp) - vParam(1)))]);
				[~, model.pb] = kstest(vb, ...
					[vb, 1 - exp(-vParam(2)*(vb./vPropSG2M(iExp) - vParam(1)))]);
				[~, model.pc] = kstest(vc, ...
					[vc, 1 - exp(-vParam(2)*(vc - vParam(1)))]);
				
				models = [models, model];
			end
			
			model = []; model.name = 'exp+gauss';
			if (0)
				model.clr = 'r';
				model.style = '-';
				
%   				m = mean(vc);
%   				s = std(vc);
%   				y = skewness(vc);
%   				uHat = m - s*nthroot(y/2, 3);
%   				sHat = s*sqrt(1 - nthroot(y/2, 3).^2);
%   				lHat = 1 / (s*nthroot(y/2, 3));
				
				lHat = 1 / mean(va);
				uHat = mean(vb);
				sHat = std(vb, 1);
				
				model.vParam = [uHat, sHat, lHat];
				
				u = lHat*(vt - uHat);
				v = lHat*sHat;
				
				model.vt = vt;
				model.vPred = 1 - normcdf(u, 0, v) ...
					+ exp(-u + v*v/2 + log(normcdf(u, v*v, v)));
				
				model.vtA = vt;
				model.vPredA = exp(-lHat*vt);
				
				model.vtB = vt;
				model.vPredB = 1 - normcdf(vt, uHat, sHat);
				
				x = vc;
				
				[~, model.pa] = kstest(va, [va, 1 - exp(-lHat*va)]);
				[~, model.pb] = kstest(vb, [vb, normcdf(vb, uHat, sHat)]);
				
				u = lHat*(vc - uHat);
				v = lHat*sHat;
				[~, model.pc] = kstest(vc, [vc, ...
					normcdf(u, 0, v) - exp(-u + v*v/2 + log(normcdf(u, v*v, v)))]);
				
				%round(100.*[lHat uHat sHat] )./100
				models = [models, model];
			end
			
			model.name = 'stretched lognormal';
			if (0)
				model.clr = 'g';
				model.style = '-';
				vParam = lognfit(vc);
				model.vParam = [vParam, vPropG1(iExp), vPropSG2M(iExp)];
				model.vt = vt;
				model.vPred = 1 - logncdf(vt, vParam(1), vParam(2));
				model.vtA = vt.*vPropG1(iExp);
				model.vtB = vt.*vPropSG2M(iExp);
				model.vPredA = model.vPred;
				model.vPredB = model.vPred;
				
				
				%like = prod(lognpdf(va./vPropG1(iExp), vParam(1), vParam(2)));
				%model.aic = 2*3 - 2*log(like);
				
				%fThr = lognpdf(xi./vPropG1(iExp), vParam(1), vParam(2));
				%b = fThr;
				%hold on; plot(xi, fThr, 'color', 'g');
				% 				trapz(xi, sqrt(fEmpA.*fThr))
				% 				trapz(xi, log(fEmpA./fThr).*fEmpA)
				[~, model.pa] = kstest(va, ...
					[va, logncdf(va./vPropG1(iExp), vParam(1), vParam(2))]);
				[~, model.pb] = kstest(vb, ...
					[vb, logncdf(vb./vPropSG2M(iExp), vParam(1), vParam(2))]);
				[~, model.pc] = kstest(vc, ...
					[vc, logncdf(vc, vParam(1), vParam(2))]);
				
				%[a, b] = lognstat(vParam(1), vParam(2)); round(100.*[a, sqrt(b)])./100
				models = [models, model];
			end
			
			model.name = 'sum of lognormals';
			if (0)
				model.clr = 'c';
				model.style = '-';
				
				vParam = lognfit(va);
				m1 = vParam(1);
				s1 = vParam(2);
				
				vParam = lognfit(vb);
				m2 = vParam(1);
				s2 = vParam(2);
				
				model.vtA = vt;
				model.vPredA = 1 - logncdf(vt, m1, s1);
				
				model.vtB = vt;
				model.vPredB = 1 - logncdf(vt, m2, s2);
				
				model.vt = 0:1:50;
				nPts = numel(model.vt);
				v = conv( ...
					lognpdf(model.vt, m1, s1), lognpdf(model.vt, m2, s2));
				v = v(1:nPts);
				for i = 1:nPts
					model.vPred(i) = trapz(v(i:end));
				end
				model.vPred = round(10000*model.vPred)/10000;
				%like = prod(lognpdf(vPartA, m1, s1));
				%like = like*prod(lognpdf(vb, m2, s2));
				%like = prod(interp1(model.vt, v, vc));
				%model.aic = 2*4 - 2*log(like);
				[~, model.pa] = kstest(va, [va, logncdf(va, m1, s1)]);
				[~, model.pb] = kstest(vb, [vb, logncdf(vb, m2, s2)]);
				[~, model.pc] = kstest(vc, [model.vt(:), 1 - model.vPred(:)]);
				
				%[a, b] = lognstat(m1, s1); round(100.*[a, b])./100
				%[a, b] = lognstat(m2, s2); round(100.*[a, b])./100
				models = [models, model];
			end
			
			model.name = 'stretched inverse Gaussian';
			if (0)
				model.clr = 'k';
				model.style = '-';
				
				vParam = mle(vc, 'distribution', 'InverseGaussian');
				invgcdf = @(x, u, l) ...
					normcdf(sqrt(l./x).*(x/u - 1), 0, 1) ...
					+ exp(2*l/u)*normcdf(-sqrt(l./x).*(x/u + 1), 0, 1);
				model.vt = vt;
				model.vPred = 1 - invgcdf(vt, vParam(1), vParam(2));
				model.vPred = round(10000*model.vPred)/10000;
				model.vtA = vt.*vPropG1(iExp);
				model.vtB = vt.*vPropSG2M(iExp);
				model.vPredA = model.vPred;
				model.vPredB = model.vPred;
				
				% 				invgpdf = @(x, u, l) ...
				% 					sqrt(l./(2*pi*x.^3)).*exp(-(l*(x - u).^2)./(2*u*u*x));
				% 				like = prod(invgpdf(vc, vParam(1), vParam(2)));
				% 				model.aic = 2*3 - 2*log(like);
				f = invgcdf(va./vPropG1(iExp), vParam(1), vParam(2));
				f = round(10000*f)/10000;
				[~, model.pa] = kstest(va, [va, f]);
				tx = vt.*vPropSG2M(iExp);
				ty = invgcdf(vt, vParam(1), vParam(2));
				[~, model.pb] = kstest(vb, [tx(:), ty(:)]);
				[~, model.pc] = kstest(vc, ...
					[model.vt(:), 1 - model.vPred(:)]);
				
				%round(100.*vParam)./100
				models = [models, model];
			end
			
			model.name = 'stretched reciprocal normal';
			if (0)
				model.clr = 'k';
				model.style = '-';
				
				u = mean(1./vc);
				s = std(1./vc);
				
				model.vt = vt;
				model.vPred = 1 - invgcdf(vt, vParam(1), vParam(2));
				model.vPred = round(10000*model.vPred)/10000;
				model.vtA = vt.*vPropG1(iExp);
				model.vtB = vt.*vPropSG2M(iExp);
				model.vPredA = model.vPred;
				model.vPredB = model.vPred;
				
				% 				invgpdf = @(x, u, l) ...
				% 					sqrt(l./(2*pi*x.^3)).*exp(-(l*(x - u).^2)./(2*u*u*x));
				% 				like = prod(invgpdf(vc, vParam(1), vParam(2)));
				% 				model.aic = 2*3 - 2*log(like);
				f = invgcdf(va./vPropG1(iExp), vParam(1), vParam(2));
				f = round(10000*f)/10000;
				[~, model.pa] = kstest(va, [va, f]);
				tx = vt.*vPropSG2M(iExp);
				ty = invgcdf(vt, vParam(1), vParam(2));
				[~, model.pb] = kstest(vb, [tx(:), ty(:)]);
				[~, model.pc] = kstest(vc, ...
					[model.vt(:), 1 - model.vPred(:)]);
				
				%round(100.*vParam)./100
				models = [models, model];
			end
			
			lw = 1;
			msz = 2.3;
			width = 4.5;
			ratio = 0.5;
			f = figure('Position', [50 130 900 900*ratio]);
			set(gcf, 'PaperUnits', 'inches', 'PaperPosition', ...
				[0 0 width width*ratio]);
			vxA = sort(va);
			vxB = sort(vb);
			vx = sort(vc);
			%vxB = vxB + 10;
			%vx = vx + 30;
			
			subplot(1, 3, 1);
 			semilogy( ...
 				vxA, 1 - (vProbStep:vProbStep:1), ...
 				'LineStyle', 'none','MarkerSize', msz, 'Marker', 'o');
			subplot(1, 3, 2);
 			semilogy( ...
 				vxB, 1 - (vProbStep:vProbStep:1), ...
 				'LineStyle', 'none','MarkerSize', msz, 'Marker', 'o');
			subplot(1, 3, 3);
 			semilogy( ...
 				vx, 1 - (vProbStep:vProbStep:1), ...
 				'LineStyle', 'none','MarkerSize', msz, 'Marker', 'o');
			
			nModels = numel(models);
			handles = [];
			names = {};
			%vAIC = [];
			
			mResults = NaN(nModels, 3);
			for i = 1:nModels
				model = models{i};
				%model.vtB = model.vtB + 10;
				%model.vt = model.vt + 30;
				
				subplot(1, 3, 1);
				hold on;
				vSelect = (model.vtA <= xLimA);
				semilogy(model.vtA(vSelect), model.vPredA(vSelect), ...
					'MarkerSize', msz, 'LineWidth', lw, ...
					'Color', model.clr, 'LineStyle', model.style);
				
				hold on;
				text(7, 0.7, par.sTitle);
				%text(7, 0.5, sprintf('N = %d', numel(vc)));
				xlabel('time (h)');
				ylabel('proporiton');
				
				subplot(1, 3, 2);
				hold on;
				vSelect = (model.vtB <= xLimB);
				semilogy(model.vtB(vSelect), model.vPredB(vSelect), ...
					'MarkerSize', msz, 'LineWidth', lw, ...
					'Color', model.clr, 'LineStyle', model.style);
				
				subplot(1, 3, 3);
				hold on;
				h = semilogy(model.vt, model.vPred, ...
					'MarkerSize', msz, 'LineWidth', lw, ...
					'Color', model.clr, 'LineStyle', model.style);
				
				handles = [handles, h];
				names = [names, model.name];
				%vAIC = [vAIC, model.aic];
				if (isfield(model, 'pa'))
					fprintf('KS: %0.5f, %0.5f, %0.5f (%s)\n', ...
						model.pa, model.pb, model.pc, model.name);
					mResults(i,:) = [model.pa, model.pb, model.pc];
				end
			end
			mResults = round(10000*mResults)/10000;
			
			subplot(1, 3, 1);
			hold on;
			xlim([0 10]);
			ylim([0.001 1]);
			subplot(1, 3, 2);
			hold on;
			xlim([0 20]);
			ylim([0.001 1]);
			subplot(1, 3, 3);
			hold on;
			xlim([0 30]);
			ylim([0.001 1]);
			%xlim([0 60]);
			%set(gca, 'XTick', [0:10:60]);
			%set(gca, 'XTickLabel', {'0', '10', '10', '20', '10', '20', '30', '40'});
			legend(handles, names, 'Location', 'best');
			print(f, '-dpdf', '-r900', ...
				fullfile('../data/', sFigNames{iExp}));
			
			plotAlphaHist(vAxesAHist, mCellData, mCycleMeas, models, par);
		end
		% analyse data (for the supplementary section)
		if (1)
			T = mCellData(:,idxData.stop) - mCellData(:,idxData.start);
			fprintf('mean time to div = %0.2f\n', mean(T));
			fprintf('st. dev. time to div = %0.2f\n', std(T));
			
			vParHat = lognfit(T);
			[uLog, vLog] = lognstat(vParHat(1), vParHat(2));
			sLog = sqrt(vLog);
			fprintf('MLE lognormal uLog = %0.2f; sLog = %0.2f\n', uLog, sLog);
			
			b2 = T - mCycleMeas(:,idxCycle.timeGrnOn);
			ab1 = mCycleMeas(:,idxCycle.timeGrnOn);
			
			vSelect = ~isnan(b2);
			T = T(vSelect);
			b2 = b2(vSelect);
			ab1 = ab1(vSelect);
			fprintf('analytics: n = %d\n', sum(vSelect));
			
			s = std(T);
			sZ = std(b2);
			fprintf('s(b2) = %.02f; s(T) = %.02f; ratio = %.02f\n', ...
				sZ, s, sZ/s);
			
			uaUpBnd = sqrt(s*s - sZ*sZ);
			fprintf('u(A) <= %.02f (hours)\n', uaUpBnd);
			
			m = mean(T);
			y = skewness(T);
			
			uHat = m - s*nthroot(y/2, 3);
			sHat = s*s*(1 - nthroot(y/2, 3).^2);
			kHat = s*nthroot(y/2, 3);
				
			fprintf('EMG estimates: u = %.02f; s = %.02f; k = %.02f\n', ...
				uHat, sHat, kHat);
			fprintf('recall that: s >= %.02f; k <= %.02f\n', ...
				sZ, uaUpBnd);
			
			cv = cov(ab1(:), b2(:));
			fprintf('cov(ab1, b2) = %.02f\n', cv(1, 2));
		end
		
		% produce other figures
		if (~isempty(vAxesFits))
			par.hAxes = vAxesFits(iDataset);
			plotFigureFits( ...
				mCellDataA, mCellDataB, mCellDataX, ...
				mCycleMeasA, mCycleMeasB, mCycleMeasX, par);
			set(par.hAxes, 'XTick', [0 10 20 30]);
			set(par.hAxes, 'YTick', [0 5 10 15 20]);
			if (iDataset == 3)
				axes(vAxesFits(iDataset));
				xlabel('total time to division, T_{div} (h)');
				ylabel('duration of S/G2/M phases, T_{grn} (h)');
			end
		end
		if (~isempty(vAxesFitsSupp))
			par.hAxes = vAxesFitsSupp(iDataset);
			plotFigureFitsSupp( ...
				mCellDataA, mCellDataB, mCellDataX, ...
				mCycleMeasA, mCycleMeasB, mCycleMeasX, par);
			set(par.hAxes, 'XTick', [0 10 20 30]);
			set(par.hAxes, 'YTick', [0 5 10 15 20]);
			if (iDataset == 3)
				axes(vAxesFitsSupp(iDataset));
				xlabel('total time to division, T_{div} (h)');
				ylabel('duration of S/G2/M phases, T_{grn} (h)');
			end
		end
		if (~isempty(vAxesFitsG1))
			par.hAxes = vAxesFitsG1(iDataset);
			plotFigureFitG1( ...
				mCellDataA, mCellDataB, mCellDataX, ...
				mCycleMeasA, mCycleMeasB, mCycleMeasX, par);
			set(par.hAxes, 'XTick', [0 10 20 30]);
			set(par.hAxes, 'YTick', [0 5 10 15]);
			if (iDataset == 3)
				axes(vAxesFitsG1(iDataset));
				xlabel('total time to division, T_{div} (h)');
				ylabel('time to red off, T_{red} (h)');
			end
		end
		if (~isempty(vAxesCorrRG))
			par.hAxes = vAxesCorrRG(iDataset);
			plotFigureCorrRG( ...
				mCellDataA, mCellDataB, mCellDataX, ...
				mCycleMeasA, mCycleMeasB, mCycleMeasX, par);
			set(par.hAxes, 'XTick', [0 5 10 15 20]);
			set(par.hAxes, 'YTick', [0 5 10 15]);
			if (iDataset == 3)
				axes(vAxesCorrRG(iDataset));
				xlabel('duration of S/G_{2}/M phases, T_{grn} (h)');
				ylabel('time to red off, T_{red} (h)');
			end
		end
		if (~isempty(vAxesSisTop))
			par.hAxesTop = vAxesSisTop(iDataset);
			par.hAxesBot = vAxesSisBot(iDataset);
			plotFigureSis( ...
				mCellDataA, mCellDataB, mCycleMeasA, mCycleMeasB, par);
			set(par.hAxesTop, 'XTick', 0:5:25);
			set(par.hAxesTop, 'YTick', 0:5:25);
			set(par.hAxesBot, 'XTick', 0:5:25);
			set(par.hAxesBot, 'YTick', 0:5:25);
			if (iDataset == 1)
				axes(par.hAxesTop);
				xlabel('sib. 1, total time to divide, T_{div} (h)');
				ylabel('sib. 2, total time to divide, T_{div} (h)');
				axes(par.hAxesBot);
				xlabel('sib. 1, time in S/G2/M phases, T_{grn} (h)');
				ylabel('sib. 2, time in S/G2/M phases, T_{grn} (h)');
			end
		end
		if (~isempty(vAxesSisRed))
			par.hAxes = vAxesSisRed(iDataset);
			plotFigureSisRed(mCycleMeasA, mCycleMeasB, par);
			set(par.hAxes, 'XTick', 0:5:25);
			set(par.hAxes, 'YTick', 0:5:25);
			if (iDataset == 1)
				if (~isempty(vAxesSisRed))
					axes(par.hAxes);
					xlabel('sib. 1, time to red off, T_{red} (h)');
					ylabel('sib. 2, time to red off, T_{red} (h)');
				end
			end
		end
		if (~isempty(vAxesBars))
			par.hAxes = vAxesBars{iDataset};
			plotFigureBars( ...
				mCellDataA, mCellDataB, mCycleMeasA, mCycleMeasB, par);
			set(par.hAxes(1), 'XTick', []);
			set(par.hAxes(2), 'XTick', []);
			
			set(par.hAxes(1), 'YTick', []);
			set(par.hAxes(2), 'YTick', []);
			set(par.hAxes(3), 'YTick', []);
			set(par.hAxes(4), 'YTick', []);
			
			axes(par.hAxes(1)); axis off;
			axes(par.hAxes(2)); axis off;
			axes(par.hAxes(3)); axis off;
			axes(par.hAxes(4)); axis off;
			
			if (iDataset == 4)
				yAxes = -0.7;
				yTick = 0.08;
				yText = 0.25;
				axes(par.hAxes(3));
				vLimY = ylim();
				ylim([-1.5, vLimY(end)]);
				line([-20, 20], [yAxes, yAxes], 'Color', 'k');
				if (1)
					line([-20, -20], [yAxes, yAxes + yTick], 'Color', 'k');
					line([20, 20], [yAxes, yAxes + yTick], 'Color', 'k');
					text(-20 - 1.6, yAxes - yText, '20');
					text(20 - 1.6, yAxes - yText, '20');
					
					line([-10, -10], [yAxes, yAxes + yTick], 'Color', 'k');
					line([10, 10], [yAxes, yAxes + yTick], 'Color', 'k');
					text(-10 - 1.6, yAxes - yText, '10');
					text(10 - 1.6, yAxes - yText, '10');
					
					line([0 0], [yAxes, yAxes + yTick], 'Color', 'k');
					text(-0.6, yAxes - yText, '0');
				end
				text(-3, yAxes - 0.7, 'time (h)');
				
				axes(par.hAxes(4));
				vLimY = ylim();
				ylim([-1.5, vLimY(end)]);
				line([-1, 1], [yAxes, yAxes], 'Color', 'k');
				if (1)
					line([-1, -1], [yAxes, yAxes + yTick], 'Color', 'k');
					line([1, 1], [yAxes, yAxes + yTick], 'Color', 'k');
					text(-1 - 0.05, yAxes - yText, '1');
					text(1 - 0.05, yAxes - yText, '1');
					
					line([0 0], [yAxes, yAxes + yTick], 'Color', 'k');
					text(-0.05, yAxes - yText, '0');
				end
				text(-0.7, yAxes - 0.7, 'normalized time');
			end
		end
		iDataset = iDataset + 1;
	end
end

% save figures
if (1)
	if (~isempty(vAxesFits))
		print(fFits, par.sFigFormat, '-r600', fullfile(par.sFigPath, 'f2'));
	end
	if (~isempty(vAxesBars))
		print(fBars, par.sFigFormat, '-r600', fullfile(par.sFigPath, 'f4'));
	end
	if (~isempty(vAxesSisTop))
		print(fSis, par.sFigFormat, '-r600', fullfile(par.sFigPath, 'f5'));
	end
	if (~isempty(vAxesFitsSupp))
		print(fSupp, par.sFigFormat, '-r600', fullfile(par.sFigPath, 'fs1'));
	end
	if (~isempty(vAxesFitsG1))
		print(fG1, par.sFigFormat, '-r600', ...
			fullfile(par.sFigPath, 'fits-g1'));
	end
	if (~isempty(vAxesCorrRG))
		print(fRG, par.sFigFormat, '-r600', ...
			fullfile(par.sFigPath, 'corr-rg'));
	end
	if (~isempty(vAxesSisRed))
		print(fSRed, par.sFigFormat, '-r600', ...
			fullfile(par.sFigPath, 'sis-red'));
	end
	if (~isempty(vAxesAHist))
		print(fAHist, par.sFigFormat, '-r600', ...
			fullfile(par.sFigPath, 'alpha-hist'));
	end
end