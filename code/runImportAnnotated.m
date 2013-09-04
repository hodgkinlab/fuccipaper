%{
% Import annotated data (doesn't measure fluorescence peaks)
%
% Andrey Kan, et al.
% akan@wehi.edu.au
% Walter and Eliza Hall Institute, 2012--2013
%
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
clc, clearvars, close all; addpath(genpath('.'));

% specify paths and file names
if (1)
	par.sExpName = '20111118'; bType1 = true;
	%par.sExpName = '20120309'; bType1 = true;
	%par.sExpName = '20120316'; bType1 = true;
	%par.sExpName = '20120413'; bType1 = false;
	
	par.sDataPath = '..\data';
	par.sCurrPath = '..\curr';
	
	par.sCellFile = [par.sExpName, '_division_data.csv'];
	par.sDataFile = [par.sExpName, '-cycle-meas.mat'];
end

% <idxCell> defines the table of cells (row <--> cell)
if (bType1)
	idxCell.start	= 14;
	idxCell.redOn	= 15;
	idxCell.redOff	= 16;
	idxCell.grnOn	= 17;
	idxCell.grnOff	= 18;
	idxCell.end		= 19;
	idxCell.pos		= 22;
	idxCell.col		= 23;
	idxCell.row		= 24;
	idxCell.last	= 24;
	timeOffset		= 0;
else
	idxCell.start	= 14;
	idxCell.redOn	= 15;
	idxCell.redOff	= 16;
	idxCell.grnOn	= 17;
	idxCell.grnOff	= 18;
	idxCell.end		= 19;
	idxCell.pos		= 41;
	idxCell.col		= 42;
	idxCell.row		= 43;
	idxCell.last	= 43;
	timeOffset		= 24;
end

% read table of cells (row <--> cell)
if (1)
	fID = fopen(fullfile(par.sDataPath, par.sCellFile));
	raw = textscan(fID, repmat('%s', 1, idxCell.last), ...
		'delimiter', ',', 'HeaderLines', 1);
	fclose(fID);
	
	% times are measured in hours
	vStart	= str2double(raw{idxCell.start});
	vStop	= str2double(raw{idxCell.end});
	vRedOn	= str2double(raw{idxCell.redOn});
	vRedOff	= str2double(raw{idxCell.redOff});
	vGrnOn	= str2double(raw{idxCell.grnOn});
	vGrnOff	= str2double(raw{idxCell.grnOff});
	
	% mark missed measurements
	vRedOff(vRedOn == 0)	= NaN;
	vRedOff(vRedOff == 0)	= NaN;
	vRedOn(vRedOn == 0)		= NaN;
	
	vStart	= vStart - timeOffset;
	vStop	= vStop - timeOffset;
	vRedOn	= vRedOn - timeOffset;
	vRedOff	= vRedOff - timeOffset;
	vGrnOn	= vGrnOn - timeOffset;
	vGrnOff	= vGrnOff - timeOffset;
	
	vPos = str2double(raw{idxCell.pos});
	vCol = str2double(raw{idxCell.col});
	vRow = str2double(raw{idxCell.row});
	
	maxNumCols = max(vCol);
	maxNumRows = max(vRow);
	vWells = vPos.*(maxNumCols*maxNumRows) + vRow.*maxNumCols + vCol;
	nCells = numel(vStart);
end

% form data structures
if (1)
	idxData		= initTableCellData();
	idxCycle	= initTableCycleMeas();
	mCellData	= NaN(nCells, idxData.last);
	mCycleMeas	= NaN(nCells, idxCycle.last);
	
	mCellData(:,idxData.start)	= vStart;
	mCellData(:,idxData.stop)	= vStop;
	mCellData(:,idxData.dvsn)	= 1;
	mCellData(:,idxData.fate)	= 1;
	mCellData(:,idxData.pos)	= vPos;
	mCellData(:,idxData.well)	= vWells;
	mCellData(:,idxData.mom)	= 1;
	
	mCycleMeas(:,idxCycle.timeRedOn)	= vRedOn - vStart;
	mCycleMeas(:,idxCycle.timeRedOff)	= vRedOff - vStart;
	mCycleMeas(:,idxCycle.timeGrnOn)	= vGrnOn - vStart;
	mCycleMeas(:,idxCycle.timeGrnOff)	= vGrnOff - vStart;
end

% remove wrong records
vMakeSense = (vGrnOff >= vGrnOn) & (vGrnOff <= vStop);
if (~all(vMakeSense))
	warning('there were inconsistent records');
end
mCellData = mCellData(vMakeSense,:);
mCycleMeas = mCycleMeas(vMakeSense,:);

% check and save
testFaceValidity(mCellData, mCycleMeas);
save(fullfile(par.sDataPath, par.sDataFile), 'mCellData', 'mCycleMeas');