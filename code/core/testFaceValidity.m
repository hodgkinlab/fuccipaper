%{
% Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
%}
function testFaceValidity(mCellData, mCycleMeas)
idxData = initTableCellData();
idxCycle = initTableCycleMeas();

if (any(mCellData(:,idxData.stop) ...
		< mCellData(:,idxData.start)))
	error('smells wrong');
end
if (any(mCycleMeas(:,idxCycle.timeRedOff) ...
		< mCycleMeas(:,idxCycle.timeRedOn)))
	error('smells wrong');
end
if (any(mCycleMeas(:,idxCycle.timeGrnOff) ...
		< mCycleMeas(:,idxCycle.timeGrnOn)))
	error('smells wrong');
end

vTimeStop = mCellData(:,idxData.stop) - mCellData(:,idxData.start);
if (any(vTimeStop < mCycleMeas(:,idxCycle.timeRedOff)))
	error('smells wrong');
end
if (any(vTimeStop < mCycleMeas(:,idxCycle.timeGrnOff)))
	error('smells wrong');
end

end