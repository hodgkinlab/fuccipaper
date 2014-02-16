clc; clear variables; close all;
load('../data/bounds.mat', 'bounds');
load('../data/moments.mat', 'moments');

for i = 1:4
	a = bounds{i};
	b = 1./moments{i};
	[p, h] = ranksum(a, b);
	fprintf('u(a) = %0.3f; u(b) = %0.3f\n', nanmean(a), nanmean(b));
	fprintf('p = %0.3f; h = %d\n', p, h);
end