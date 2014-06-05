% Example 1. Rigid CPD point-set registration. No options are set, so the
% default ones are used. 2D fish point-set.
clear; close all; clc;

fp = fopen('target-affine.dat', 'rb');
inp = fread(fp, Inf, 'float32');
fclose(fp);
X = zeros(length(inp)/2, 2)
X(:,1) = inp(1:2:end);
X(:,2) = inp(2:2:end);

fp = fopen('input-affine.dat', 'rb');
inp = fread(fp, Inf, 'float32');
fclose(fp);
Y = zeros(length(inp)/2, 2)
Y(:,1) = inp(1:2:end);
Y(:,2) = inp(2:2:end);

opt.method = 'affine'

Transform=cpd_register(X,Y,opt);

% Initial point-sets
figure,cpd_plot_iter(X(:,1:2), Y(:,1:2)); title('Before');

% Registered point-sets
figure,cpd_plot_iter(X(:,1:2), Transform.Y(:,1:2));  title('After');

if exist('output-affine.dat', 'file') ~= 0,
    fp = fopen('output-affine.dat', 'rb');
    inp = fread(fp, Inf, 'float32');
    fclose(fp);
    O = zeros(length(inp)/3, 3)
    O(:,1) = inp(1:3:end);
    O(:,2) = inp(2:3:end);
    O(:,3) = inp(3:3:end);

    figure,cpd_plot_iter(X(:,1:2), O(:,1:2));  title('Output');
end;