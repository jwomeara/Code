% Example 1. Rigid CPD point-set registration. No options are set, so the
% default ones are used. 2D fish point-set.
clear; close all; clc;

fp = fopen('target-rigid.dat', 'rb');
inp = fread(fp, Inf, 'float32');
fclose(fp);
X = zeros(length(inp)/2, 3)
X(:,1) = inp(1:2:end);
X(:,2) = inp(2:2:end);

fp = fopen('input-rigid.dat', 'rb');
inp = fread(fp, Inf, 'float32');
fclose(fp);
Y = zeros(length(inp)/2, 3)
Y(:,1) = inp(1:2:end);
Y(:,2) = inp(2:2:end);

R = [0.9615, -0.2749, 0; 0.2749, 0.9615, 0; 0, 0, 0]
s = 0.5469

Transform=cpd_register(X,Y);

% Initial point-sets
figure,cpd_plot_iter(X(:,1:2), Y(:,1:2)); title('Before');

% Registered point-sets
figure,cpd_plot_iter(X(:,1:2), Transform.Y(:,1:2));  title('After');

% Rotation and scaling errors after the registration
E_R=norm(R-Transform.R)
E_s=norm(s-Transform.s)

if exist('output-rigid.dat', 'file') ~= 0,
    fp = fopen('output-rigid.dat', 'rb');
    inp = fread(fp, Inf, 'float32');
    fclose(fp);
    O = zeros(length(inp)/3, 3)
    O(:,1) = inp(1:3:end);
    O(:,2) = inp(2:3:end);
    O(:,3) = inp(3:3:end);

    figure,cpd_plot_iter(X(:,1:2), O(:,1:2));  title('Output');
end;
