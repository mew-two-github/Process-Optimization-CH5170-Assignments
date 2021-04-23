clear; 
%% Given info
minor = [9;14.99];
major = [27.48;-13];
centre = [1.5;2];
c0 = 15;
% minor = [27.99;9.5];
% major = [27.5;-19.65];
% centre = [15;2];
% c0 = 15;
%% Part 1)
% Eigen vectors
v1 = major-centre;
v2 = minor-centre;
% Eigen values
beta = (norm(v2,2)^2)/(norm(v1,2)^2);
V = [v2 v1];
% Minor axis corresponds to the larger eigen value
H_unscaled = V*[1 0;0 beta]*inv(V);
f_min = c0 - 0.5;
d = 1/(v1'*H_unscaled*v1);
% Get H
H = d*H_unscaled;
% Sanity Check
val1 = (v1'*H*v1)/2 + f_min;
val2 = (v2'*H*v2)/2 + f_min;
% Gradient at any point = Hx + a
a = -H*centre;
% At contour function value is c0
c = -0.5*major'*H*major - a'*major + c0;
%% Part 2) and 3)
x = major;
grad = H*x + a;
% Steepest descent: along -1*gradient
sd = -(H*x + a);
% Steepest ascent: along gradient
sa = -sd;
% No change: tangential (perpendicular to gradient)
nc = [grad(2); -grad(1)];
% Will decrease till centre
alpha = grad'*grad/(grad'*H*grad);
%% Part 4)
x = [2;2];
grad = H*x + a;
sd2 = -grad;
sa2 = grad;
nc2 = [grad(2); -grad(1)];
alpha2 = grad'*grad/(grad'*H*grad);
%% Part 5)
H_nd = -H;
a_nd = -H_nd*centre;
c_nd = -0.5*major'*H_nd*major - a_nd'*major + c0;