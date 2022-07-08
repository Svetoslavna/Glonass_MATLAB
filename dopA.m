function [output3] = dopA(s1)

tb = 306000;
dt = 45000 - 306000; % ti - tb

% ti = 45000 - round(dt/86400)*86400;

% t = tb:1:ti

ad_x = [-1.3642421e-12,-1.6237011735e-13 ,1.7485470537e-16, -1.0455562943e-20, 5.3011452831e-26];
ad_y = [1.1368684e-12, 1.2870815524e-12, 2.6054733458e-17, -2.2786344334e-20, 1.0112818152e-24];
ad_z = [-1.5916158e-12, -1.3594680937e-13, -1.5930995672e-17, 1.1662419456e-20, -5.5518137243e-25];

a_test_point = struct ('ad_x',ad_x,'ad_y',ad_y,'ad_z',ad_z);
% dopa.ad_x = s1.ad_x(1)+s1.ad_x(2)*(t-tb)+s1.ad_x(3)*(t-tb)^2+s1.ad_x(4)*(t-tb)^3+s1.ad_x(5)*(t-tb)^4;
% dopa.ad_y = s1.ad_y(1)+s1.ad_y(2)*(t-tb)+s1.ad_y(3)*(t-tb)^2+s1.ad_y(4)*(t-tb)^3+s1.ad_y(5)*(t-tb)^4;
% dopa.ad_z = s1.ad_z(1)+s1.ad_z(2)*(t-tb)+s1.ad_z(3)*(t-tb)^2+s1.ad_z(4)*(t-tb)^3+s1.ad_z(5)*(t-tb)^4;

dopa.ax = a_test_point.ad_x(1)+a_test_point.ad_x(2)*(t-tb)+a_test_point.ad_x(3)*(t-tb).^2+a_test_point.ad_x(4)*(t-tb).^3+a_test_point.ad_x(5)*(t-tb).^4;
dopa.ay = a_test_point.ad_y(1)+a_test_point.ad_y(2)*(t-tb)+a_test_point.ad_y(3)*(t-tb).^2+a_test_point.ad_y(4)*(t-tb).^3+a_test_point.ad_y(5)*(t-tb).^4;
dopa.az = a_test_point.ad_z(1)+a_test_point.ad_z(2)*(t-tb)+a_test_point.ad_z(3)*(t-tb).^2+a_test_point.ad_z(4)*(t-tb).^3+a_test_point.ad_z(5)*(t-tb).^4;

output3 = dopa;
end

