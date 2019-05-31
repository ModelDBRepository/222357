% compute_centroid.m
% a renamed centroid function because a variable centroid was assigned in
% RPSTATIC.m
%
% [centx, centy]=compute_centroid(theta, rho);
% Returns the cartesian centroid of a polar plot

function [centx, centy] = compute_centroid(xcenters, FR_breath_angles)

% create x, y vectors for the histogram
u=zeros(1,length(xcenters));
v=zeros(1,length(xcenters));

for i=1:length(xcenters)
    u(i)=cos(xcenters(i))*FR_breath_angles(i);
    v(i)=sin(xcenters(i))*FR_breath_angles(i);
end
centx=mean(u);
centy=mean(v);
