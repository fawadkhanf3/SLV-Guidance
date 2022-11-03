function [R,V] = Orb2State(h,e,i,om,w,theta)
%Converts the classical Keplerian orbital elements into state vectors in
%the inertial frame of reference.

global mu

%Compute position and velocity in perifocal frame of reference:
R_per = h^2/mu * 1/(1+e*cosd(theta)) * [cosd(theta);sind(theta);0];
V_per = mu/h * [-sin(theta);e+cos(theta);0];
%Calculate transformation matrix to convert from perifocal to ECI:
QXx = [cos(w),sin(w),0;-sin(w),cos(w),0;0,0,1]*[1,0,0;0, cos(i),sin(i); 0,-sin(i),cos(i)]...
    *[cos(om),sin(om),0; -sin(om),cos(om),0; 0,0,1];
%Apply transformation to R and V to convert from perifocal to ECI:
R = QXx*R_per;
V = QXx*V_per;
end

