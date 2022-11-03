function Vrot = Rodrigues(V,rot,theta)

% RODRIGUES ROTATION EQUATION
%
% Rotates an initial vector V through an angle theta about a specified axis
% (rot).

Vrot = V*cosd(theta) + cross(rot,V)*sind(theta) + rot*dot(rot,V)*(1-cosd(theta));


end

