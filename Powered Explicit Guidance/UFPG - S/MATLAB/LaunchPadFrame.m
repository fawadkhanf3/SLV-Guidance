function Out = LaunchPadFrame(X, Y, Z, r0,ph,th)
% Plot trajectory in the pad reference frame

% Latitude and longitude of the pad:



% ECI position vector:

R = [X;Y;Z];

% Compute matrix to transform from inertial to topocentric:

Q = [-sind(th), cosd(th), 0; 
    -sind(ph)*cosd(th), -sind(ph)*sind(th), cosd(ph);
    cosd(ph)*cosd(th), cosd(ph)*sind(th), sind(ph)];

% Position vector of vehicle relative to pad:

q = R - r0;
r_top = zeros(3, length(R));

for i=1:length(R)
   
    r_top(:,i) = Q*q(:,i);
    
end
%r_top = Q*q;

plot3(r_top(1,:)/1e3, r_top(2,:)/1e3, r_top(3,:)/1e3,'k');
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
grid on
%axis equal
%axis([0 200 0 200 0 200])
end

