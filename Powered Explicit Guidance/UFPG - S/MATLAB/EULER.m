function S2 = EULER(t0,x0,y0,z0,vx0,vy0,vz0,ax0,ay0,az0,T,dt)
%Uses Euler's method to integrate the system of ODE's, for a given
%step size dt and integration period T. Returns a matrix of t, r, a, p, q

%Pre-allocate the output arrays:
t = zeros(1,length(T/dt));
x = zeros(1,length(T/dt));
y = zeros(1,length(T/dt));
z = zeros(1,length(T/dt));
vx = zeros(1,length(T/dt));
vy = zeros(1,length(T/dt));
vz = zeros(1,length(T/dt));
ax = zeros(1,length(T/dt));
ay = zeros(1,length(T/dt));
az = zeros(1,length(T/dt));
Vars2 = zeros(length(T/dt),7);
UPFG_Array = zeros(length(T/dt),7);

%Initial conditions:
t(1) = t0;
x(1) = x0;
y(1) = y0;
z(1) = z0;
vx(1) = vx0;
vy(1) = vy0;
vz(1) = vz0;
ax(1) = ax0;
ay(1) = ay0;
az(1) = az0;
Vars2(1,:) = 0;
n = 1;

Rd = 0;
Rbias = 0;
Rgrav = 0;
Vd = 0;
Vgo = 0;
Tgo = 0;
spre = 0;
csetc = 0; csexc = 0; dtc = 0; xc = 0;
UPFG_vars = 0;

while t < T
    
    t(n+1) = t(n) + dt;                                                     %The time [s]

    %Call the equations of motion function once in order to capture the
    %debug array, and secondly to multiply by dt to add to the current
    %state variables to propagate forwards:
    
    Y = [x(n),y(n),z(n),vx(n),vy(n),vz(n),ax(n),ay(n),az(n)];
    State = EqOfMotion(t(n),Y, UPFG_vars);
    
    St = cell2mat(State(1:2));
    Vars = cell2mat(State(3));
    UPFG_vars = State(4);
    
    
    dx = St(1);
    dy = St(2);
    dz = St(3);
    dvx = St(4);
    dvy = St(5);
    dvz = St(6);
    Vars2(n,:) = Vars;
    
    x(n+1) = x(n) + dx*dt;
    y(n+1) = y(n) + dy*dt;
    z(n+1) = z(n) + dz*dt;
    vx(n+1) = vx(n) + dvx*dt;
    vy(n+1) = vy(n) + dvy*dt;
    vz(n+1) = vz(n) + dvz*dt;
    ax(n+1) = dvx;
    ay(n+1) = dvy;
    az(n+1) = dvz;

    n = n+1;                                                                %Update n for next solution
    
end

S2 = {t;x;y;z;vx;vy;vz;ax;ay;az;Vars2};

end
