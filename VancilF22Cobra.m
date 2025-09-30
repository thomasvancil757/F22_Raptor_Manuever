
clear; clc;
syms t real
syms t0 t1 t2 t3 t4 real
syms x0 y0 z0 real
syms v1 v2 v3 v4 v5 real
syms chi0 chi1 chi2 chi3 chi4 real
syms gamma0 gamma1 gamma2 gamma3 gamma4 real

Vx1 = v1*cos(chi0)*cos(gamma0);  Vy1 = v1*sin(chi0)*cos(gamma0);  Vz1 = v1*sin(gamma0);
Vx2 = v2*cos(chi1)*cos(gamma1);  Vy2 = v2*sin(chi1)*cos(gamma1);  Vz2 = v2*sin(gamma1);
Vx3 = v3*cos(chi2)*cos(gamma2);  Vy3 = v3*sin(chi2)*cos(gamma2);  Vz3 = v3*sin(gamma2);
Vx4 = v4*cos(chi3)*cos(gamma3);  Vy4 = v4*sin(chi3)*cos(gamma3);  Vz4 = v4*sin(gamma3);

x_0=x0; y_0=y0; z_0=z0;
x_1 = x_0 + Vx1*(t1 - t0);  y_1 = y_0 + Vy1*(t1 - t0);  z_1 = z_0 + Vz1*(t1 - t0);
x_2 = x_1 + Vx2*(t2 - t1);  y_2 = y_1 + Vy2*(t2 - t1);  z_2 = z_1 + Vz2*(t2 - t1);
x_3 = x_2 + Vx3*(t3 - t2);  y_3 = y_2 + Vy3*(t3 - t2);  z_3 = z_2 + Vz3*(t3 - t2);
x_4 = x_3 + Vx4*(t4 - t3);  y_4 = y_3 + Vy4*(t4 - t3);  z_4 = z_3 + Vz4*(t4 - t3);

x_s1 = x_0 + Vx1*(t - t0);  y_s1 = y_0 + Vy1*(t - t0);  z_s1 = z_0 + Vz1*(t - t0);
x_s2 = x_1 + Vx2*(t - t1);  y_s2 = y_1 + Vy2*(t - t1);  z_s2 = z_1 + Vz2*(t - t1);
x_s3 = x_2 + Vx3*(t - t2);  y_s3 = y_2 + Vy3*(t - t2);  z_s3 = z_2 + Vz3*(t - t2);
x_s4 = x_3 + Vx4*(t - t3);  y_s4 = y_3 + Vy4*(t - t3);  z_s4 = z_3 + Vz4*(t - t3);

NaNsym = sym(NaN);
x_pw = piecewise(t<t0,NaNsym,(t>=t0)&(t<=t1),x_s1,(t>t1)&(t<=t2),x_s2,(t>t2)&(t<=t3),x_s3,(t>t3)&(t<=t4),x_s4,NaNsym);
y_pw = piecewise(t<t0,NaNsym,(t>=t0)&(t<=t1),y_s1,(t>t1)&(t<=t2),y_s2,(t>t2)&(t<=t3),y_s3,(t>t3)&(t<=t4),y_s4,NaNsym);
z_pw = piecewise(t<t0,NaNsym,(t>=t0)&(t<=t1),z_s1,(t>t1)&(t<=t2),z_s2,(t>t2)&(t<=t3),z_s3,(t>t3)&(t<=t4),z_s4,NaNsym);
r_pw = [x_pw; y_pw; z_pw];

times = [0 10 15 20 30];              
v0 = 10;                              

chis_seg = deg2rad([0 0 22.5 45 45]);
gammas_seg = deg2rad([0 0 75   0  0 ]);

x0n=0; y0n=0; z0n=0;

vars = {t0,t1,t2,t3,t4, x0,y0,z0, ...
        v1,chi0,gamma0, v2,chi1,gamma1, v3,chi2,gamma2, ...
        v4,chi3,gamma3, v5,chi4,gamma4};

vals = {times(1),times(2),times(3),times(4),times(5), ...
        x0n,y0n,z0n, ...
        v0,chis_seg(1),gammas_seg(1), ...
        v0,chis_seg(2),gammas_seg(2), ...
        v0,chis_seg(3),gammas_seg(3), ...
        v0,chis_seg(4),gammas_seg(4), ... 
        v0,chis_seg(5),gammas_seg(5)};

r_num = simplify(subs(r_pw, vars, vals));
matlabFunction(r_num(1), r_num(2), r_num(3), ...
    'Vars', t, ...
    'File', 'r_fun.m', ...
    'Outputs', {'x','y','z'});

chi_nodes_deg = [0 0 22.5 45 45 ];
gamma_nodes_deg = [0 0 75  0  0  ];

[times_sorted, idx] = sort(times(:));
chi_nodes_deg_sorted = chi_nodes_deg(idx);
gamma_nodes_deg_sorted = gamma_nodes_deg(idx);

F_chi_deg = griddedInterpolant(times_sorted, chi_nodes_deg_sorted,   'linear','linear');
F_gamma_deg = griddedInterpolant(times_sorted, gamma_nodes_deg_sorted, 'linear','linear');
chi_rad = @(t) deg2rad(F_chi_deg(t));
gamma_rad = @(t) deg2rad(F_gamma_deg(t));

vx = @(t) v0 .* cos(chi_rad(t)) .* cos(gamma_rad(t));
vy = @(t) v0 .* sin(chi_rad(t)) .* cos(gamma_rad(t));
vz = @(t) v0 .* sin(gamma_rad(t));

tt = times(:);           
vx_t = vx(tt);  vy_t = vy(tt);  vz_t = vz(tt);

tt = linspace(times(1), times(end), 5000);
x_ramp = x0n + cumtrapz(tt, vx(tt));
y_ramp = y0n + cumtrapz(tt, vy(tt));
z_ramp = z0n + cumtrapz(tt, vz(tt));

Vx = gradient(x_ramp, tt);
Vy = gradient(y_ramp, tt);
Vz = gradient(z_ramp, tt);

figure('Name','V0, Heading Angle, Flight Path Angle'); 
tiledlayout(3,1,'TileSpacing','compact');
nexttile; plot(tt, chi_rad(tt), 'LineWidth', 2); grid on; ylabel('Chi'); xlabel('Time');title('Chi(t)');
nexttile; plot(tt, gamma_rad(tt), 'LineWidth', 2); grid on; ylabel('Gamma'); xlabel('Time'); title('Gamma(t)');
nexttile; plot(times, v0*ones(size(times)), 'LineWidth', 2); grid on; ylabel('V0'); xlabel('Time'); title('Constant V0');

figure('Name','Velocity from Position (ramp)'); 
tiledlayout(3,1,'TileSpacing','compact');
nexttile; plot(tt, Vx, 'LineWidth', 2); grid on; ylabel('V_x'); xlabel('Time'); title('V_x(t)');
nexttile; plot(tt, Vy, 'LineWidth', 2); grid on; ylabel('V_y'); xlabel('Time'); title('V_y(t)');
nexttile; plot(tt, Vz, 'LineWidth', 2); grid on; ylabel('V_z'); xlabel('Time'); title('V_z(t)');

figure('Name','Position vs Time'); 
tiledlayout(3,1,'TileSpacing','compact');
nexttile; plot(tt, x_ramp, 'LineWidth', 2); grid on; ylabel('X(t)'); xlabel('Time'); title('x(t)');
nexttile; plot(tt, y_ramp, 'LineWidth', 2); grid on; ylabel('Y(t)'); xlabel('Time'); title('y(t)');
nexttile; plot(tt, z_ramp, 'LineWidth', 2); grid on; ylabel('Z(t)'); xlabel('Time'); title('z(t)');

figure('Name', 'Y vs X, Z vs X');
tiledlayout(2,1,'TileSpacing','Compact');
nexttile; plot(z_ramp, x_ramp,'LineWidth', 2); grid on; ylabel('Z(t)'); xlabel('X(t)'); title('Z(t) vs X(t)');
nexttile; plot(y_ramp, x_ramp,'LineWidth', 2); grid on; ylabel('Y(t)'); xlabel('X(t)'); title('Y(t) vs X(t)');

title('V_x(t)'); figure('Name','3D trajectory: chi & gamma');
plot3(x_ramp, y_ramp, z_ramp, 'LineWidth', 2);
grid on; box on;
axis([-50 250 -150 150 -50 250]); 
view(-115, 25);                   
axis vis3d; 
xlabel('X'); ylabel('Y'); zlabel('Z');

title('Trajectory from V0=10, chi(t), gamma(t)');
