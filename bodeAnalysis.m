%{
SharkSat ADCS Gain and Phase Margin
Author: Anthony Camarillo
The purpose of this script is to create a linear model from the dynamic attitude equations and
perform controller analysis from this model.

Linearization can be done by taking the dynamic equation x_dot = f(x,u) and taking the jacobian about some
equilibrium point 
%}
clear
clc


% Calculate Spacecraft Inertia Tensor 
mass_s = 5.125;
height_s = 0.3;
depth_s = 0.1;
width_s = 0.1;
J_x = (1/12) * mass_s * (depth_s^2 + height_s ^ 2);
J_y = (1/12) * mass_s * (width_s^2 + height_s ^ 2);
J_z = (1/12) * mass_s * (width_s^2 + depth_s ^ 2);
J_sat = diag([J_x, J_y, J_z]);

% Calculate Reaction Wheel Inertia Tensor


% Attitude Dynamics
% Source: Fundamentals of Spacecraft Attitude Determination and Control - Crassidis and Markley

% Euler's rotational equation with Reaction Wheels
% omega_dot: Time derivate of angular velocity
% J: Body Moment of Inertia
% T: Net Torque acting on body
% omega: angular velocity in the body frame
%
% omega_dot = inv(J_sat)*[T - omega x (J_sat*omega + J_wheel^(-1)*T_wheel)]

% Define our symbolic variables
syms T_wheel_x T_wheel_y T_wheel_z omegaX omegaY omegaZ h_wheel_x h_wheel_y h_wheel_z q1 q2 q3 q4

% spacecraft angular velocity in inertial frame
omega_sat = [omegaX; omegaY; omegaZ];

% quaternion vector

q = [q1; q2; q3; q4];

% applied torque from reaction wheels
T = [T_wheel_x; T_wheel_y; T_wheel_z];

% angular momentum of reaction wheels in body frame
H_wheel = [h_wheel_x; h_wheel_y; h_wheel_z];

THETA_q = [q4 -q3 q2; q3 q4 -q1; -q2 q1 q4; -q1 -q2 -q3];

% Euler's rotational equation of motion with a reaction wheel
%omega_dot = inv(J_sat)*(-cross(omega_sat,(J_sat*omega_sat + H_wheel)) - T);
%omega_dot_no_torque = inv(J_sat)*(-cross(omega_sat,(J_sat*omega_sat + H_wheel)));

% Finding control inputs(T) for f(omega, H, T) = 0
% T_wheel_x = omegaZ(h_wheel_x + (41/960)omegaY) - omegaY(h_wheel_z + (41/4800)omegaZ)
% T_wheel_y = omegaX(h_wheel_z + (41/4800)omegaZ) - omegaZ(h_wheel_x + (41/960)omegaX)
% T_wheel_z = omegaY(h_wheel_x + (41/960)omegaX) - omegaX(h_wheel_y + (41/960)omegaY)

%omega_dot_jacobian = jacobian(omega_dot, [omega_sat; H_wheel]);
%input_jacobian = jacobian(omega_dot, T);
%
%u1 = omegaZ*(h_wheel_x + (41/960)*omegaY) - omegaY*(h_wheel_z + (41/4800)*omegaZ);
%u2 = omegaX*(h_wheel_z + (41/4800)*omegaZ) - omegaZ*(h_wheel_x + (41/960)*omegaX);
%u3 = omegaY*(h_wheel_x + (41/960)*omegaX) - omegaX*(h_wheel_y + (41/960)*omegaY);
%u = [u1; u2; u3];

%omega_dot_linear = omega_dot_jacobian*omega_sat + input_jacobian*u;


 
% From Section 7.2 in Fundamentals of Spacecraft Attitude Determination and Control - Crassidis and Markley
% Euler's rotational equation can be separated into two equations

% To get the open loop transfer function we look at the system without feedback
omega_dot = inv(J_sat)*(-cross(omega_sat,(J_sat*omega_sat)) + T);
h_dot = -cross(omega_sat,H_wheel) - T;

% From Section 2.7 in Fundamentals of Spacecraft Attitude Determination and Control - Crassidis and Markley
% Notation for quaternions
% THETA(q) = [q4 -q3 q2; q3 q4 -q1; -q2 q1 q4; -q1 -q2 -q3]
% then q_dot = (1/2)*THETA(q)*omega

q_dot = (1/2)*THETA_q*omega_sat;

% Finding the equilibrium points
%{
    For a function f(x,u) the equilibrium points will occur where the time
    derivatives of this function are equal to zero.
    So setting 
    omega_dot = 0 = inv(J_sat)*(-cross(omega_sat,(J_sat*omega_sat)) + T)
    h_dot = 0 = -cross(omega_sat,H_wheel) - T;
    we can find the values of our states and input that correspond to an
    equilibrium point. 
%}

spacecraft_system = [omega_dot; h_dot; q_dot];

system_jacobian = jacobian(spacecraft_system, [omega_sat; H_wheel; q]);
input_jacobian = jacobian(spacecraft_system, T);

% Linearize the point around [wX, wY, wZ, hX, hY, hZ, q1, q2, q3, q4] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
% double() around the A and B matrices to convert them to numeric arrays or else system
% functions will throw an error
A = double(subs(system_jacobian,[omega_sat; H_wheel; q], [0; 0; 0; 0; 0; 0; 0; 0; 0; 1]));
A_mod = [A(1:3,1:6);A(7:9,1:6)];
B = double(input_jacobian);
B_mod = [B(1:3,:); B(7:9,:)];
C = [zeros(6,10); zeros(4,6), eye(4,4)];
C_mod = eye(6);
D = zeros(10,3);
D_mod = zeros(6,3);

poles = [-1-1i, -1+1i, -2-1i, -2+1i, -3-1i, -3+1i];

K = place(A_mod, B_mod, poles);
L = place(A_mod', C_mod', 5*real(poles));

Q = [
     1, 0, 0, 0 ,0, 0;
     0, 1, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0;
     0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 1;
];

R = [
     1, 0, 0;
     0, 1, 0;
     0, 0, 1;
];


lqrGain = lqr(A_mod, B_mod, Q, R);

A_feedback_observer = [(A_mod - B_mod*K), B_mod*K; zeros(6), (A_mod - L*C_mod)];
B_feedback_observer = [B_mod; zeros(6,3)];

systemObs = obsv(A_mod, C_mod);
systemCtrb = ctrb(A_mod,B_mod);

spacecraftLQR = ss(A_mod-B_mod*lqrGain, B_mod, C_mod, D_mod);

spacecraftSS = ss(A_feedback_observer,B_feedback_observer,eye(12),zeros(12,3));

% From control theory the transfer function G(s) is given by
% G(s) = C*((s*I-A)^-1)*B + D
s = tf('s');
spacecraftTF = spacecraftSS.C*inv(s*eye(12) - spacecraftSS.A)*spacecraftSS.B + spacecraftSS.D;

[sys_response, tOut] = initial(spacecraftSS, [0,0,0,0,0,0,1,1,1,1,1,1]);
[sysResponseLQR, tOutLQR] = initial(spacecraftLQR,[0,0,0,1,1,1]);


figure(1)
plot(tOutLQR, sysResponseLQR)
title('System Initial Response')


%states = ["{\omega}_x", "{\omega}_y", "{\omega}_z", "q_1", "q_2", "q_3"];
%
%for j=1:3
%    for i=1:6
%        switch j
%            case 1
%                figure(i)
%            case 2
%                figure(6+i)
%            case 3
%                figure(12+i)
%        end
%        nyquistplot(spacecraftTF(i,j))
%        title(['Output: ', char(states(i)), 'for Input: ', num2str(j)])
%    end
%end

%figure(2)
%bode(spacecraftTF)
%title('Bode Plot for Transfer Function Response')