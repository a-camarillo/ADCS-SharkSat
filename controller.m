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
jSat = diag([J_x, J_y, J_z]);

% Calculate Reaction Wheel Inertia Tensor

% Attitude Dynamics
% Source: Fundamentals of Spacecraft Attitude Determination and Control - Crassidis and Markley

% Euler's rotational equation with Reaction Wheels
% omegaDot: Time derivate of angular velocity
% J: Body Moment of Inertia
% T: Net Torque acting on body
% omega: angular velocity in the body frame
%
% omegaDot = inv(jSat)*[T - omega x (jSat*omega + J_wheel^(-1)*T_wheel)]

% Define our symbolic variables
syms tWheelX tWheelY tWheelZ omegaX omegaY omegaZ hWheelX hWheelY hWheelZ q1 q2 q3 q4

% spacecraft angular velocity in inertial frame
omegaSat = [omegaX; omegaY; omegaZ];

% quaternion vector

q = [q1; q2; q3; q4];

% applied torque from reaction wheels
tWheel = [tWheelX; tWheelY; tWheelZ];

% angular momentum of reaction wheels in body frame
hWheel = [hWheelX; hWheelY; hWheelZ];

thetaQ = [q4 -q3 q2; q3 q4 -q1; -q2 q1 q4; -q1 -q2 -q3];

% Euler's rotational equation of motion with a reaction wheel
%omegaDot = inv(jSat)*(-cross(omegaSat,(jSat*omega_sat + hWheel)) - T);
%omegaDot_no_torque = inv(jSat)*(-cross(omegaSat,(jSat*omega_sat + hWheel)));

% Finding control inputs(T) for f(omega, H, T) = 0
% tWheelX = omegaZ(h_wheel_x + (41/960)omegaY) - omegaY(h_wheel_z + (41/4800)omegaZ)
% tWheelY = omegaX(h_wheel_z + (41/4800)omegaZ) - omegaZ(h_wheel_x + (41/960)omegaX)
% tWheelZ = omegaY(h_wheel_x + (41/960)omegaX) - omegaX(h_wheel_y + (41/960)omegaY)

%omegaDotJacobian = jacobian(omegaDot, [omegaSat; hWheel]);
%inputJacobian = jacobian(omegaDot, T);
%
%u1 = omegaZ*(h_wheel_x + (41/960)*omegaY) - omegaY*(h_wheel_z + (41/4800)*omegaZ);
%u2 = omegaX*(h_wheel_z + (41/4800)*omegaZ) - omegaZ*(h_wheel_x + (41/960)*omegaX);
%u3 = omegaY*(h_wheel_x + (41/960)*omegaX) - omegaX*(h_wheel_y + (41/960)*omegaY);
%u = [u1; u2; u3];

%omegaDotLinear = omegaDotJacobian*omegaSat + input_jacobian*u;


 
% From Section 7.2 in Fundamentals of Spacecraft Attitude Determination and Control - Crassidis and Markley
% Euler's rotational equation can be separated into two equations

% To get the open loop transfer function we look at the system without feedback
% omegaDot = inv(jSat)*(-cross(omegaSat,(jSat*omegaSat)) + tWheel);
% hDot = -cross(omegaSat,hWheel) - tWheel;

%{ 
    From Section 2.7 in Fundamentals of Spacecraft Attitude Determination and Control - Crassidis and Markley
    Notation for quaternions
    THETA(q) = [q4 -q3 q2; q3 q4 -q1; -q2 q1 q4; -q1 -q2 -q3]
    OMEGA(omega) = [0 w3 -w2 w1; -w3 0 w1 w2; w2 -w1 0 w3; w2 -w1 0 w3; -w1 -w2 -w3 0]
    then q_dot = (1/2)*THETA(q)*omega

    Finding attitude quaternion error
    With qDesired as a constant desired quaternion
    deltaQ = [ deltaQDesired(1:3); deltaQDesired(4) ]
    where
    deltaQ(1:3) = THETA(qDesired)' * q
    deltaQ(4) = q' * qDesired

    Then
    deltaQDot = (1/2) * OMEGA(omega) * deltaQ
    where
    deltaQDot(1:3) = (1/2) * [0 -deltaQ(3) deltaQ(2); deltaQ(3) 0 -deltaQ(1); -deltaQ(2) deltaQ(1) 0]*omega + (1/2) * deltaQ(4)*omega
    deltaQDot(4) = -(1/2) * deltaQ(1:3)' * omega
%}

qDesired = [0; 0; 0; 1];
thetaQDesired = [
            qDesired(4) -qDesired(3) qDesired(2);
            qDesired(3) qDesired(4) -qDesired(1);
            -qDesired(2) qDesired(1) qDesired(4);
            -qDesired(1) -qDesired(2) -qDesired(3);
            ];

% q.' is used to return the nonconjugate transpose
deltaQ = [thetaQDesired'*q; q.'*qDesired];
deltaQDot = sym(zeros(4,1));
deltaQDot(1:3) = (0.5)*([0 -deltaQ(3) deltaQ(2); deltaQ(3) 0 -deltaQ(1); -deltaQ(2) deltaQ(1) 0])*omegaSat + ((0.5)*deltaQ(4)*omegaSat);
deltaQDot(4) = -(1/2) * deltaQ(1:3).' * omegaSat;

% Source: Fundamentals of Spacecraft Attitude Determination and Control - Crassidis and Markley
% Section 7.2 pg. 290
% From here we have the feedback controller
% L = -kP*deltaQ(1:3) - kD*omegaSat
kP = 50;
kD = 300;
tWheel = -kP*deltaQ(1:3) - kD*omegaSat;

qDot = (1/2)*thetaQ*omegaSat;

omegaDot = inv(jSat)*(-cross(omegaSat,(jSat*omegaSat)) + tWheel);
hDot = -cross(omegaSat,hWheel) - tWheel;

% Finding the equilibrium points
%{
    For a function f(x,u) the equilibrium points will occur where the time
    derivatives of this function are equal to zero.
    So setting 
    omegaDot = 0 = inv(jSat)*(-cross(omegaSat,(jSat*omega_sat)) + T)
    h_dot = 0 = -cross(omegaSat,hWheel) - T;
    we can find the values of our states and input that correspond to an
    equilibrium point. 

    From setting wZ = 1, uZ=hZ=0
    It is found that
    wX = (4800/164)uY
    wY = -(4800/164)uX
    uX = hY
    uY = -hX
%}

spacecraftSystem = [omegaDot; deltaQDot];

systemJacobian = jacobian(spacecraftSystem, [omegaSat; q]);

systemMatrix = double(subs(systemJacobian, [omegaSat; q], [0.53; 0.53; 0.053; 0.6853; 0.6953; 0.1531; 0.1531]));
inputMatrix = [];
outputMatrix = eye(7);
disturbanceMatrix = [];

spacecraftStateSpace = ss(systemMatrix, inputMatrix, outputMatrix, disturbanceMatrix);

x0 = [0, 0, 0, 0.1, 0.1, 0.1 0.1];

figure;
states = {'\omega_x', '\omega_y', '\omega_z', '\delta_q_1', '\delta_q_2', '\delta_q_3', '\delta_q_4'};
[~, time, stateTrajectories]  = initial(spacecraftStateSpace, x0);
title('PD Controller Response for Small Initial Quaternion Error');
for i=1:length(states)
    subplot(length(states), 1, i);
    plot(time, stateTrajectories(:, i))
    ylabel(states{i});
    grid on;
end
sgtitle('PD Controller Response for Initial Small Quaternion Error')
xlabel('Time')
