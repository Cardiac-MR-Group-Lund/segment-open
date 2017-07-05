function Olocestimate = kalman(x,objectnoise)
%Kalman filter of a position vector assuming a velocity and acceleration
%model. Inspired by code by "StudentDave". 
% - x is input vector (i.e position).
% - objectnoise is noise estimate of position estimates (i.e how much noise
%   do we think it have been added).
%
%Written by Einar Heiberg

%Define variables
dt = 1; %Delta time

% Define update equations (Coefficent matrices): A physics based model for 
%where we expect the object to be [state transition (state + velocity)] + [input control (acceleration)]
A = [1 dt; 0 1] ; % state transition matrix:  expected flight of the object (state prediction)
C = [1 0]; % measurement matrix: the expected measurement given the predicted state (likelihood)

% define main variables
Ez = 1;% Ez convert the measurement noise (stdv) into covariance matrix
Ex = objectnoise^2 * [dt^4/4 dt^3/2; dt^3/2 dt^2]; % Ex convert the process noise (stdv) into covariance matrix
P = Ex; % estimate of initial Object position variance (covariance matrix)

%%% Do kalman filtering
%initize estimation variables
Olocmeas = x; %  object position estimate
Oestimate = [x(1); (x(2)-x(1))/dt]; %initized state--it has two components: [position; velocity] of the object
Olocestimate = zeros(1,length(x)); %  object position estimate

for t = 1:length(x)
    
  % Predict next state of the object with the last state and predicted motion.
    Oestimate = A * Oestimate; % + B * u;    
    
    %predict next covariance
    P = A * P * A' + Ex;
    
    % predicted object measurement covariance
    % Kalman Gain
    %K = P*C'*inv(C*P*C'+Ez);
    K = (P*C')/(C*P*C'+Ez); %faster and more efficient than expression above
    
    % Update the state estimate.
    Oestimate = Oestimate + K * (Olocmeas(t) - C * Oestimate);
    
    % update covariance estimation.
    P =  (eye(2)-K*C)*P;

    %Store for plotting
    Olocestimate(t) = Oestimate(1);  
end