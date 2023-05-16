clc;clear;close all; 
% format short g;
[DATA,TEXT]=xlsread('EURUSD_DAX_SPX_final.xlsx');
% TIME_all=datetime(TEXT(2:end,1));
N=length(DATA(:,1));
TIME_num=zeros(N,6);
TIME=cell(N,1);
TIME_all=TEXT(2:end,1);
for i=1:N
    TIME{i,1}=split(TIME_all{i,1},{'/',' ',':'});
    for j=1:6
        TT=TIME{i,1}(j);
        TIME_num(i,j)=str2num(TT{1});
    end
end

%% Daily
L=1;
% D=1;
NL=1;
for i=1:N-1
Day{L,1}(NL,1:3)=DATA(i,1:3);
    if TIME_num(i+1,3)~=TIME_num(i,3)
        L=L+1;
        NL=1;
    else 
        NL=NL+1;
    end
end
Day{L,1}(NL,1:3)=DATA(N,1:3);
for i=1:length(Day)
    % Day_r(i,1:3)=mean(Day{i,1}(:,1:3))/100/252;
    Day_r_end(i,1:3)=Day{i,1}(end,1:3);
end
r_day_percent=mean(Day_r_end/100/252);
figure;
plot(Day_r_end(:,1),'-.r');
title('Daily(DAX)');
xlabel('day');
ylabel('Index');
figure;
plot(Day_r_end(:,2),'-.k');
title('Daily(SP500)');
xlabel('day');
ylabel('Index');
figure;
plot(Day_r_end(:,3),'-.b');
% legend('DAX','SP500','EUR-to-USD');
title('Daily(EUR-to-USD)');
xlabel('day');
ylabel('Index');
% xlswrite('Daily.xlsx',Day_r_end);
%% Weekly
N_week=1+2965/5;
Week_r_end(1,1:3)=Day_r_end(1,1:3);
for i=2:N_week
    Week_r_end(i,1:3)=Day_r_end(1+(i-1)*5,1:3);
end
r_week_percent=mean(Week_r_end/100/252);
% figure;
% plot(Week_r_end(:,1),'-.r');
% title('Weekly(DAX)');
% xlabel('week');
% ylabel('');
% figure;
% plot(Week_r_end(:,2),'-.k');
% title('Weekly(SP500)');
% xlabel('week');
% ylabel('');
% figure;
% plot(Week_r_end(:,3),'-.b');
% % legend('DAX','SP500','EUR-to-USD');
% title('Weekly(EUR-to-USD)');
% xlabel('week');
% ylabel('');
% xlswrite('Weekly.xlsx',Week_r_end);
%% hourly
L=1;
% D=1;
NL=1;
for i=1:N-1
    Hour{L,1}(NL,1:3)=DATA(i,1:3);
    if TIME_num(i+1,4)~=TIME_num(i,4)
        L=L+1;
        NL=1;
    else 
        NL=NL+1;
    end
end
Hour{L,1}(NL,1:3)=DATA(N,1:3);
for i=1:length(Hour)
Hour_r_end(i,1:3)=Hour{i,1}(end,1:3);
end
r_hour_percent=mean(Hour_r_end/100/252);
% figure;
% plot(Hour_r_end(:,1),'-.r');
% title('hourly(DAX)');
% xlabel('hour');
% ylabel('');
% figure;
% plot(Hour_r_end(:,2),'-.k');
% title('hourly(SP500)');
% xlabel('hour');
% ylabel('');
% figure;
% plot(Hour_r_end(:,3),'-.b');
% % legend('DAX','SP500','EUR-to-USD');
% title('hourly(EUR-to-USD)');
% xlabel('hour');
% ylabel('');
% xlswrite('hourly.xlsx',Hour_r_end);
%% 
% Returns simulation, estimation, std error computation
format long

% Settings
N_ret = 593; % number of return

% optimization options
options = optimoptions('fminunc','MaxFunEvals',15000,'MaxIter',3000,...
'TolFun',1e-6,'TolX',1e-6);

%% Return
% [row,column]=size(Day_r_end);
% [row,column]=size(Week_r_end);
[row,column]=size(Hour_r_end);  

% for i = 1:row-1
    for j =1:column
        
        r_annual_percent=0.45067;
        r=(r_annual_percent/100)/252;

        log_real= log(Hour_r_end(:,j));
        log_ret_real = (log_real(2:end) - log_real(1:end-1))';
        
        h_init = var(log_ret_real) ;
        
        x0 = [0,0,3,5,1];

        [x_max,f_val,flag,output] = fminunc(@(x) loss(x,r,h_init,log_ret_real),x0,options);
%         [x_max,f_val,flag,output] = fminunc(@(x) loss(x,r,h_init,ret_sim),x0,options);
        % output.message

        %std error
        gradient = numfirstderiv(@(x) hn_lik_vec(inv_rescale(x),r,h_init,log_ret_real),x_max,0.00001);
%         gradient = numfirstderiv(@(x) hn_lik_vec(inv_rescale(x),r,h_init,ret_sim),x_max,0.00001);
        std_error = sqrt(diag(inv(gradient'*gradient)));
        estimate_ret = inv_rescale(x_max);
        std_error_ret = inv_rescale(std_error);

        format short
%         return_real=table(["lambda", "omega","alpha","beta","gamma"]',estimate_ret',std_error_ret','VariableNames',["Para", "Est_ret", "StdErr_ret"]);
        return_real{1,j}=table(["lambda", "omega","alpha","beta","gamma"]',estimate_ret',std_error_ret','VariableNames',["Para", "Est_ret", "StdErr_ret"]);

        % est_ret{1,1}=table(["lambda", "omega","alpha","beta","gamma"]',estimate_ret',std_error_ret','VariableNames',["Para", "Est_ret", "StdErr_ret"]);
        format long

        % save('Output.mat','est_ret');
        % BBB=load('Output.mat');
    end
% end

% return lik vec only
function [error,z,h] = hn_lik_vec(theta,r,h_init,ret_frequency)
% This function receives parameters theta = (lambda,omega,alpha,beta,
% gamma), r, h_init and calculates z_1,...,z_N and returns 
% vector of log likelihood, implied residual and h

N = length(ret_frequency);
lambda = theta(1);
omega = theta(2);
alpha = theta(3);
beta = theta(4);
gamma = theta(5);

h = zeros(1,N+1); 
h(1) = h_init;
z = zeros(1,N);

% calculate h and z implied by stock prices
for i = 1:N
z(i) = (ret_frequency(i) - r - lambda * h(i)) / sqrt(h(i));
h(i+1) = omega + beta * h(i) + alpha * (z(i) - gamma * sqrt(h(i)))^2;
end

h = h(1:N);
% vectorized error
error = -.5 * (log(2*pi*h(1:N)) + z.^2);
% error = -.5 * (log(h(1:N)) + z.^2);
end

% scaling
function theta = inv_rescale(x)
% apply the scaling in order to make optimization efficient

lambda = x(1) * 10^2;
omega = x(2) * 10^(-6);
alpha = x(3) * 10^(-6);
beta = x(4) * 10^(-1);
gamma = x(5) * 10^(2);
theta = [lambda, omega, alpha, beta, gamma];
end

% optimize function
function error = loss(x,r,h_init,ret)
vec_ret = hn_lik_vec(inv_rescale(x),r,h_init,ret);
error = -1 * sum(vec_ret(1:end));
end

% first derivative via finite difference
function [gmatrix]=numfirstderiv(func_name, theta, eps)
% This function calculates the numerical first derivative of a function
% The function should return a Nx1 matrix and the parameter vector should
% be kx1. The returned value is a Nxk matrix

for i=1:length(theta)
theta_forward=theta;
theta_backward=theta;
% delta_theta=max(abs(theta(i)*eps),eps)
delta_theta=abs(theta(i)*eps);
theta_forward(i)=theta_forward(i)+delta_theta;
theta_backward(i)=theta_backward(i)-delta_theta;
gmatrix(:,i)=.5*(func_name(theta_forward)-func_name(theta_backward))./delta_theta;
end

if ~isreal(gmatrix)
disp('WARNING!')
disp('Complex valued gradient matrix.')
end

end