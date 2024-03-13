clc
close all
clear all
%%
theta1 = 0;
theta2 = pi;

theta = theta2 - theta1;
q_max = theta;
v_max = 11.52;
a_max = 115;

t0 = 0;
t1 = v_max / a_max;
k1 = a_max/t1;
t2 = t1 + t1;
t3 = theta/v_max;
% t4 = t3 + t1;
% t5 = t4 + t1;
%%

%%
% while t3 < t2
%     k1 = a_max/t1;
%     t3 = t2;
%     %t1 = t2/2;
%     %a_max = k1*t1;
%     %v_max = t1*a_max;
%     %v_max = sqrt(theta*a_max/2);
%     v_max = theta/t3;
%     %t1 = t2/2;
% %     t1 = v_max/a_max;
% %     t2 = t1 + t1;
% end

% t3 = t2+5/10000;
t4 = t3 + t1;
t5 = t4 + t1;
%%
n = 5000;
total_step = round(t5*n);

% a(t) = k1*t 
% a(t1) = a_max => v(t1) = 1/2 v_max => q(t1) = 
% a(t2) = 0 => v(t2) = v_max
% v(t) = (1/2)*k1*(t^2) + k2
% v(t2) = v_max => v(t2) = v_max
% q(t) = (1/6)*k1*(t^3) + (1/2)*k2*(t^2) + k3
% q(t5) = q_max => k3

%k2 = v_max - (1/2)*k1*(t2^2);
%k3 = q_max - (1/6)*k1*(t5^3) + (1/2)*k2*(t5^2);


% a(1:total_step) = 0;
% t(1:total_step) = 0;


for i = 1:total_step
    
    if i <= t1*n
        t(i) = i/n;
        a(i) = k1*t(i);  
        v(i) = (1/2)*k1*(t(i)^2);
        q(i) = (1/6)*k1*(t(i)^3);
        i1 = i;
    elseif i <= t2*n
        t(i) = i/n;
        a(i) = -k1*(t(i)-t1) + a_max;
        %v(i) = (-1/2)*k1*(t(i)^2) + (k1*t1 + a_max)*(t(i)) + k2;
        v(i) = -((1/2)*k1*((t(i)-t1)^2)) + a_max*(t(i)-t1) + v(i1);
        %q(i) = (-1/6)*k1*(t(i)^3) + (1/2)*(k1*t1 + a_max)*(t(i)^2) - v_max*t(i); 
        q(i) = -(1/6)*k1*((t(i)-t1)^3) + a_max*(1/2)*((t(i)-t1)^2) + v(i1)*(t(i)-t1) + q(i1);
        i2 = i;
        i3 = i2;
%         if t3 == t2
%             i3 = i2;
%         end
    elseif i <= t3*n
        t(i) = i/n;
        a(i) = 0;
        v(i) = v_max;
        q(i) = v_max*(t(i)-t2)+ q(i2);
        i3 = i;
    elseif i <= t4*n
        t(i) = i/n;
        a(i) = -k1*(t(i)-t3);
        v(i) = -((1/2)*k1*((t(i)-t3)^2)) + v(i3);
        q(i) = -((1/6)*k1*((t(i)-t3)^3)) + v(i3)*(t(i)-t3) + q(i3);
        i4 = i;
    elseif i<= t5*n
        t(i) = i/n;
        a(i) = k1*(t(i)-t4) - a_max;
        v(i) = ((1/2)*k1*((t(i)-t4)^2) - a_max*(t(i) - t4)) + v(i4);
        q(i) = ((1/6)*k1*((t(i)-t4)^3) - (1/2)*a_max*(t(i)-t4)^2) + v(i4)*(t(i)-t4) + q(i4);
        i5 = i;
    end
%     plot(i*tf,a(i));
%     hold on
end
subplot(3,1,1);
plot(t,a);
subplot(3,1,2);
plot(t,v);
subplot(3,1,3);
plot(t,q);