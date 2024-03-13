clear all
clc
%%
axis = [0 pi 115 660*pi/180; -pi/6 pi/4 115 11.5; 0.5 1.4 112 11.2; 0 pi/6 262 1500*pi/180]; % qmax qmin a_max v_min
qmin = [axis(1,1)  axis(2,1) axis(3,1) axis(4,1)];
qmax = [axis(1,2)  axis(2,2) axis(3,2) axis(4,2)];
qm = qmax - qmin;
for i = 1:4
    if axis(i,3) < 2*axis(i,4)^2/qm(i)
        axis(i,3) = 2*axis(i,4)^2/qm(i);
    end
end
amax = [axis(1,3)  axis(2,3) axis(3,3) axis(4,3)];
vmax = [axis(1,4)  axis(2,4) axis(3,4) axis(4,4)];
t1(1:4) = 0;
t2(1:4) = 0;
t3(1:4) = 0;
t4(1:4) = 0;
t5(1:4) = 0;
k1(1:4) = 0;

t5max = 0;
imax = 0;
%%
for i = 1:4
    t3(i) = qm(i)/vmax(i);
    t1(i) = vmax(i)/amax(i);
    k1(i) = amax(i)/t1(i);
    t2(i) = t1(i)*2;
    t4(i) = t3(i) + t1(i);
    t5(i) = t4(i) + t1(i);
    if t5(i) > t5max
        t5max = t5(i);
        imax = i;
    end
end
%%
nmax = 5000;
total_step = round(t5max*nmax);
n(1:4) = 5000;
% for i = 1:4
%     if t5(i) ~= t5max
%         n(i) = t5max*nmax/t5(i);
%     else
%         n(i) = nmax;
%     end
% end
t(1:4,1:total_step) = 0;
q(1:4,1:total_step) = 0;
v(1:4,1:total_step) = 0;
a(1:4,1:total_step) = 0;
i1(4) = 0;
i2(4) = 0;
i3(4) = 0;
i4(4) = 0;

for j = 1:4   
    for i = 1:total_step

        if i <= t1(j)*n(j)
            t(j,i) = i/n(j);
            a(j,i) = k1(j)*t(j,i);  
            v(j,i) = (1/2)*k1(j)*(t(j,i)^2);
            q(j,i) = (1/6)*k1(j)*(t(j,i)^3);
            i1(j) = i;
        elseif i <= t2(j)*n(j)
            t(j,i) = i/n(j);
            a(j,i) = -k1(j)*(t(j,i)-t1(j)) + amax(j);
            v(j,i) = -((1/2)*k1(j)*((t(j,i)-t1(j))^2)) + amax(j)*(t(j,i)-t1(j)) + v(j,i1(j));
            q(j,i) = -(1/6)*k1(j)*((t(j,i)-t1(j))^3) + amax(j)*(1/2)*((t(j,i)-t1(j))^2) + v(j,i1(j))*(t(j,i)-t1(j)) + q(j,i1(j));
            i2(j) = i;
            i3(j) = i2(j);
        elseif i <= t3(j)*n(j)
            t(j,i) = i/n(j);
            a(j,i) = 0;
            v(j,i) = vmax(j);
            q(j,i) = vmax(j)*(t(j,i)-t2(j))+ q(j,i2(j));
            i3(j) = i;
        elseif i <= t4(j)*n(j)
            t(j,i) = i/n(j);
            a(j,i) = -k1(j)*(t(j,i)-t3(j));
            v(j,i) = -((1/2)*k1(j)*((t(j,i)-t3(j))^2)) + v(j,i3(j));
            q(j,i) = -((1/6)*k1(j)*((t(j,i)-t3(j))^3)) + v(j,i3(j))*(t(j,i)-t3(j)) + q(j,i3(j));
            i4(j) = i;
        elseif i<= t5(j)*n(j)
            t(j,i) = i/n(j);
            a(j,i) = k1(j)*(t(j,i)-t4(j)) - amax(j);
            v(j,i) = ((1/2)*k1(j)*((t(j,i)-t4(j))^2) - amax(j)*(t(j,i) - t4(j))) + v(j,i4(j));
            q(j,i) = ((1/6)*k1(j)*((t(j,i)-t4(j))^3) - (1/2)*amax(j)*(t(j,i)-t4(j))^2) + v(j,i4(j))*(t(j,i)-t4(j)) + q(j,i4(j));
        end
    end
end

j=1;
for i = 1:4
    subplot(4,3,j);
    plot(t(i,:),a(i,:));
    subplot(4,3,j+1);
    plot(t(i,:),v(i,:));
    subplot(4,3,j+2);
    plot(t(i,:),q(i,:));
    j = j+3;
end
