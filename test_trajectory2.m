clc
clear all

axis = [0       pi      115     660*pi/180; 
        -pi/6   pi      115     11.5; 
        0.5     1.4     112     11.2; 
        0       pi/6    262     1500*pi/180]; % qmax qmin a_max v_min
axis = [
   -0.5236    0.5236  115.2000   11.5192;
   -1.7141    1.7141  115.2000   11.5192;
    1.4000    0.2000  112.0000   11.2000;
         0    3.1416  261.8000   26.1799];
qmin = [axis(1,1)  axis(2,1) axis(3,1) axis(4,1)];
qmax = [axis(1,2)  axis(2,2) axis(3,2) axis(4,2)];
amax = [axis(1,3)  axis(2,3) axis(3,3) axis(4,3)];
vmax = [axis(1,4)  axis(2,4) axis(3,4) axis(4,4)];


% f1 = 1;
% for j = 1:4
%     m = create_trajectory(qmin(j),qmax(j),amax(j),vmax(j));
%     for i = 1:3
%         subplot(4,3,f1);
%         plot(m(1,:),m(i+1,:));
%         f1 = f1+1;
%     end
% end

m1 = create_trajectory(qmin(1),qmax(1),amax(1),vmax(1));
t1 = m1(1,end);
m2 = create_trajectory(qmin(2),qmax(2),amax(2),vmax(2));
t2 = m2(1,end);
m3 = create_trajectory(qmin(3),qmax(3),amax(3),vmax(3));
t3 = m3(1,end);
m4 = create_trajectory(qmin(4),qmax(4),amax(4),vmax(4));
t4 = m4(1,end);

n = 500; i = 2;
f1 = 1;
while (i < t1*n || i < t2*n || i < t3*n || i < t4*n)
    if i <= t1*n
        subplot(4,3,1);
        plot(m1(1,i-1:i),m1(2,i-1:i));
        hold on
        subplot(4,3,2);
        plot(m1(1,i-1:i),m1(3,i-1:i));
        hold on
        subplot(4,3,3);
        plot(m1(1,i-1:i),m1(4,i-1:i)); 
        hold on
        end1 = i;

    else
        end1 = floor(t1*n);
    end

    if i <= t2*n
        subplot(4,3,4);
        plot(m2(1,i-1:i),m2(2,i-1:i));
        hold on
        subplot(4,3,5);
        plot(m2(1,i-1:i),m2(3,i-1:i));
        hold on
        subplot(4,3,6);
        plot(m2(1,i-1:i),m2(4,i-1:i));
        hold on
        end2 = i;
    else
        end2 = floor(t2*n);
    end

    if i <= t3*n
        subplot(4,3,7);
        plot(m3(1,i-1:i),m3(2,i-1:i));
        hold on
        subplot(4,3,8);
        plot(m3(1,i-1:i),m3(3,i-1:i));
        hold on
        subplot(4,3,9);
        plot(m3(1,i-1:i),m3(4,i-1:i)); 
        hold on
        end3 = i;
    else
        end3 = floor(t3*n);
    end

    if i <= t4*n
        subplot(4,3,10);
        plot(m4(1,i-1:i),m4(2,i-1:i));
        hold on
        subplot(4,3,11);
        plot(m4(1,i-1:i),m4(3,i-1:i));
        hold on
        subplot(4,3,12);
        plot(m4(1,i-1:i),m4(4,i-1:i));
        hold on
        end4 = i;
    else
        end4 = floor(t4*n);
    end
    %main(app,m1(4,end1),m2(4,end2),m3(4,end3),m4(4,end4));
    i = i + 1
    m1(4,end1) + qmin(1)
    m2(4,end2) + qmin(2)
    m3(4,end3) + qmin(3)
    m4(4,end4) + qmin(4)
    pause(1/1000);   
end