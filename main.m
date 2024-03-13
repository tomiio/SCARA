%function main(theta1,theta2,d3,theta4)
clf
%% D-H parameters
% Link 1 D-H
a1 = 3; alpha1 = 0; d1 = 5; theta1 = -pi/3;
A01 = DH(a1,alpha1,d1,theta1);

% Link 2 D-H
a2 = 3; alpha2 = pi; d2 = 0; theta2 = pi/3;
A12 = DH(a2,alpha2,d2,theta2);

% Link 3 D-H
a3 = 0; alpha3 = 0; d3 = 1; 
theta3 = 0;
A23 = DH(a3,alpha3,d3,theta3);

% Link 4 D-H
a4 = 0; alpha4 = 0; d4 = 1; theta4 = pi/2;
A34 = DH(a4,alpha4,d4,theta4);

%%% Get position of joints
% Position of joint 1
link1_end = A01(1:3,4);

% Position of joint 2
A02 = A01*A12;
link2_end = A02(1:3,4);

% Positon of joint 3
A03 = A02*A23;
link3_end = A03(1:3,4);

% Position of joint 4
A04 = A03*A34;
link4_end = A04(1:3,4);

% End effector position and orientation
position_ee = A04(1:3,4);
orientation_ee = A04(1:3,1:3);

%% Inverse
X = -3;
Y = -4;
Z = 3;

% theta2 = acos((X^2 + Y^2 - a1^2 - a2^2)/(2*a1*a2));
% theta1 = atan(Y/X) - atan((a2*sin(theta2))/(a1 + a2*cos(theta2)));
% 
% d3 = Z + d4;
% theta4 = pi/2;

cos2 = (X^2+Y^2-a1*a1-a2*a2)/(2*a1*a2);

sin2 = sqrt(1-cos2*cos2);


theta1 = atan(Y/X)-atan(a2*sin2/(a1+a2*cos2));
theta2 = atan(sin2/cos2);

if (X < 0) 
    theta1 = pi + theta1;
end

d3 = Z+d4 ;

A01 = DH(a1,alpha1,d1,theta1);
A12 = DH(a2,alpha2,d2,theta2);
A23 = DH(a3,alpha3,d3,theta3);
A34 = DH(a4,alpha4,d4,theta4);
% Position of joint 1
link1_end = A01(1:3,4);

% Position of joint 2
A02 = A01*A12;
link2_end = A02(1:3,4);

% Positon of joint 3
A03 = A02*A23;
link3_end = A03(1:3,4);

% Position of joint 4
A04 = A03*A34;
link4_end = A04(1:3,4);

%% Draw objects between joints

% Origin to start of link 1
create_cylinder(0,0,0,0.75,0.5,'black');
%line( [0 0], [0 0], [0 d1], 'LineWidth',5,'Color','yellow');
create_cylinder(0,0,0,0.5,d1,'r');
create_cylinder(0,0,link1_end(3)-0.2,0.6,0.5,'blue');


% start of link 1 to end of link 1
%line([0 link1_end(1)], [0 link1_end(2)], [d1 link1_end(3)], 'LineWidth',5,'Color','red');
create_box(0,0,link1_end(3),a1,0.4,0.2,0,0,theta1,0,0,'blue');
create_cylinder(link1_end(1),link1_end(2),link1_end(3)-0.3,0.6,0.6,'blue');
create_cylinder(link1_end(1),link1_end(2),link1_end(3)+0.3,0.5,0.6,'green');

% End of link 1 to end of link 2
line( [link1_end(1) link2_end(1)], [link1_end(2) link2_end(2)], [link1_end(3) link2_end(3)]+0.5, 'LineWidth',5,'Color','green');
create_box(link1_end(1),link1_end(2),link1_end(3)+0.5,a2,0.4,0.2,0,0,theta1+theta2,link1_end(1),link1_end(2),'green');
create_cylinder(link2_end(1),link2_end(2),link2_end(3)+0.2,0.5,0.6,'green');

% End of link 2 to end of link 3
line( [link2_end(1) link3_end(1)], [link2_end(2) link3_end(2)], [link2_end(3)+5-d3 link3_end(3)], 'LineWidth',5,'Color','black');
create_cylinder(link3_end(1),link3_end(2),link3_end(3),0.2,5,'black');

% End of link 3 to end of link 4
create_box(link4_end(1)-0.5,link4_end(2),link3_end(3),1,0.1,0.1,0,0,theta4,link4_end(1),link4_end(2),'black');

%% PLOT
xlabel('x');
ylabel('y');
zlabel('z');
grid on;
axis('equal');
xlim([-7 7]);
ylim([-7 7]);
zlim([0 10]);

view(3);

