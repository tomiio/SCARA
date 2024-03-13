x = 0; y =0; z = 1;

theta = 0;
theta_x = 0;
theta_y = 0;
line( [x x+cos(theta)], [y y+sin(theta)], [z z], 'LineWidth',5,'Color','red');
theta = theta + pi/2;
line( [x x+cos(theta)], [y y+sin(theta)], [z z], 'LineWidth',5,'Color','green');
line( [x x], [y y], [z z-1], 'LineWidth',5,'Color','blue');


x = 5; y =0; z = 1;
Lx = 1 - x;
Ly = 1 - y;
theta = 0;
line( [x x+cos(theta)], [y y+sin(theta)], [z z], 'LineWidth',5,'Color','yellow');
theta = theta + pi/2;
line( [x x+cos(theta)], [y y+sin(theta)], [z z], 'LineWidth',5,'Color','black');
line( [x x], [y y], [z z+   1], 'LineWidth',5,'Color','blue');

grid on;
view(3);