function create_box(x0,y0,z0,d,r,h,theta_x,theta_y,theta_z,rotx,roty,color)
%x0 = 1;y0 = 4; z0 = 2; xr = 1; yr = 1; zr = 1; 

% Define the initial coordinates of the box corners
x = [x0 x0+d x0+d x0 x0 x0+d x0+d x0];
y = [y0-r y0-r y0+r y0+r y0-r y0-r y0+r y0+r];
z = [z0-h z0-h z0-h z0-h z0+h z0+h z0+h z0+h];

% Create rotation matrices
Rx = [1 0 0; 0 cos(theta_x) -sin(theta_x); 0 sin(theta_x) cos(theta_x)];
Ry = [cos(theta_y) 0 sin(theta_y); 0 1 0; -sin(theta_y) 0 cos(theta_y)];
Rz = [cos(theta_z) -sin(theta_z) 0; sin(theta_z) cos(theta_z) 0; 0 0 1];

% Rotate the coordinates of the box corners
%coords = [x + d*(cos(theta_z)-1); y - d*sin(theta_z); z];
coords = [x - rotx; y - roty ; z];

rotated_coords = Rz * Ry * Rx * coords;

% Extract the rotated coordinates
x_rotated = rotated_coords(1, :) + rotx;
y_rotated = rotated_coords(2, :) + roty;
z_rotated = rotated_coords(3, :);

faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
% Plot the box
patch('Faces', faces, 'Vertices', [x_rotated' y_rotated' z_rotated'], 'FaceColor', color);