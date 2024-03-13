classdef scara_app_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                 matlab.ui.Figure
        COORDINATEPanel          matlab.ui.container.Panel
        OpacityEditField         matlab.ui.control.NumericEditField
        OpacityEditFieldLabel    matlab.ui.control.Label
        CLEARButton              matlab.ui.control.Button
        VIEWALLButton            matlab.ui.control.Button
        Coor4CheckBox            matlab.ui.control.CheckBox
        Coor3CheckBox            matlab.ui.control.CheckBox
        Coor2CheckBox            matlab.ui.control.CheckBox
        Coor1CheckBox            matlab.ui.control.CheckBox
        Coor0CheckBox            matlab.ui.control.CheckBox
        TRAJECTORYPLANNINGPanel  matlab.ui.container.Panel
        timeAxis4                matlab.ui.control.NumericEditField
        timeEditField_4Label     matlab.ui.control.Label
        sLabel_4                 matlab.ui.control.Label
        timeAxis3                matlab.ui.control.NumericEditField
        timeEditField_3Label     matlab.ui.control.Label
        sLabel_3                 matlab.ui.control.Label
        timeAxis2                matlab.ui.control.NumericEditField
        timeEditField_2Label     matlab.ui.control.Label
        sLabel_2                 matlab.ui.control.Label
        timeAxis1                matlab.ui.control.NumericEditField
        timeEditFieldLabel       matlab.ui.control.Label
        sLabel                   matlab.ui.control.Label
        dmsLabel                 matlab.ui.control.Label
        dm2sLabel                matlab.ui.control.Label
        rad2sLabel_3             matlab.ui.control.Label
        rad2sLabel_2             matlab.ui.control.Label
        rad2sLabel               matlab.ui.control.Label
        radsLabel_3              matlab.ui.control.Label
        radsLabel_2              matlab.ui.control.Label
        radsLabel                matlab.ui.control.Label
        AXIS4Label               matlab.ui.control.Label
        AXIS3Label               matlab.ui.control.Label
        AXIS2Label               matlab.ui.control.Label
        AXIS1Label               matlab.ui.control.Label
        MaxVel_ax4               matlab.ui.control.NumericEditField
        MaxVelEditField_4Label   matlab.ui.control.Label
        MaxAcc_ax4               matlab.ui.control.NumericEditField
        MaxAccEditField_4Label   matlab.ui.control.Label
        MaxVel_ax3               matlab.ui.control.NumericEditField
        MaxVelEditField_3Label   matlab.ui.control.Label
        MaxAcc_ax3               matlab.ui.control.NumericEditField
        MaxAccEditField_3Label   matlab.ui.control.Label
        MaxVel_ax2               matlab.ui.control.NumericEditField
        MaxVelEditField_2Label   matlab.ui.control.Label
        MaxAcc_ax2               matlab.ui.control.NumericEditField
        MaxAccEditField_2Label   matlab.ui.control.Label
        Panel                    matlab.ui.container.Panel
        GOButton                 matlab.ui.control.Button
        B_Yaw                    matlab.ui.control.NumericEditField
        B_Z                      matlab.ui.control.NumericEditField
        B_Y                      matlab.ui.control.NumericEditField
        B_X                      matlab.ui.control.NumericEditField
        A_Yaw                    matlab.ui.control.NumericEditField
        XEditFieldLabel_19       matlab.ui.control.Label
        A_Z                      matlab.ui.control.NumericEditField
        XEditFieldLabel_18       matlab.ui.control.Label
        A_X                      matlab.ui.control.NumericEditField
        XEditFieldLabel_17       matlab.ui.control.Label
        BLabel                   matlab.ui.control.Label
        ALabel                   matlab.ui.control.Label
        XEditFieldLabel_2        matlab.ui.control.Label
        A_Y                      matlab.ui.control.NumericEditField
        MaxAcc_ax1               matlab.ui.control.NumericEditField
        MaxAccEditFieldLabel     matlab.ui.control.Label
        MaxVel_ax1               matlab.ui.control.NumericEditField
        MaxVelEditFieldLabel     matlab.ui.control.Label
        INVERSEKINEMATICSPanel   matlab.ui.container.Panel
        GOTOButton               matlab.ui.control.Button
        inverse_Z                matlab.ui.control.NumericEditField
        inverse_Y                matlab.ui.control.NumericEditField
        inverse_X                matlab.ui.control.NumericEditField
        inverse_Yaw              matlab.ui.control.NumericEditField
        X_inverseLabel           matlab.ui.control.Label
        Y_inverseLabel           matlab.ui.control.Label
        Z_inverseLabel           matlab.ui.control.Label
        Yaw_inverseLabel         matlab.ui.control.Label
        FORDWARDKINEMATICSPanel  matlab.ui.container.Panel
        RESETButton              matlab.ui.control.Button
        theta4SliderLabel        matlab.ui.control.Label
        theta4Slider             matlab.ui.control.Slider
        d3SliderLabel            matlab.ui.control.Label
        d3Slider                 matlab.ui.control.Slider
        theta2SliderLabel        matlab.ui.control.Label
        theta2Slider             matlab.ui.control.Slider
        theta1SliderLabel        matlab.ui.control.Label
        theta1Slider             matlab.ui.control.Slider
        EditField_theta1         matlab.ui.control.NumericEditField
        EditField_theta2         matlab.ui.control.NumericEditField
        EditField_d3             matlab.ui.control.NumericEditField
        EditField_theta4         matlab.ui.control.NumericEditField
        ENDEFFECTORPanel         matlab.ui.container.Panel
        Yawee                    matlab.ui.control.NumericEditField
        Yawee_label              matlab.ui.control.Label
        Zee                      matlab.ui.control.NumericEditField
        Yee                      matlab.ui.control.NumericEditField
        Xee                      matlab.ui.control.NumericEditField
        Zee_label                matlab.ui.control.Label
        Yee_label                matlab.ui.control.Label
        Xee_label                matlab.ui.control.Label
        WORKSPACEPanel           matlab.ui.container.Panel
        WS_ONButton              matlab.ui.control.Button
        WS_OFFButton             matlab.ui.control.Button
        UIAxes_eep               matlab.ui.control.UIAxes
        UIAxes_eea               matlab.ui.control.UIAxes
        UIAxes_eev               matlab.ui.control.UIAxes
        UIAxes_theta1dd          matlab.ui.control.UIAxes
        UIAxes_theta2dd          matlab.ui.control.UIAxes
        UIAxes_d3dd              matlab.ui.control.UIAxes
        UIAxes_theta4dd          matlab.ui.control.UIAxes
        UIAxes                   matlab.ui.control.UIAxes
        UIAxes_theta1            matlab.ui.control.UIAxes
        UIAxes_theta1d           matlab.ui.control.UIAxes
        UIAxes_theta2            matlab.ui.control.UIAxes
        UIAxes_theta2d           matlab.ui.control.UIAxes
        UIAxes_d3                matlab.ui.control.UIAxes
        UIAxes_d3d               matlab.ui.control.UIAxes
        UIAxes_theta4            matlab.ui.control.UIAxes
        UIAxes_theta4d           matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        a1 = 1.25;  alpha1 = 0;     d1 = 1.6;  theta1 = 0;
        a2 = 1.75;  alpha2 = pi;    d2 = 0;    theta2 = 0;
        a3 = 0;     alpha3 = 0;     d3 = 1;    theta3 = 0;
        a4 = 0;     alpha4 = 0;     d4 = 0;    theta4 = 0;
        
        link1_end; link2_end; link3_end; link4_end;
        EEx, EEy, EEz, EEyaw;
        
        X = 6; Y = 0; Z = 2; yaw = 0;
        ws_en = false;
        plot_coor = false;
        opacity = 1;
        showEE = false;
        matrixEE;
        theta1_mat; theta2_mat; d3_mat; theta4_mat;
        draw_inverse = true;
    end
    
    methods (Access = public)
        
        function A = DH(~,a,alpha,d,theta)
            Rot_z_theta = [cos(theta)   -sin(theta)     0       0;
                           sin(theta)    cos(theta)     0       0;
                                    0             0     1       0;
                                    0             0     0       1];
            Trans_z_d = [1 0 0 0;
                         0 1 0 0;
                         0 0 1 d;
                         0 0 0 1];
                     
            Trans_x_a = [1 0 0 a;
                         0 1 0 0;
                         0 0 1 0;
                         0 0 0 1];
                     
            Rot_x_alpha = [1             0              0     0;
                           0    cos(alpha)    -sin(alpha)     0;
                           0    sin(alpha)     cos(alpha)     0;
                           0             0              0     1];
             
            A = Rot_z_theta*Trans_z_d*Trans_x_a*Rot_x_alpha;            
        end
        
        function create_box(app,x0,y0,z0,d,r,h,theta_x,theta_y,theta_z,rotx,roty,color)
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
            hold(app.UIAxes,'on')
            a = patch(app.UIAxes,'Faces', faces, 'Vertices', [x_rotated' y_rotated' z_rotated'], 'FaceColor', color);
            a.FaceAlpha = app.opacity;
            
        end
        
        function create_cylinder(app,x,y,z,r,d,color)
            [X0,Y0,Z0] = cylinder(r,100);
            %hold(app.UIAxes,'on')
            a = surf(app.UIAxes,X0+x,Y0+y,Z0*d+z,'facecolor',color,'LineStyle','none');
            a.FaceAlpha = app.opacity;
            a = fill3(app.UIAxes,X0(1,:) + x, Y0(1,:) + y, Z0(1,:) + z, color);
            a.FaceAlpha = app.opacity;            
            a = fill3(app.UIAxes,X0(2,:) + x, Y0(2,:) + y, Z0(2,:)*d + z, color);  
            a.FaceAlpha = app.opacity;            
        end
        
        function out = main(app,theta1,theta2,d3,theta4)
            
            app.theta1 = theta1;
            app.theta2 = theta2;
            app.d3 = d3;
            app.theta4 = theta4;
            
            %% D-H parameters
            % Link 1 D-H
            A01 = DH(app,app.a1,app.alpha1,app.d1,theta1);
            
            % Link 2 D-H
            A12 = DH(app,app.a2,app.alpha2,app.d2,theta2);
            
            % Link 3 D-H
            A23 = DH(app,app.a3,app.alpha3,d3,app.theta3);
            
            % Link 4 D-H
            A34 = DH(app,app.a4,app.alpha4,app.d4,theta4);
            
            %%% Get position of joints
            % Position of joint 1
            app.link1_end = A01(1:3,4);
            
            % Position of joint 2
            A02 = A01*A12;
            app.link2_end = A02(1:3,4);
            
            % Positon of joint 3
            A03 = A02*A23;
            app.link3_end = A03(1:3,4);
            
            % Position of joint 4
            A04 = A03*A34;
            app.EEyaw = app.theta1 + app.theta2 + app.theta3 + app.theta4;            
            app.link4_end = A04(1:3,4);
            app.EEx = app.link4_end(1,1);
            app.EEy = app.link4_end(2,1);
            app.EEz = app.d1 - app.d3;
            
            p = ((app.EEx - app.A_X.Value)^2 + (app.EEy - app.A_Y.Value)^2 + (app.EEz - app.A_Z.Value)^2)^(1/2);
            pd = ((app.EEx - app.A_X.Value) + (app.EEy - app.A_Y.Value) + (app.EEz - app.A_Z.Value))*((app.EEx - app.A_X.Value)^2 + (app.EEy - app.A_Y.Value)^2 + (app.EEz - app.A_Z.Value)^2)^(-1/2);
            pdd = (-1)*((app.EEx - app.A_X.Value) + (app.EEy - app.A_Y.Value) + (app.EEz - app.A_Z.Value))*((app.EEx - app.A_X.Value)^2 + (app.EEy - app.A_Y.Value)^2 + (app.EEz - app.A_Z.Value)^2)^(-3/2);
            
            out = [p;pd;pdd];

            
            app.matrixEE(end+1,:) = [app.EEx app.EEy app.EEz];
%             app.inverse_X.Value = app.link4_end(1,1);
%             app.inverse_Y.Value = app.link4_end(2,1);
%             app.inverse_Z.Value = app.d1 - app.d3;
%             app.inverse_Yaw.Value = app.Yawee_label;

            app.Xee.Value = app.EEx;
            app.Yee.Value = app.EEy;
            app.Zee.Value = app.EEz;
            app.Yawee.Value = app.EEyaw;
            
            % End effector position and orientation
            %position_ee = A04(1:3,4);
            %orientation_ee = A04(1:3,1:3);
            draw(app,app.link1_end,app.link2_end,app.link3_end,app.link4_end); 
            
            %view(3);            
            
        end
        
        function out = inverse(app,X,Y,Z,yaw)
            
            
            
            app.theta2 = acos((X^2 + Y^2 - app.a1^2 - app.a2^2)/(2*app.a1*app.a2));            
            app.theta1 = atan(Y/X) - atan((app.a2*sin(app.theta2))/(app.a1 + app.a2*cos(app.theta2)));
                      
            
            if ((X<0 && Y>0))
                if app.theta1 >= 0
                    if app.theta1 - pi >= -125*pi/180 || app.theta1 - pi <= 125*pi/180
                        app.theta1 = app.theta1 - pi;
                    else
                        app.theta2 = -app.theta2;
                        app.theta1 = atan(Y/X) - atan((app.a2*sin(app.theta2))/(app.a1 + app.a2*cos(app.theta2)));
                        app.theta1 = app.theta1 - pi;
                    end
                elseif app.theta1 < 0
                    if app.theta1 + pi >= -125*pi/180 || app.theta1 + pi <= 125*pi/180
                        app.theta1 = app.theta1 + pi;
                    else
                        app.theta2 = -app.theta2;
                        app.theta1 = atan(Y/X) - atan((app.a2*sin(app.theta2))/(app.a1 + app.a2*cos(app.theta2)));
                        app.theta1 = app.theta1 + pi;
                    end
                end
            end
            
            if ((X<0 && Y<0))
                app.theta1 = app.theta1 - pi;
                if app.theta1 < -125*pi/180 || app.theta1 > 125*pi/180
                    app.theta2 = -app.theta2;
                    app.theta1 = atan(Y/X) - atan((app.a2*sin(app.theta2))/(app.a1 + app.a2*cos(app.theta2)));
                    app.theta1 = app.theta1 - pi;
                end
            end
            
            if (app.theta1 < -125*pi/180 || app.theta1 > 125*pi/180)
                app.theta2 = -app.theta2;
                app.theta1 = atan(Y/X) - atan((app.a2*sin(app.theta2))/(app.a1 + app.a2*cos(app.theta2)));
            end     
            
            if (Y==0 && X<0)
                if app.theta1 < 0
                    app.theta1 = app.theta1 + pi;
                else
                    app.theta1 = app.theta1 - pi;
                end
            end
            

            app.d3 = app.d1 - Z;
            app.theta4 = yaw - app.theta1 - app.theta2;
            
            A01 = DH(app,app.a1,app.alpha1,app.d1,app.theta1);
            A12 = DH(app,app.a2,app.alpha2,app.d2,app.theta2);
            A23 = DH(app,app.a3,app.alpha3,app.d3,app.theta3);
            A34 = DH(app,app.a4,app.alpha4,app.d4,app.theta4);
            % Position of joint 1
            app.link1_end = A01(1:3,4);
            
            % Position of joint 2
            A02 = A01*A12;
            app.link2_end = A02(1:3,4);
            
            % Positon of joint 3
            A03 = A02*A23;
            app.link3_end = A03(1:3,4);
            
            % Position of joint 4
            A04 = A03*A34;
            app.link4_end = A04(1:3,4);
            
            out = [app.theta1 app.theta2 app.d3 app.theta4];
            
            app.EditField_theta1.Value = app.theta1;
            app.EditField_theta2.Value = app.theta2;
            app.EditField_d3.Value     = app.d3;
            app.EditField_theta4.Value = app.theta4;
            
            app.theta1Slider.Value = app.theta1;
            app.theta2Slider.Value = app.theta2;
            app.d3Slider.Value     = app.d3;
            app.theta4Slider.Value = app.theta4;
            
%             if app.draw_inverse == true
%                 draw(app,app.link1_end,app.link2_end,app.link3_end,app.link4_end); 
%             end
            
        end
        
        function draw(app,link1_end,link2_end,link3_end,link4_end)
            cla(app.UIAxes) 
            if app.ws_en == true
                app.WS_ONButtonPushed;
            end
            if app.Coor0CheckBox.Value == 1
                plot_coordinate(app,0,0,app.d1,app.theta1,1);
            end
            if app.Coor1CheckBox.Value == 1
                plot_coordinate(app,app.link1_end(1),app.link1_end(2),app.d1 + app.d2,app.theta2 + app.theta1,1);
            end
            if app.Coor2CheckBox.Value == 1
                plot_coordinate(app,app.link2_end(1),app.link2_end(2),app.d1 + app.d2,app.theta3 + app.theta2 + app.theta1,1);
            end
            if app.Coor3CheckBox.Value == 1
                plot_coordinate(app,app.link2_end(1),app.link2_end(2),app.link3_end(3),app.theta3 + app.theta2 + app.theta1,1);
            end
            if app.Coor4CheckBox.Value == 1
                plot_coordinate(app,app.link3_end(1),app.link3_end(2),app.link4_end(3),app.EEyaw+pi,-1);
            end
            
            if app.showEE == true
                plot3(app.UIAxes,app.matrixEE(:,1),app.matrixEE(:,2),app.matrixEE(:,3));
            end

            %% Draw objects between joints
            % Origin to start of link 1
            create_cylinder(app,0,0,0,0.5,0.2,'black');
            %line( [0 0], [0 0], [0 d1], 'LineWidth ',5,'Color','yellow');
            create_cylinder(app,0,0,0,0.2,app.d1,'r');
            create_cylinder(app,0,0,link1_end(3)-0.2,0.25,0.35,'blue');
            
            % start of link 1 to end of link 1
            %line([0 link1_end(1)], [0 link1_end(2)], [d1 link1_end(3)], 'LineWidth',5,'Color','red');
            create_box(app,0,0,link1_end(3),app.a1,0.2,0.1,0,0,app.theta1,0,0,'blue');
            create_cylinder(app,link1_end(1),link1_end(2),link1_end(3)-0.1,0.25,0.25,'blue');
            create_cylinder(app,link1_end(1),link1_end(2),link1_end(3)+0.1,0.225,0.3,'green');
            
            % End of link 1 to end of link 2
            %line( [link1_end(1) link2_end(1)], [link1_end(2) link2_end(2)], [link1_end(3) link2_end(3)]+0.5, 'LineWidth',5,'Color','green');
            create_box(app,link1_end(1),link1_end(2),link1_end(3)+0.25,app.a2,0.2,0.1,0,0,app.theta1+app.theta2,link1_end(1),link1_end(2),'green');
            create_cylinder(app,link2_end(1),link2_end(2),link2_end(3)+0.1,0.225,0.3,'green');
            
            % End of link 2 to end of link 3
            %line( [link2_end(1) link3_end(1)], [link2_end(2) link3_end(2)], [link2_end(3)+5-d3 link3_end(3)], 'LineWidth',5,'Color','black');
            create_cylinder(app,link3_end(1),link3_end(2),link3_end(3),0.125,2.5,'black');
            
            % End of link 3 to end of link 4
            %create_cylinder(app,link4_end(1),link4_end(2),link4_end(3),0.15,0.2,'yellow');
            
            create_box(app,link4_end(1)-0.25,link4_end(2),link3_end(3),0.5,0.05,0.1,0,0,app.EEyaw,link4_end(1),link4_end(2),'black');
            create_cylinder(app,link4_end(1)-0.25*cos(app.EEyaw),link4_end(2)-0.25*sin(app.EEyaw),link4_end(3)-0.3,0.05,0.4,'black');
            create_cylinder(app,link4_end(1)+0.25*cos(app.EEyaw),link4_end(2)+0.25*sin(app.EEyaw),link4_end(3)-0.3,0.05,0.4,'black');
            
            %% PLOT
            xlabel(app.UIAxes,'x');
            ylabel(app.UIAxes,'y');
            zlabel(app.UIAxes,'z');
            grid(app.UIAxes,"on");
            axis(app.UIAxes,'equal');
            xlim(app.UIAxes,[-3.5 3.5]);
            ylim(app.UIAxes,[-3.5 3.5]);
            zlim(app.UIAxes,[0 4]);            
        end
        
        function plot_coordinate(app,x,y,z,theta,z1)
            line(app.UIAxes,[x x+cos(theta)], [y y+sin(theta)], [z z], 'LineWidth',5,'Color','red');
            theta = theta + pi/2;
            line(app.UIAxes,[x x+cos(theta)], [y y+sin(theta)], [z z], 'LineWidth',5,'Color','green');
            line(app.UIAxes,[x x], [y y], [z z+z1], 'LineWidth',5,'Color','blue');
        end
        
        function out = create_trajectory(~,q1, q2, a_max, v_max)
            
            q_max = q2 - q1;
            reverse = false;
            
            if q_max < 0
                q_max = -q_max;
                reverse = true;
            end
            
            if a_max < 2*v_max^2/q_max
                a_max = 2*v_max^2/q_max;
            end
            
            t1 = v_max / a_max;
            t2 = t1 + t1;
            t3 = q_max/v_max;
            k1 = a_max/t1;
            t4 = t3 + t1;
            t5 = t4 + t1;
        
            n = 500;
            total_step = floor(t5*n);
            q(1:total_step) = 0;
            v(1:total_step) = 0;
            a(1:total_step) = 0;
            t(1:total_step) = 0;
        
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
                    v(i) = -((1/2)*k1*((t(i)-t1)^2)) + a_max*(t(i)-t1) + v(i1);
                    q(i) = -(1/6)*k1*((t(i)-t1)^3) + a_max*(1/2)*((t(i)-t1)^2) + v(i1)*(t(i)-t1) + q(i1);
                    i2 = i;
                    i3 = i2;
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
                end
            end
            if reverse == true
                out = [t;-a;-v;-q];
            else
                out = [t;a;v;q];
            end
        end

    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.a1 = 1.25;  app.alpha1 = 0;  app.d1 = 1.6;  app.theta1 = 0;
            app.a2 = 1.75;  app.alpha2 = pi; app.d2 = 0;    app.theta2 = 0;
            app.a3 = 0;     app.alpha3 = 0;  app.d3 = 1;    app.theta3 = 0;
            app.a4 = 0;     app.alpha4 = 0;  app.d4 = 0;    app.theta4 = 0;
            
            %app.X = 6; app.Y = 0; app.Z = 2; app.yaw = 0;    
            
            main(app,app.theta1,app.theta2,app.d3,app.theta4);
            main(app,app.theta1,app.theta2,app.d3,app.theta4);            
        end

        % Value changing function: theta1Slider
        function theta1SliderValueChanging(app, event)
            changingValue = event.Value;
            app.theta1 = changingValue;
            app.EditField_theta1.Value = app.theta1;
            main(app,app.theta1,app.theta2,app.d3,app.theta4);
        end

        % Value changing function: theta2Slider
        function theta2SliderValueChanging(app, event)
            changingValue = event.Value;
            app.theta2 = changingValue;
            app.EditField_theta2.Value = app.theta2;
            main(app,app.theta1,app.theta2,app.d3,app.theta4);
        end

        % Value changing function: d3Slider
        function d3SliderValueChanging(app, event)
            changingValue = event.Value;
            app.d3 = changingValue;
            app.EditField_d3.Value = app.d3;
            main(app,app.theta1,app.theta2,app.d3,app.theta4);
        end

        % Value changing function: theta4Slider
        function theta4SliderValueChanging(app, event)
            changingValue = event.Value;
            app.theta4 = changingValue;
            app.EditField_theta4.Value = app.theta4;
            main(app,app.theta1,app.theta2,app.d3,app.theta4);
        end

        % Value changed function: inverse_X
        function inverse_XValueChanged(app, event)
            value = app.inverse_X.Value;
            app.X = value;
        end

        % Value changed function: inverse_Y
        function inverse_YValueChanged(app, event)
            value = app.inverse_Y.Value;
            app.Y = value;
        end

        % Value changed function: inverse_Z
        function inverse_ZValueChanged(app, event)
            value = app.inverse_Z.Value;
            app.Z = value;
        end

        % Value changed function: inverse_Yaw
        function inverse_YawValueChanged(app, event)
            value = app.inverse_Yaw.Value;
            app.yaw = value;
        end

        % Button pushed function: GOTOButton
        function GOTOButtonPushed(app, event)
            inverse(app,app.inverse_X.Value,app.inverse_Y.Value,app.inverse_Z.Value,app.inverse_Yaw.Value)
            draw(app,app.link1_end,app.link2_end,app.link3_end,app.link4_end)
        end

        % Value changed function: EditField_theta1
        function EditField_theta1ValueChanged(app, event)
            value = app.EditField_theta1.Value;
            app.theta1 = value;
            app.theta1Slider.Value = value;
        end

        % Value changed function: EditField_theta2
        function EditField_theta2ValueChanged(app, event)
            value = app.EditField_theta2.Value;
            app.theta2 = value;
            app.theta2Slider.Value = value;
        end

        % Value changed function: EditField_d3
        function EditField_d3ValueChanged(app, event)
            value = app.EditField_d3.Value;
            app.d3 = value;
            app.d3Slider.Value = value;
        end

        % Value changed function: EditField_theta4
        function EditField_theta4ValueChanged(app, event)
            value = app.EditField_theta4.Value;
            app.theta4 = value;
            app.theta4Slider.Value = value;
        end

        % Button pushed function: WS_ONButton
        function WS_ONButtonPushed(app, event)
            app.ws_en = true;
            
            x = 0; y = 0; z = 0;
            [X0,Y0,Z0] = cylinder(3,1000);
            op = 0.2;
            
            p = surf(app.UIAxes,X0(:,1:347) + x, Y0(:,1:347) + y, Z0(:,1:347)*1.6 + z,"EdgeColor","none");
            p.FaceAlpha = op;
            %hold on;
            p = surf(app.UIAxes,X0(:,654:1001) + x, Y0(:,654:1001) + y, Z0(:,654:1001)*1.6 + z,"EdgeColor","none");
            p.FaceAlpha = op;
            
            t = linspace(-125*pi/180,125*pi/180);
            a = 1.25; b = 1.75;
            x0=0; y0=0;              % circles center
            c = sqrt(a^2 + b^2 -2*a*b*cos(pi-145*pi/180));
            rin = c;%sqrt(1.25^2 + 1.75^2 - 2*1.25*1.75*cos(pi-2)); 
            rout = 3;   % radii sizes= 
            p = patch(app.UIAxes,[x0+rout*cos(t),x0+flip(rin*cos(t))],[y0+rout*sin(t),y0+flip(rin*sin(t))],(1.6)*ones(1,2*length(rout*cos(t))),'r','linestyle','none');
            p.FaceAlpha = op;
            p = patch(app.UIAxes,[x0+rout*cos(t),x0+flip(rin*cos(t))],[y0+rout*sin(t),y0+flip(rin*sin(t))],(0)*ones(1,2*length(rout*cos(t))),'r','linestyle','none');
            p.FaceAlpha = op;
            
            c = sqrt(a^2 + b^2 -2*a*b*cos(pi-145*pi/180));
            [X1,Y1,Z1] = cylinder(c,1000);
            p = surf(app.UIAxes, X1 + x, Y1 + y, Z1*1.6 + z,"EdgeColor","none");
            p.FaceAlpha = op;
            
            y = sin(55*pi/180)*1.25; x = -cos(55*pi/180)*1.25; z = 0;
            h = sin(55*pi/180)*1.25;
            alpha = acos(h/1.75) + (35*pi/180);
            beta = pi - alpha;
            [X0,Y0,Z0] = cylinder(1.75,1000);
            j = 601;
            i = 345;
            p = surf(app.UIAxes, X0(:,i:j) + x, Y0(:,i:j) + y, Z0(:,i:j)*1.6 + z,"EdgeColor","none");
            p.FaceAlpha = op;
            
            t = linspace(125*pi/180,beta+125*pi/180);
            t1 = linspace(125*pi/180,pi);
            x0=-sqrt(1.25^2 - h^2); y0=h;              % circles center
            rin = c;%sqrt(1.25^2 + 1.75^2 - 2*1.25*1.75*cos(pi-2)); 
            rout = 1.75;   % radii sizes
            p = patch(app.UIAxes, [x0+rout*cos(t),0+flip(rin*cos(t1))],[y0+rout*sin(t),0+flip(rin*sin(t1))],(1.6)*ones(1,2*length(rout*cos(t))),'red','linestyle','none');
            p.FaceAlpha = op;
            p = patch(app.UIAxes, [x0+rout*cos(t),0+flip(rin*cos(t1))],[y0+rout*sin(t),0+flip(rin*sin(t1))],(0)*ones(1,2*length(rout*cos(t))),'red','linestyle','none');
            p.FaceAlpha = op;
            
            y = -y;
            [X0,Y0,Z0] = cylinder(1.75,1000);
            j = 670;
            i = 400;
            p = surf(app.UIAxes,X0(:,i:j) + x, Y0(:,i:j) + y, Z0(:,i:j)*1.6 + z,"EdgeColor","none");
            p.FaceAlpha = op;
            
            t = linspace(alpha+55*pi/180,alpha+beta+55*pi/180);
            t1 = linspace(pi,alpha+beta+55*pi/180);
            x0=-sqrt(1.25^2 - h^2); y0= -h;              % circles center
            rin = c;%sqrt(1.25^2 + 1.75^2 - 2*1.25*1.75*cos(pi-2)); 
            rout = 1.75;   % radii sizes
            tx = t;
            p = patch(app.UIAxes,[x0+rout*cos(tx),0+flip(rin*cos(t1))],[y0+rout*sin(tx),0+flip(rin*sin(t1))],(1.6)*ones(1,2*length(rout*cos(tx))),'red','linestyle','none');
            p.FaceAlpha = op;
            p = patch(app.UIAxes,[x0+rout*cos(tx),0+flip(rin*cos(t1))],[y0+rout*sin(tx),0+flip(rin*sin(t1))],(0)*ones(1,2*length(rout*cos(tx))),'red','linestyle','none');
            p.FaceAlpha = op;
        end

        % Button pushed function: WS_OFFButton
        function WS_OFFButtonPushed(app, event)
            cla(app.UIAxes) 
            app.ws_en = false;
            draw(app,app.link1_end,app.link2_end,app.link3_end,app.link4_end); 
        end

        % Callback function
        function ButtonGroupSelectionChanged(app, event)
            selectedButton = app.ButtonGroup.SelectedObject;
            if (selectedButton == app.LinearPathButton)
                app.Panel.Visible = 'on';
                app.CircularPathPanel.Visible = 'off';
                app.SingularityPathPanel.Visible = 'off';
            end
            if (selectedButton == app.CircularPathButton)
                app.Panel.Visible = 'off';
                app.CircularPathPanel.Visible = 'on';
                app.SingularityPathPanel.Visible = 'off';
                
            end
            if (selectedButton == app.SingularityPathButton)
                app.Panel.Visible = 'off';
                app.CircularPathPanel.Visible = 'off';
                app.SingularityPathPanel.Visible = 'on';            
            end            
        end

        % Value changed function: OpacityEditField
        function OpacityEditFieldValueChanged(app, event)
            value = app.OpacityEditField.Value;
            app.opacity = value;
            draw(app,app.link1_end,app.link2_end,app.link3_end,app.link4_end);
        end

        % Value changed function: Coor0CheckBox
        function Coor0CheckBoxValueChanged(app, event)
            value = app.Coor0CheckBox.Value;
            if value == 1
                plot_coordinate(app,0,0,app.d1,app.theta1,1);
            else
                draw(app,app.link1_end,app.link2_end,app.link3_end,app.link4_end); 
            end
        end

        % Value changed function: Coor1CheckBox
        function Coor1CheckBoxValueChanged(app, event)
            value = app.Coor1CheckBox.Value;
            if value == 1
                plot_coordinate(app,app.link1_end(1),app.link1_end(2),app.d1 + app.d2,app.theta2 + app.theta1,1);
            else
                draw(app,app.link1_end,app.link2_end,app.link3_end,app.link4_end); 
            end            
        end

        % Value changed function: Coor2CheckBox
        function Coor2CheckBoxValueChanged(app, event)
            value = app.Coor2CheckBox.Value;
            if value == 1
                plot_coordinate(app,app.link2_end(1),app.link2_end(2),app.d1 + app.d2,app.theta3 + app.theta2 + app.theta1,1);
            else
                draw(app,app.link1_end,app.link2_end,app.link3_end,app.link4_end); 
            end              
        end

        % Value changed function: Coor3CheckBox
        function Coor3CheckBoxValueChanged(app, event)
            value = app.Coor3CheckBox.Value;
            if value == 1
                plot_coordinate(app,app.link2_end(1),app.link2_end(2),app.link3_end(3),app.theta3 + app.theta2 + app.theta1,1);
            else
                draw(app,app.link1_end,app.link2_end,app.link3_end,app.link4_end); 
            end              
        end

        % Value changed function: Coor4CheckBox
        function Coor4CheckBoxValueChanged(app, event)
            value = app.Coor4CheckBox.Value;
            if value == 1
                plot_coordinate(app,app.link3_end(1),app.link3_end(2),app.link4_end(3),app.EEyaw+pi,-1);
            else
                draw(app,app.link1_end,app.link2_end,app.link3_end,app.link4_end); 
            end              
        end

        % Button pushed function: VIEWALLButton
        function VIEWALLButtonPushed(app, event)
            app.Coor0CheckBox.Value = 1;
            plot_coordinate(app,0,0,app.d1,app.theta1,1);
            app.Coor1CheckBox.Value = 1;
            plot_coordinate(app,app.link1_end(1),app.link1_end(2),app.d1 + app.d2,app.theta2 + app.theta1,1);
            app.Coor2CheckBox.Value = 1;
            plot_coordinate(app,app.link2_end(1),app.link2_end(2),app.d1 + app.d2,app.theta3 + app.theta2 + app.theta1,1);
            app.Coor3CheckBox.Value = 1;
            plot_coordinate(app,app.link2_end(1),app.link2_end(2),app.link3_end(3),app.theta3 + app.theta2 + app.theta1,1);
            app.Coor4CheckBox.Value = 1;
            plot_coordinate(app,app.link3_end(1),app.link3_end(2),app.link4_end(3),app.EEyaw+pi,-1);
        end

        % Button pushed function: CLEARButton
        function CLEARButtonPushed(app, event)
            app.Coor0CheckBox.Value = 0;
            app.Coor1CheckBox.Value = 0;
            app.Coor2CheckBox.Value = 0;
            app.Coor3CheckBox.Value = 0;
            app.Coor4CheckBox.Value = 0;   
            draw(app,app.link1_end,app.link2_end,app.link3_end,app.link4_end); 
        end

        % Button pushed function: RESETButton
        function RESETButtonPushed(app, event)
            app.theta1 = 0;
            app.theta2 = 0;
            app.d3 = 1;
            app.theta4 = 0;
            app.EditField_theta1.Value = 0;
            app.EditField_theta2.Value = 0;
            app.EditField_d3.Value = 0;
            app.EditField_theta4.Value = 0;
            app.theta1Slider.Value = 0;
            app.theta2Slider.Value = 0;
            app.d3Slider.Value = 01;
            app.theta4Slider.Value = 0;
            main(app,app.theta1,app.theta2,app.d3,app.theta4);
        end

        % Callback function
        function EEButtonValueChanged(app, event)
            value = app.EEButton.Value;
            if value == 1
                app.showEE = true;
            else
                app.showEE = false;
            end
        end

        % Button pushed function: GOButton
        function GOButtonPushed(app, event)
            cla(app.UIAxes_theta1);
            cla(app.UIAxes_theta1d);
            cla(app.UIAxes_theta1dd);
            
            cla(app.UIAxes_theta2);
            cla(app.UIAxes_theta2d);
            cla(app.UIAxes_theta2dd);
            
            cla(app.UIAxes_d3);
            cla(app.UIAxes_d3d);
            cla(app.UIAxes_d3dd);
            
            cla(app.UIAxes_theta4);
            cla(app.UIAxes_theta4d);
            cla(app.UIAxes_theta4dd);
            
            cla(app.UIAxes_eep);
            cla(app.UIAxes_eev);
            cla(app.UIAxes_eea);
            
            app.matrixEE = [];
            
            
            
            value1 = inverse(app,app.A_X.Value,app.A_Y.Value,app.A_Z.Value,app.A_Yaw.Value);
            app.draw_inverse = false;
            value2 = inverse(app,app.B_X.Value,app.B_Y.Value,app.B_Z.Value,app.B_Yaw.Value);
            app.draw_inverse = true;
            
            axis = [value1(1) value2(1) app.MaxAcc_ax1.Value app.MaxVel_ax1.Value; 
                    value1(2) value2(2) app.MaxAcc_ax2.Value app.MaxVel_ax2.Value; 
                    value1(3) value2(3) app.MaxAcc_ax3.Value app.MaxVel_ax3.Value; 
                    value1(4) value2(4) app.MaxAcc_ax4.Value app.MaxVel_ax4.Value];
            %axis = [0 pi 115 660*pi/180; -pi/6 pi 115 11.5; 0.5 1.4 112 11.2; 0 pi/6 262 1500*pi/180]; % qmax qmin a_max v_min
            
            qmin = [axis(1,1)  axis(2,1) axis(3,1) axis(4,1)];
            qmax = [axis(1,2)  axis(2,2) axis(3,2) axis(4,2)];
            amax = [axis(1,3)  axis(2,3) axis(3,3) axis(4,3)];
            vmax = [axis(1,4)  axis(2,4) axis(3,4) axis(4,4)];
    

%             j = 1;
%             m = app.create_trajectory(qmin(j),qmax(j),amax(j),vmax(j));
%             plot(app.UIAxes_theta1,m(1,:),m(2,:));
%             plot(app.UIAxes_theta1d,m(1,:),m(3,:));
%             plot(app.UIAxes_theta1dd,m(1,:),m(4,:));
%             
%             j = 2;
%             m = app.create_trajectory(qmin(j),qmax(j),amax(j),vmax(j));
%             plot(app.UIAxes_theta2,m(1,:),m(2,:));
%             plot(app.UIAxes_theta2d,m(1,:),m(3,:));
%             plot(app.UIAxes_theta2dd,m(1,:),m(4,:));
%             
%             j = 3;
%             m = app.create_trajectory(qmin(j),qmax(j),amax(j),vmax(j));
%             plot(app.UIAxes_d3,m(1,:),m(2,:));
%             plot(app.UIAxes_d3d,m(1,:),m(3,:));
%             plot(app.UIAxes_d3dd,m(1,:),m(4,:));   
%             
%             j = 4;
%             m = app.create_trajectory(qmin(j),qmax(j),amax(j),vmax(j));
%             plot(app.UIAxes_theta4,m(1,:),m(2,:));
%             plot(app.UIAxes_theta4d,m(1,:),m(3,:));
%             plot(app.UIAxes_theta4dd,m(1,:),m(4,:));

            m1 = app.create_trajectory(qmin(1),qmax(1),amax(1),vmax(1));
            t1 = m1(1,end);
%             tmax = t1;
            app.timeAxis1.Value = t1;
            m2 = app.create_trajectory(qmin(2),qmax(2),amax(2),vmax(2));
            t2 = m2(1,end);
%             if t2 > tmax
%                 tmax = t2;
%             end
            app.timeAxis2.Value = t2;
            m3 = app.create_trajectory(qmin(3),qmax(3),amax(3),vmax(3));
            t3 = m3(1,end);
%             if t3 > tmax
%                 tmax = t3;
%             end
            app.timeAxis3.Value = t3;
            m4 = app.create_trajectory(qmin(4),qmax(4),amax(4),vmax(4));
            t4 = m4(1,end);
%             if t4 > tmax
%                 tmax = t4;
%             end
            app.timeAxis4.Value = t4;
            
            pp = [];
            
            n = 500; i = 2;
            while (i < t1*n || i < t2*n || i < t3*n || i < t4*n)
                if i <= t1*n
                    plot(app.UIAxes_theta1dd,m1(1,i-1:i),m1(2,i-1:i));
                    hold(app.UIAxes_theta1dd,"on")
                    plot(app.UIAxes_theta1d,m1(1,i-1:i),m1(3,i-1:i));
                    hold(app.UIAxes_theta1d,"on")
                    plot(app.UIAxes_theta1,m1(1,i-1:i),m1(4,i-1:i)); 
                    hold(app.UIAxes_theta1,"on")
                    end1 = i;
                    %app.timeAxis1.Value = t1*i;
                else
                    end1 = floor(t1*n);
                end
                
                if i <= t2*n
                    plot(app.UIAxes_theta2dd,m2(1,i-1:i),m2(2,i-1:i));
                    hold(app.UIAxes_theta2dd,"on")
                    plot(app.UIAxes_theta2d,m2(1,i-1:i),m2(3,i-1:i));
                    hold(app.UIAxes_theta2d,"on")
                    plot(app.UIAxes_theta2,m2(1,i-1:i),m2(4,i-1:i));
                    hold(app.UIAxes_theta2,"on")
                    end2 = i;
                    %app.timeAxis2.Value = t2*(i/n);
                else
                    end2 = floor(t2*n);
                end
                
                if i <= t3*n
                    plot(app.UIAxes_d3dd,m3(1,i-1:i),m3(2,i-1:i));
                    hold(app.UIAxes_d3dd,"on")
                    plot(app.UIAxes_d3d,m3(1,i-1:i),m3(3,i-1:i));
                    hold(app.UIAxes_d3d,"on")
                    plot(app.UIAxes_d3,m3(1,i-1:i),m3(4,i-1:i)); 
                    hold(app.UIAxes_d3,"on")
                    %app.timeAxis3.Value = t3*(i/n);
                    end3 = i;
                else
                    end3 = floor(t3*n);
                end
                
                if i <= t4*n
                    plot(app.UIAxes_theta4dd,m4(1,i-1:i),m4(2,i-1:i));
                    hold(app.UIAxes_theta4dd,"on")
                    plot(app.UIAxes_theta4d,m4(1,i-1:i),m4(3,i-1:i));
                    hold(app.UIAxes_theta4d,"on")
                    plot(app.UIAxes_theta4,m4(1,i-1:i),m4(4,i-1:i));
                    hold(app.UIAxes_theta4,"on")
                    %app.timeAxis4.Value = t4*i;
                    end4 = i;
                else
                    end4 = floor(t4*n);
                end
                app.showEE = true;
                
                pp(:,2) = main(app,m1(4,end1) + qmin(1), m2(4,end2) + qmin(2), m3(4,end3) + qmin(3), m4(4,end4) + qmin(4));
                
                plot(app.UIAxes_eep,((i-1:i)/n),pp(1,1:2));
                hold(app.UIAxes_eep,"on")
%                 plot(app.UIAxes_eev,(tmax/n-(i-1))*(i-1:i),pp(2,1:2));
%                 hold(app.UIAxes_eev,"on")                
%                 plot(app.UIAxes_eea,(tmax/n-(i-1))*(i-1:i),pp(3,1:2));
%                 hold(app.UIAxes_eea,"on") 
                
                pp(:,1) = pp(:,2);
                i = i + 1;
                pause(1/1000);
            end
            app.showEE = false;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.9686 0.9686 0.9686];
            app.UIFigure.Position = [100 100 1917 1029];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.WindowStyle = 'modal';
            app.UIFigure.WindowState = 'maximized';

            % Create UIAxes_theta4d
            app.UIAxes_theta4d = uiaxes(app.UIFigure);
            title(app.UIAxes_theta4d, 'theta4 dot')
            xlabel(app.UIAxes_theta4d, 'X')
            ylabel(app.UIAxes_theta4d, 'Y')
            zlabel(app.UIAxes_theta4d, 'Z')
            app.UIAxes_theta4d.Position = [1497 265 205 166];

            % Create UIAxes_theta4
            app.UIAxes_theta4 = uiaxes(app.UIFigure);
            title(app.UIAxes_theta4, 'theta4')
            xlabel(app.UIAxes_theta4, 'X')
            ylabel(app.UIAxes_theta4, 'Y')
            zlabel(app.UIAxes_theta4, 'Z')
            app.UIAxes_theta4.Position = [1289 265 205 166];

            % Create UIAxes_d3d
            app.UIAxes_d3d = uiaxes(app.UIFigure);
            title(app.UIAxes_d3d, 'd3 dot')
            xlabel(app.UIAxes_d3d, 'X')
            ylabel(app.UIAxes_d3d, 'Y')
            zlabel(app.UIAxes_d3d, 'Z')
            app.UIAxes_d3d.Position = [1497 457 205 166];

            % Create UIAxes_d3
            app.UIAxes_d3 = uiaxes(app.UIFigure);
            title(app.UIAxes_d3, 'd3')
            xlabel(app.UIAxes_d3, 'X')
            ylabel(app.UIAxes_d3, 'Y')
            zlabel(app.UIAxes_d3, 'Z')
            app.UIAxes_d3.Position = [1289 457 205 166];

            % Create UIAxes_theta2d
            app.UIAxes_theta2d = uiaxes(app.UIFigure);
            title(app.UIAxes_theta2d, 'theta2 dot')
            xlabel(app.UIAxes_theta2d, 'X')
            ylabel(app.UIAxes_theta2d, 'Y')
            zlabel(app.UIAxes_theta2d, 'Z')
            app.UIAxes_theta2d.Position = [1497 648 205 166];

            % Create UIAxes_theta2
            app.UIAxes_theta2 = uiaxes(app.UIFigure);
            title(app.UIAxes_theta2, 'theta2')
            xlabel(app.UIAxes_theta2, 'X')
            ylabel(app.UIAxes_theta2, 'Y')
            zlabel(app.UIAxes_theta2, 'Z')
            app.UIAxes_theta2.Position = [1289 648 205 166];

            % Create UIAxes_theta1d
            app.UIAxes_theta1d = uiaxes(app.UIFigure);
            title(app.UIAxes_theta1d, 'theta1 dot')
            xlabel(app.UIAxes_theta1d, 'X')
            ylabel(app.UIAxes_theta1d, 'Y')
            zlabel(app.UIAxes_theta1d, 'Z')
            app.UIAxes_theta1d.Position = [1497 839 205 166];

            % Create UIAxes_theta1
            app.UIAxes_theta1 = uiaxes(app.UIFigure);
            title(app.UIAxes_theta1, 'theta1')
            xlabel(app.UIAxes_theta1, 'X')
            ylabel(app.UIAxes_theta1, 'Y')
            zlabel(app.UIAxes_theta1, 'Z')
            app.UIAxes_theta1.Position = [1289 839 205 166];

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'MY SCARA')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.FontWeight = 'bold';
            app.UIAxes.FontSize = 12;
            app.UIAxes.Position = [612 429 694 577];

            % Create UIAxes_theta4dd
            app.UIAxes_theta4dd = uiaxes(app.UIFigure);
            title(app.UIAxes_theta4dd, 'theta4 dot dot')
            xlabel(app.UIAxes_theta4dd, 'X')
            ylabel(app.UIAxes_theta4dd, 'Y')
            zlabel(app.UIAxes_theta4dd, 'Z')
            app.UIAxes_theta4dd.Position = [1708 265 205 166];

            % Create UIAxes_d3dd
            app.UIAxes_d3dd = uiaxes(app.UIFigure);
            title(app.UIAxes_d3dd, 'd3 dot dot')
            xlabel(app.UIAxes_d3dd, 'X')
            ylabel(app.UIAxes_d3dd, 'Y')
            zlabel(app.UIAxes_d3dd, 'Z')
            app.UIAxes_d3dd.Position = [1708 457 205 166];

            % Create UIAxes_theta2dd
            app.UIAxes_theta2dd = uiaxes(app.UIFigure);
            title(app.UIAxes_theta2dd, 'theta2 dot dot dot')
            xlabel(app.UIAxes_theta2dd, 'X')
            ylabel(app.UIAxes_theta2dd, 'Y')
            zlabel(app.UIAxes_theta2dd, 'Z')
            app.UIAxes_theta2dd.Position = [1708 648 205 166];

            % Create UIAxes_theta1dd
            app.UIAxes_theta1dd = uiaxes(app.UIFigure);
            title(app.UIAxes_theta1dd, 'theta1 dot dot')
            xlabel(app.UIAxes_theta1dd, 'X')
            ylabel(app.UIAxes_theta1dd, 'Y')
            zlabel(app.UIAxes_theta1dd, 'Z')
            app.UIAxes_theta1dd.Position = [1708 839 205 166];

            % Create UIAxes_eev
            app.UIAxes_eev = uiaxes(app.UIFigure);
            title(app.UIAxes_eev, 'ee v(t)')
            xlabel(app.UIAxes_eev, 'X')
            ylabel(app.UIAxes_eev, 'Y')
            zlabel(app.UIAxes_eev, 'Z')
            app.UIAxes_eev.Position = [1599 63 205 166];

            % Create UIAxes_eea
            app.UIAxes_eea = uiaxes(app.UIFigure);
            title(app.UIAxes_eea, 'ee a(t)')
            xlabel(app.UIAxes_eea, 'X')
            ylabel(app.UIAxes_eea, 'Y')
            zlabel(app.UIAxes_eea, 'Z')
            app.UIAxes_eea.Position = [1624 55 205 166];

            % Create UIAxes_eep
            app.UIAxes_eep = uiaxes(app.UIFigure);
            title(app.UIAxes_eep, 'ee p(t)')
            xlabel(app.UIAxes_eep, 'X')
            ylabel(app.UIAxes_eep, 'Y')
            zlabel(app.UIAxes_eep, 'Z')
            app.UIAxes_eep.Position = [1294 39 624 227];

            % Create WORKSPACEPanel
            app.WORKSPACEPanel = uipanel(app.UIFigure);
            app.WORKSPACEPanel.TitlePosition = 'centertop';
            app.WORKSPACEPanel.Title = 'WORKSPACE';
            app.WORKSPACEPanel.BackgroundColor = [1 1 1];
            app.WORKSPACEPanel.FontWeight = 'bold';
            app.WORKSPACEPanel.FontSize = 15;
            app.WORKSPACEPanel.Position = [25 855 169 150];

            % Create WS_OFFButton
            app.WS_OFFButton = uibutton(app.WORKSPACEPanel, 'push');
            app.WS_OFFButton.ButtonPushedFcn = createCallbackFcn(app, @WS_OFFButtonPushed, true);
            app.WS_OFFButton.Position = [61 40 48 22];
            app.WS_OFFButton.Text = 'OFF';

            % Create WS_ONButton
            app.WS_ONButton = uibutton(app.WORKSPACEPanel, 'push');
            app.WS_ONButton.ButtonPushedFcn = createCallbackFcn(app, @WS_ONButtonPushed, true);
            app.WS_ONButton.Position = [61 72 47 22];
            app.WS_ONButton.Text = 'ON';

            % Create ENDEFFECTORPanel
            app.ENDEFFECTORPanel = uipanel(app.UIFigure);
            app.ENDEFFECTORPanel.TitlePosition = 'centertop';
            app.ENDEFFECTORPanel.Title = 'END EFFECTOR';
            app.ENDEFFECTORPanel.BackgroundColor = [1 1 1];
            app.ENDEFFECTORPanel.FontWeight = 'bold';
            app.ENDEFFECTORPanel.FontSize = 15;
            app.ENDEFFECTORPanel.Position = [943 37 280 325];

            % Create Xee_label
            app.Xee_label = uilabel(app.ENDEFFECTORPanel);
            app.Xee_label.FontSize = 15;
            app.Xee_label.FontWeight = 'bold';
            app.Xee_label.Position = [66 226 68 33];
            app.Xee_label.Text = 'X';

            % Create Yee_label
            app.Yee_label = uilabel(app.ENDEFFECTORPanel);
            app.Yee_label.FontSize = 15;
            app.Yee_label.FontWeight = 'bold';
            app.Yee_label.Position = [66 185 68 33];
            app.Yee_label.Text = 'Y';

            % Create Zee_label
            app.Zee_label = uilabel(app.ENDEFFECTORPanel);
            app.Zee_label.FontSize = 15;
            app.Zee_label.FontWeight = 'bold';
            app.Zee_label.Position = [66 143 68 33];
            app.Zee_label.Text = 'Z';

            % Create Xee
            app.Xee = uieditfield(app.ENDEFFECTORPanel, 'numeric');
            app.Xee.Limits = [-6 6];
            app.Xee.ValueDisplayFormat = '%.3f';
            app.Xee.Editable = 'off';
            app.Xee.FontSize = 15;
            app.Xee.Position = [90 226 143 33];

            % Create Yee
            app.Yee = uieditfield(app.ENDEFFECTORPanel, 'numeric');
            app.Yee.Limits = [-6 6];
            app.Yee.ValueDisplayFormat = '%.3f';
            app.Yee.Editable = 'off';
            app.Yee.FontSize = 15;
            app.Yee.Position = [90 185 143 33];

            % Create Zee
            app.Zee = uieditfield(app.ENDEFFECTORPanel, 'numeric');
            app.Zee.ValueDisplayFormat = '%.3f';
            app.Zee.Editable = 'off';
            app.Zee.FontSize = 15;
            app.Zee.Position = [90 143 143 33];
            app.Zee.Value = 2;

            % Create Yawee_label
            app.Yawee_label = uilabel(app.ENDEFFECTORPanel);
            app.Yawee_label.FontSize = 15;
            app.Yawee_label.FontWeight = 'bold';
            app.Yawee_label.Position = [56 100 71 33];
            app.Yawee_label.Text = 'Yaw';

            % Create Yawee
            app.Yawee = uieditfield(app.ENDEFFECTORPanel, 'numeric');
            app.Yawee.ValueDisplayFormat = '%.3f';
            app.Yawee.Editable = 'off';
            app.Yawee.FontSize = 15;
            app.Yawee.Position = [90 100 143 33];

            % Create FORDWARDKINEMATICSPanel
            app.FORDWARDKINEMATICSPanel = uipanel(app.UIFigure);
            app.FORDWARDKINEMATICSPanel.TitlePosition = 'centertop';
            app.FORDWARDKINEMATICSPanel.Title = 'FORDWARD KINEMATICS';
            app.FORDWARDKINEMATICSPanel.BackgroundColor = [1 1 1];
            app.FORDWARDKINEMATICSPanel.FontWeight = 'bold';
            app.FORDWARDKINEMATICSPanel.FontSize = 15;
            app.FORDWARDKINEMATICSPanel.Position = [25 429 311 366];

            % Create EditField_theta4
            app.EditField_theta4 = uieditfield(app.FORDWARDKINEMATICSPanel, 'numeric');
            app.EditField_theta4.ValueDisplayFormat = '%.3f';
            app.EditField_theta4.ValueChangedFcn = createCallbackFcn(app, @EditField_theta4ValueChanged, true);
            app.EditField_theta4.Position = [260 63 42 22];

            % Create EditField_d3
            app.EditField_d3 = uieditfield(app.FORDWARDKINEMATICSPanel, 'numeric');
            app.EditField_d3.ValueDisplayFormat = '%.3f';
            app.EditField_d3.ValueChangedFcn = createCallbackFcn(app, @EditField_d3ValueChanged, true);
            app.EditField_d3.Position = [260 138 42 22];
            app.EditField_d3.Value = 1;

            % Create EditField_theta2
            app.EditField_theta2 = uieditfield(app.FORDWARDKINEMATICSPanel, 'numeric');
            app.EditField_theta2.ValueDisplayFormat = '%.3f';
            app.EditField_theta2.ValueChangedFcn = createCallbackFcn(app, @EditField_theta2ValueChanged, true);
            app.EditField_theta2.Position = [260 213 42 22];

            % Create EditField_theta1
            app.EditField_theta1 = uieditfield(app.FORDWARDKINEMATICSPanel, 'numeric');
            app.EditField_theta1.ValueDisplayFormat = '%.3f';
            app.EditField_theta1.ValueChangedFcn = createCallbackFcn(app, @EditField_theta1ValueChanged, true);
            app.EditField_theta1.Position = [260 287 42 22];

            % Create theta1Slider
            app.theta1Slider = uislider(app.FORDWARDKINEMATICSPanel);
            app.theta1Slider.Limits = [-2.18166156499291 2.18166156499291];
            app.theta1Slider.ValueChangingFcn = createCallbackFcn(app, @theta1SliderValueChanging, true);
            app.theta1Slider.Position = [70 296 150 3];

            % Create theta1SliderLabel
            app.theta1SliderLabel = uilabel(app.FORDWARDKINEMATICSPanel);
            app.theta1SliderLabel.HorizontalAlignment = 'right';
            app.theta1SliderLabel.Position = [10 287 39 22];
            app.theta1SliderLabel.Text = 'theta1';

            % Create theta2Slider
            app.theta2Slider = uislider(app.FORDWARDKINEMATICSPanel);
            app.theta2Slider.Limits = [-2.53072741539178 2.53072741539178];
            app.theta2Slider.ValueChangingFcn = createCallbackFcn(app, @theta2SliderValueChanging, true);
            app.theta2Slider.Position = [70 223 150 3];

            % Create theta2SliderLabel
            app.theta2SliderLabel = uilabel(app.FORDWARDKINEMATICSPanel);
            app.theta2SliderLabel.HorizontalAlignment = 'right';
            app.theta2SliderLabel.Position = [10 214 39 22];
            app.theta2SliderLabel.Text = 'theta2';

            % Create d3Slider
            app.d3Slider = uislider(app.FORDWARDKINEMATICSPanel);
            app.d3Slider.Limits = [0 1.6];
            app.d3Slider.ValueChangingFcn = createCallbackFcn(app, @d3SliderValueChanging, true);
            app.d3Slider.Position = [70 147 150 3];
            app.d3Slider.Value = 1;

            % Create d3SliderLabel
            app.d3SliderLabel = uilabel(app.FORDWARDKINEMATICSPanel);
            app.d3SliderLabel.HorizontalAlignment = 'right';
            app.d3SliderLabel.Position = [10 138 25 22];
            app.d3SliderLabel.Text = 'd3';

            % Create theta4Slider
            app.theta4Slider = uislider(app.FORDWARDKINEMATICSPanel);
            app.theta4Slider.Limits = [-3.14159265358979 3.14159265358979];
            app.theta4Slider.ValueChangingFcn = createCallbackFcn(app, @theta4SliderValueChanging, true);
            app.theta4Slider.Position = [70 72 150 3];

            % Create theta4SliderLabel
            app.theta4SliderLabel = uilabel(app.FORDWARDKINEMATICSPanel);
            app.theta4SliderLabel.HorizontalAlignment = 'right';
            app.theta4SliderLabel.Position = [10 63 39 22];
            app.theta4SliderLabel.Text = 'theta4';

            % Create RESETButton
            app.RESETButton = uibutton(app.FORDWARDKINEMATICSPanel, 'push');
            app.RESETButton.ButtonPushedFcn = createCallbackFcn(app, @RESETButtonPushed, true);
            app.RESETButton.FontSize = 14;
            app.RESETButton.FontWeight = 'bold';
            app.RESETButton.Position = [123 9 65 24];
            app.RESETButton.Text = 'RESET';

            % Create INVERSEKINEMATICSPanel
            app.INVERSEKINEMATICSPanel = uipanel(app.UIFigure);
            app.INVERSEKINEMATICSPanel.TitlePosition = 'centertop';
            app.INVERSEKINEMATICSPanel.Title = 'INVERSE KINEMATICS';
            app.INVERSEKINEMATICSPanel.BackgroundColor = [1 1 1];
            app.INVERSEKINEMATICSPanel.FontWeight = 'bold';
            app.INVERSEKINEMATICSPanel.FontSize = 15;
            app.INVERSEKINEMATICSPanel.Position = [354 429 200 366];

            % Create Yaw_inverseLabel
            app.Yaw_inverseLabel = uilabel(app.INVERSEKINEMATICSPanel);
            app.Yaw_inverseLabel.Position = [33 63 28 22];
            app.Yaw_inverseLabel.Text = 'Yaw';

            % Create Z_inverseLabel
            app.Z_inverseLabel = uilabel(app.INVERSEKINEMATICSPanel);
            app.Z_inverseLabel.Position = [36 138 25 22];
            app.Z_inverseLabel.Text = 'Z';

            % Create Y_inverseLabel
            app.Y_inverseLabel = uilabel(app.INVERSEKINEMATICSPanel);
            app.Y_inverseLabel.Position = [36 213 25 22];
            app.Y_inverseLabel.Text = 'Y';

            % Create X_inverseLabel
            app.X_inverseLabel = uilabel(app.INVERSEKINEMATICSPanel);
            app.X_inverseLabel.Position = [36 287 25 22];
            app.X_inverseLabel.Text = 'X';

            % Create inverse_Yaw
            app.inverse_Yaw = uieditfield(app.INVERSEKINEMATICSPanel, 'numeric');
            app.inverse_Yaw.ValueDisplayFormat = '%.3f';
            app.inverse_Yaw.ValueChangedFcn = createCallbackFcn(app, @inverse_YawValueChanged, true);
            app.inverse_Yaw.Position = [60 63 100 22];

            % Create inverse_X
            app.inverse_X = uieditfield(app.INVERSEKINEMATICSPanel, 'numeric');
            app.inverse_X.Limits = [-6 6];
            app.inverse_X.ValueDisplayFormat = '%.3f';
            app.inverse_X.ValueChangedFcn = createCallbackFcn(app, @inverse_XValueChanged, true);
            app.inverse_X.Position = [60 287 100 22];

            % Create inverse_Y
            app.inverse_Y = uieditfield(app.INVERSEKINEMATICSPanel, 'numeric');
            app.inverse_Y.Limits = [-6 6];
            app.inverse_Y.ValueDisplayFormat = '%.3f';
            app.inverse_Y.ValueChangedFcn = createCallbackFcn(app, @inverse_YValueChanged, true);
            app.inverse_Y.Position = [60 213 100 22];

            % Create inverse_Z
            app.inverse_Z = uieditfield(app.INVERSEKINEMATICSPanel, 'numeric');
            app.inverse_Z.ValueDisplayFormat = '%.3f';
            app.inverse_Z.ValueChangedFcn = createCallbackFcn(app, @inverse_ZValueChanged, true);
            app.inverse_Z.Position = [60 138 100 22];

            % Create GOTOButton
            app.GOTOButton = uibutton(app.INVERSEKINEMATICSPanel, 'push');
            app.GOTOButton.ButtonPushedFcn = createCallbackFcn(app, @GOTOButtonPushed, true);
            app.GOTOButton.FontSize = 14;
            app.GOTOButton.FontWeight = 'bold';
            app.GOTOButton.Position = [67 9 65 24];
            app.GOTOButton.Text = 'GO TO';

            % Create TRAJECTORYPLANNINGPanel
            app.TRAJECTORYPLANNINGPanel = uipanel(app.UIFigure);
            app.TRAJECTORYPLANNINGPanel.TitlePosition = 'centertop';
            app.TRAJECTORYPLANNINGPanel.Title = 'TRAJECTORY PLANNING';
            app.TRAJECTORYPLANNINGPanel.BackgroundColor = [1 1 1];
            app.TRAJECTORYPLANNINGPanel.FontWeight = 'bold';
            app.TRAJECTORYPLANNINGPanel.FontSize = 15;
            app.TRAJECTORYPLANNINGPanel.Position = [25 37 821 325];

            % Create MaxVelEditFieldLabel
            app.MaxVelEditFieldLabel = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.MaxVelEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxVelEditFieldLabel.Position = [84 258 52 22];
            app.MaxVelEditFieldLabel.Text = ' Max Vel';

            % Create MaxVel_ax1
            app.MaxVel_ax1 = uieditfield(app.TRAJECTORYPLANNINGPanel, 'numeric');
            app.MaxVel_ax1.Position = [153 258 56 22];
            app.MaxVel_ax1.Value = 11.5191730631626;

            % Create MaxAccEditFieldLabel
            app.MaxAccEditFieldLabel = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.MaxAccEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxAccEditFieldLabel.Position = [85 228 51 22];
            app.MaxAccEditFieldLabel.Text = 'Max Acc';

            % Create MaxAcc_ax1
            app.MaxAcc_ax1 = uieditfield(app.TRAJECTORYPLANNINGPanel, 'numeric');
            app.MaxAcc_ax1.Position = [153 228 56 22];
            app.MaxAcc_ax1.Value = 115.2;

            % Create Panel
            app.Panel = uipanel(app.TRAJECTORYPLANNINGPanel);
            app.Panel.BackgroundColor = [1 1 1];
            app.Panel.Position = [469 16 318 264];

            % Create A_Y
            app.A_Y = uieditfield(app.Panel, 'numeric');
            app.A_Y.HorizontalAlignment = 'center';
            app.A_Y.FontSize = 15;
            app.A_Y.Position = [121 152 54 29];
            app.A_Y.Value = -2;

            % Create XEditFieldLabel_2
            app.XEditFieldLabel_2 = uilabel(app.Panel);
            app.XEditFieldLabel_2.HorizontalAlignment = 'center';
            app.XEditFieldLabel_2.FontSize = 15;
            app.XEditFieldLabel_2.Position = [121 181 54 29];
            app.XEditFieldLabel_2.Text = 'Y';

            % Create ALabel
            app.ALabel = uilabel(app.Panel);
            app.ALabel.HorizontalAlignment = 'center';
            app.ALabel.FontSize = 15;
            app.ALabel.FontWeight = 'bold';
            app.ALabel.Position = [7 149 47 35];
            app.ALabel.Text = 'A';

            % Create BLabel
            app.BLabel = uilabel(app.Panel);
            app.BLabel.HorizontalAlignment = 'center';
            app.BLabel.FontSize = 15;
            app.BLabel.FontWeight = 'bold';
            app.BLabel.Position = [7 106 47 35];
            app.BLabel.Text = 'B';

            % Create XEditFieldLabel_17
            app.XEditFieldLabel_17 = uilabel(app.Panel);
            app.XEditFieldLabel_17.HorizontalAlignment = 'center';
            app.XEditFieldLabel_17.FontSize = 15;
            app.XEditFieldLabel_17.Position = [58 181 54 29];
            app.XEditFieldLabel_17.Text = 'X';

            % Create A_X
            app.A_X = uieditfield(app.Panel, 'numeric');
            app.A_X.HorizontalAlignment = 'center';
            app.A_X.FontSize = 15;
            app.A_X.Position = [58 152 54 29];

            % Create XEditFieldLabel_18
            app.XEditFieldLabel_18 = uilabel(app.Panel);
            app.XEditFieldLabel_18.HorizontalAlignment = 'center';
            app.XEditFieldLabel_18.FontSize = 15;
            app.XEditFieldLabel_18.Position = [184 181 54 29];
            app.XEditFieldLabel_18.Text = 'Z';

            % Create A_Z
            app.A_Z = uieditfield(app.Panel, 'numeric');
            app.A_Z.HorizontalAlignment = 'center';
            app.A_Z.FontSize = 15;
            app.A_Z.Position = [184 152 54 29];
            app.A_Z.Value = 0.2;

            % Create XEditFieldLabel_19
            app.XEditFieldLabel_19 = uilabel(app.Panel);
            app.XEditFieldLabel_19.HorizontalAlignment = 'center';
            app.XEditFieldLabel_19.FontSize = 15;
            app.XEditFieldLabel_19.Position = [246 181 54 29];
            app.XEditFieldLabel_19.Text = 'Yaw';

            % Create A_Yaw
            app.A_Yaw = uieditfield(app.Panel, 'numeric');
            app.A_Yaw.HorizontalAlignment = 'center';
            app.A_Yaw.FontSize = 15;
            app.A_Yaw.Position = [246 152 54 29];

            % Create B_X
            app.B_X = uieditfield(app.Panel, 'numeric');
            app.B_X.HorizontalAlignment = 'center';
            app.B_X.FontSize = 15;
            app.B_X.Position = [58 109 54 29];

            % Create B_Y
            app.B_Y = uieditfield(app.Panel, 'numeric');
            app.B_Y.HorizontalAlignment = 'center';
            app.B_Y.FontSize = 15;
            app.B_Y.Position = [121 109 54 29];
            app.B_Y.Value = 2;

            % Create B_Z
            app.B_Z = uieditfield(app.Panel, 'numeric');
            app.B_Z.HorizontalAlignment = 'center';
            app.B_Z.FontSize = 15;
            app.B_Z.Position = [184 109 54 29];
            app.B_Z.Value = 1.4;

            % Create B_Yaw
            app.B_Yaw = uieditfield(app.Panel, 'numeric');
            app.B_Yaw.HorizontalAlignment = 'center';
            app.B_Yaw.FontSize = 15;
            app.B_Yaw.Position = [246 109 54 29];
            app.B_Yaw.Value = 3.14159265358979;

            % Create GOButton
            app.GOButton = uibutton(app.Panel, 'push');
            app.GOButton.ButtonPushedFcn = createCallbackFcn(app, @GOButtonPushed, true);
            app.GOButton.FontSize = 15;
            app.GOButton.Position = [109 51 119 29];
            app.GOButton.Text = 'GO';

            % Create MaxAccEditField_2Label
            app.MaxAccEditField_2Label = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.MaxAccEditField_2Label.HorizontalAlignment = 'right';
            app.MaxAccEditField_2Label.Position = [85 162 51 22];
            app.MaxAccEditField_2Label.Text = 'Max Acc';

            % Create MaxAcc_ax2
            app.MaxAcc_ax2 = uieditfield(app.TRAJECTORYPLANNINGPanel, 'numeric');
            app.MaxAcc_ax2.Position = [153 162 56 22];
            app.MaxAcc_ax2.Value = 115.2;

            % Create MaxVelEditField_2Label
            app.MaxVelEditField_2Label = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.MaxVelEditField_2Label.HorizontalAlignment = 'right';
            app.MaxVelEditField_2Label.Position = [84 192 52 22];
            app.MaxVelEditField_2Label.Text = ' Max Vel';

            % Create MaxVel_ax2
            app.MaxVel_ax2 = uieditfield(app.TRAJECTORYPLANNINGPanel, 'numeric');
            app.MaxVel_ax2.Position = [153 192 56 22];
            app.MaxVel_ax2.Value = 11.5191730631626;

            % Create MaxAccEditField_3Label
            app.MaxAccEditField_3Label = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.MaxAccEditField_3Label.HorizontalAlignment = 'right';
            app.MaxAccEditField_3Label.Position = [85 93 51 22];
            app.MaxAccEditField_3Label.Text = 'Max Acc';

            % Create MaxAcc_ax3
            app.MaxAcc_ax3 = uieditfield(app.TRAJECTORYPLANNINGPanel, 'numeric');
            app.MaxAcc_ax3.Position = [153 93 56 22];
            app.MaxAcc_ax3.Value = 112;

            % Create MaxVelEditField_3Label
            app.MaxVelEditField_3Label = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.MaxVelEditField_3Label.HorizontalAlignment = 'right';
            app.MaxVelEditField_3Label.Position = [84 123 52 22];
            app.MaxVelEditField_3Label.Text = ' Max Vel';

            % Create MaxVel_ax3
            app.MaxVel_ax3 = uieditfield(app.TRAJECTORYPLANNINGPanel, 'numeric');
            app.MaxVel_ax3.Position = [153 123 56 22];
            app.MaxVel_ax3.Value = 11.2;

            % Create MaxAccEditField_4Label
            app.MaxAccEditField_4Label = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.MaxAccEditField_4Label.HorizontalAlignment = 'right';
            app.MaxAccEditField_4Label.Position = [85 26 51 22];
            app.MaxAccEditField_4Label.Text = 'Max Acc';

            % Create MaxAcc_ax4
            app.MaxAcc_ax4 = uieditfield(app.TRAJECTORYPLANNINGPanel, 'numeric');
            app.MaxAcc_ax4.Position = [153 26 56 22];
            app.MaxAcc_ax4.Value = 261.8;

            % Create MaxVelEditField_4Label
            app.MaxVelEditField_4Label = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.MaxVelEditField_4Label.HorizontalAlignment = 'right';
            app.MaxVelEditField_4Label.Position = [84 56 52 22];
            app.MaxVelEditField_4Label.Text = ' Max Vel';

            % Create MaxVel_ax4
            app.MaxVel_ax4 = uieditfield(app.TRAJECTORYPLANNINGPanel, 'numeric');
            app.MaxVel_ax4.Position = [153 56 56 22];
            app.MaxVel_ax4.Value = 26.1799387799149;

            % Create AXIS1Label
            app.AXIS1Label = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.AXIS1Label.FontWeight = 'bold';
            app.AXIS1Label.Position = [36 233 44 41];
            app.AXIS1Label.Text = 'AXIS 1';

            % Create AXIS2Label
            app.AXIS2Label = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.AXIS2Label.FontWeight = 'bold';
            app.AXIS2Label.Position = [36 168 44 41];
            app.AXIS2Label.Text = 'AXIS 2';

            % Create AXIS3Label
            app.AXIS3Label = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.AXIS3Label.FontWeight = 'bold';
            app.AXIS3Label.Position = [36 99 44 41];
            app.AXIS3Label.Text = 'AXIS 3';

            % Create AXIS4Label
            app.AXIS4Label = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.AXIS4Label.FontWeight = 'bold';
            app.AXIS4Label.Position = [36 36 44 41];
            app.AXIS4Label.Text = 'AXIS 4';

            % Create radsLabel
            app.radsLabel = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.radsLabel.Position = [217 257 32 22];
            app.radsLabel.Text = 'rad/s';

            % Create radsLabel_2
            app.radsLabel_2 = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.radsLabel_2.Position = [217 191 32 22];
            app.radsLabel_2.Text = 'rad/s';

            % Create radsLabel_3
            app.radsLabel_3 = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.radsLabel_3.Position = [217 56 32 22];
            app.radsLabel_3.Text = 'rad/s';

            % Create rad2sLabel
            app.rad2sLabel = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.rad2sLabel.Position = [217 228 44 22];
            app.rad2sLabel.Text = 'rad^2/s';

            % Create rad2sLabel_2
            app.rad2sLabel_2 = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.rad2sLabel_2.Position = [217 162 44 22];
            app.rad2sLabel_2.Text = 'rad^2/s';

            % Create rad2sLabel_3
            app.rad2sLabel_3 = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.rad2sLabel_3.Position = [217 26 44 22];
            app.rad2sLabel_3.Text = 'rad^2/s';

            % Create dm2sLabel
            app.dm2sLabel = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.dm2sLabel.Position = [217 93 44 22];
            app.dm2sLabel.Text = 'dm^2/s';

            % Create dmsLabel
            app.dmsLabel = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.dmsLabel.Position = [217 123 32 22];
            app.dmsLabel.Text = 'dm/s';

            % Create sLabel
            app.sLabel = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.sLabel.Position = [400 242 25 22];
            app.sLabel.Text = 's';

            % Create timeEditFieldLabel
            app.timeEditFieldLabel = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.timeEditFieldLabel.HorizontalAlignment = 'right';
            app.timeEditFieldLabel.Position = [291 243 28 22];
            app.timeEditFieldLabel.Text = 'time';

            % Create timeAxis1
            app.timeAxis1 = uieditfield(app.TRAJECTORYPLANNINGPanel, 'numeric');
            app.timeAxis1.Editable = 'off';
            app.timeAxis1.Position = [336 243 56 22];

            % Create sLabel_2
            app.sLabel_2 = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.sLabel_2.Position = [398 176 25 22];
            app.sLabel_2.Text = 's';

            % Create timeEditField_2Label
            app.timeEditField_2Label = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.timeEditField_2Label.HorizontalAlignment = 'right';
            app.timeEditField_2Label.Position = [289 177 28 22];
            app.timeEditField_2Label.Text = 'time';

            % Create timeAxis2
            app.timeAxis2 = uieditfield(app.TRAJECTORYPLANNINGPanel, 'numeric');
            app.timeAxis2.Editable = 'off';
            app.timeAxis2.Position = [334 177 56 22];

            % Create sLabel_3
            app.sLabel_3 = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.sLabel_3.Position = [398 107 25 22];
            app.sLabel_3.Text = 's';

            % Create timeEditField_3Label
            app.timeEditField_3Label = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.timeEditField_3Label.HorizontalAlignment = 'right';
            app.timeEditField_3Label.Position = [289 108 28 22];
            app.timeEditField_3Label.Text = 'time';

            % Create timeAxis3
            app.timeAxis3 = uieditfield(app.TRAJECTORYPLANNINGPanel, 'numeric');
            app.timeAxis3.Editable = 'off';
            app.timeAxis3.Position = [334 108 56 22];

            % Create sLabel_4
            app.sLabel_4 = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.sLabel_4.Position = [398 43 25 22];
            app.sLabel_4.Text = 's';

            % Create timeEditField_4Label
            app.timeEditField_4Label = uilabel(app.TRAJECTORYPLANNINGPanel);
            app.timeEditField_4Label.HorizontalAlignment = 'right';
            app.timeEditField_4Label.Position = [289 44 28 22];
            app.timeEditField_4Label.Text = 'time';

            % Create timeAxis4
            app.timeAxis4 = uieditfield(app.TRAJECTORYPLANNINGPanel, 'numeric');
            app.timeAxis4.Editable = 'off';
            app.timeAxis4.Position = [334 44 56 22];

            % Create COORDINATEPanel
            app.COORDINATEPanel = uipanel(app.UIFigure);
            app.COORDINATEPanel.TitlePosition = 'centertop';
            app.COORDINATEPanel.Title = 'COORDINATE';
            app.COORDINATEPanel.BackgroundColor = [1 1 1];
            app.COORDINATEPanel.FontWeight = 'bold';
            app.COORDINATEPanel.FontSize = 15;
            app.COORDINATEPanel.Position = [233 855 320 150];

            % Create Coor0CheckBox
            app.Coor0CheckBox = uicheckbox(app.COORDINATEPanel);
            app.Coor0CheckBox.ValueChangedFcn = createCallbackFcn(app, @Coor0CheckBoxValueChanged, true);
            app.Coor0CheckBox.Text = 'Coor - 0';
            app.Coor0CheckBox.Position = [42 95 88 22];

            % Create Coor1CheckBox
            app.Coor1CheckBox = uicheckbox(app.COORDINATEPanel);
            app.Coor1CheckBox.ValueChangedFcn = createCallbackFcn(app, @Coor1CheckBoxValueChanged, true);
            app.Coor1CheckBox.Text = 'Coor - 1';
            app.Coor1CheckBox.Position = [42 74 88 22];

            % Create Coor2CheckBox
            app.Coor2CheckBox = uicheckbox(app.COORDINATEPanel);
            app.Coor2CheckBox.ValueChangedFcn = createCallbackFcn(app, @Coor2CheckBoxValueChanged, true);
            app.Coor2CheckBox.Text = 'Coor - 2';
            app.Coor2CheckBox.Position = [42 53 88 22];

            % Create Coor3CheckBox
            app.Coor3CheckBox = uicheckbox(app.COORDINATEPanel);
            app.Coor3CheckBox.ValueChangedFcn = createCallbackFcn(app, @Coor3CheckBoxValueChanged, true);
            app.Coor3CheckBox.Text = 'Coor - 3';
            app.Coor3CheckBox.Position = [42 31 88 22];

            % Create Coor4CheckBox
            app.Coor4CheckBox = uicheckbox(app.COORDINATEPanel);
            app.Coor4CheckBox.ValueChangedFcn = createCallbackFcn(app, @Coor4CheckBoxValueChanged, true);
            app.Coor4CheckBox.Text = 'Coor - 4';
            app.Coor4CheckBox.Position = [42 10 88 22];

            % Create VIEWALLButton
            app.VIEWALLButton = uibutton(app.COORDINATEPanel, 'push');
            app.VIEWALLButton.ButtonPushedFcn = createCallbackFcn(app, @VIEWALLButtonPushed, true);
            app.VIEWALLButton.Position = [176 51 70 23];
            app.VIEWALLButton.Text = 'VIEW ALL';

            % Create CLEARButton
            app.CLEARButton = uibutton(app.COORDINATEPanel, 'push');
            app.CLEARButton.ButtonPushedFcn = createCallbackFcn(app, @CLEARButtonPushed, true);
            app.CLEARButton.Position = [176 20 70 23];
            app.CLEARButton.Text = 'CLEAR';

            % Create OpacityEditFieldLabel
            app.OpacityEditFieldLabel = uilabel(app.COORDINATEPanel);
            app.OpacityEditFieldLabel.HorizontalAlignment = 'right';
            app.OpacityEditFieldLabel.Position = [169 91 46 22];
            app.OpacityEditFieldLabel.Text = 'Opacity';

            % Create OpacityEditField
            app.OpacityEditField = uieditfield(app.COORDINATEPanel, 'numeric');
            app.OpacityEditField.ValueChangedFcn = createCallbackFcn(app, @OpacityEditFieldValueChanged, true);
            app.OpacityEditField.Position = [222 91 31 22];
            app.OpacityEditField.Value = 1;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = scara_app_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end