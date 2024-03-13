function inv = inverse(a, alpha, d, theta, end_effector)
   
    c2 = (end_effector(1)^2 + end_effector(2)^2 - a(2)^2 - a(3)^2)/(2*a(2)*a(3));  %costheta2
    if c2>1
         warndlg('Ngoai workspace', 'Warning');
         main_gui.reset_Callback(hObject, eventdata, handles);
         return;
    end
    s2 = sqrt(1 - c2^2);  %sintheta2
    inv(2) = atan2(s2, c2); %theta2
    s1 = ((a(2)+a(3)*c2)*end_effector(2) - a(3)*s2*end_effector(1))/(end_effector(1)^2 + end_effector(2)^2); %sintheta1
    c1 = ((a(2)+a(3)*c2)*end_effector(1) + a(3)*s2*end_effector(2))/(end_effector(1)^2 + end_effector(2)^2);
    inv(1) = atan2(s1, c1); %theta1
    inv(4) = wrapToPi((inv(1) + inv(2) - end_effector(6))); %theta4, do no di tu -360 -> 360 nen ta gioi han lai tu -pi -> pi bang ham wraptoPi
    inv(3) = +d(1) - end_effector(3);   %d3 = z - d1
    
    %%
    %neu nam ngoai workspace thi cho gia tri bang gia tri max tai workspace
    %xet d3 co nam ngoai workspace hay khong
    if (inv(3) > 2.1)
        inv(3) = 2.1;
        warndlg('Check Inverse d3 > 2.1', 'Warning');
    elseif (inv(3) < 0)
        inv(3) = 0;
        warndlg('Check Inverse d3 < 0', 'Warning');
    end
           
    %xet theta1 co nam ngoai workspace ko
    if (inv(1) > deg2rad(148))
        inv(1) = deg2rad(148);
        warndlg('Check Inverse theta1 > 148', 'Warning');
    elseif (inv(1) < deg2rad(-148))
        inv(1) = deg2rad(-148);
        warndlg('Check Inverse theta1 < 148', 'Warning');
    end

    %xet theta2 co nam ngoai workspace ko
    if (inv(2) > deg2rad(150))
        inv(2) = deg2rad(150);
        warndlg('Check Inverse theta2 > 150', 'Warning');
    elseif (inv(2) < deg2rad(-150))
        inv(2) = deg2rad(-150);
        warndlg('Check Inverse theta2 < 150', 'Warning');
    end

end