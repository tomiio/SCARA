function out = create_trajectory(q1, q2, a_max, v_max)
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
    %total_step = floor(t5 / 0.01);
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
