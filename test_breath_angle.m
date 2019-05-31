function [ in_is_true ] = test_breath_angle( angle_to_test, breath_angle, breath_range )
%test_breath_angle tests to see if the first argument is a radians angle within the breath_range of breath_angel
%    test_breath_angle( angle_to_test, breath_angle, breath_range ) returns
%    true (1) if angle_to_test lies between breath_angle-breath_range and
%    breath_angle+breath_range.  If breath_angle-breath_range < 0 or
%    breath_angle+breath_range is greater than 2*pi then special wrap
%    around conditions are included in the tests
if test_0_to_2pi(angle_to_test)&&test_0_to_2pi(breath_angle)&&test_0_to_2pi(breath_range)
    
    in_is_true = 0;  % start false and set to true when find in range
    
    if breath_angle-breath_range<=0+10*eps
        % in this case check between wrap around value and 2*pi
        wrap_around_min=mod(breath_angle-breath_range, 2*pi);
        if angle_to_test>=wrap_around_min-10*eps
            in_is_true = 1;
            return
        end
        % check to see if in 0 to max part of range
        if angle_to_test<=breath_angle+breath_range+10*eps
            in_is_true = 1;
            return
        end
        return % with false
    end
    if breath_angle+breath_range>= 2*pi-10*eps
        % in this case check wrap around between 0 and value
        wrap_around_max = mod(breath_angle+breath_range,2*pi);
        if angle_to_test<=wrap_around_max+10*eps
            in_is_true = 1;
            return
        end
    end
    %ordinary case follows
    if ((breath_angle-breath_range)<=10*eps+angle_to_test) && (angle_to_test<=10*eps+(breath_angle+breath_range))
        in_is_true = 1;
        return
    end
    return % returns false otherwise
end
disp([num2str(angle_to_test-2*pi)])
error(['one of the input angles to test_breath_angle(' num2str(angle_to_test) ', ' num2str(breath_angle) ', ' num2str(breath_range) ') is not in interval [0, 2pi]'])
end

