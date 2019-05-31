function [ in_range ] = test_0_to_2pi( angle )
% test_0_to_2pi(angle) returns true (1) if in range [0,2pi] and false (0)
% otherwise.
in_range = 0;
if (angle+10*eps>=0) && (angle<=2*pi+10*eps)
    in_range = 1;
end
return

end

