function a_out = sat(a_in)
% sat is an approximation of the saturation of a function between the values [-1,1]

 a_out = (atan(5.*a_in) - atan(-5.*a_in))/pi;
 
end
 