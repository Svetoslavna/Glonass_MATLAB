function [output2] = StepMove (in, step)
%%метод рунге-кутта
    temp = SmallStepMove(in, step);
    temp11 = struct2array(temp);
    temp1.xyz = temp11(1:3); temp1.v_xyz = temp11(4:6);temp1.a_xyz = temp11(7:9);
    in1 = sum_struct(in,mult(temp1,0.5));
    temp22 = struct2array(SmallStepMove(in1,step));temp2.xyz = temp22(1:3); temp2.v_xyz = temp22(4:6);temp2.a_xyz = temp22(7:9);
    in2 = sum_struct(in,mult(temp2,0.5));
    temp33 = struct2array(SmallStepMove(in2,step));temp3.xyz = temp33(1:3); temp3.v_xyz = temp33(4:6);temp3.a_xyz = temp33(7:9);
    in3 = sum_struct(in,temp3);
    temp44 = struct2array(SmallStepMove(in3,step));temp4.xyz = temp44(1:3); temp4.v_xyz = temp44(4:6);temp4.a_xyz = temp44(7:9);
    output2 = sum_struct(in,div(summing(temp1,mult(temp2,2),mult(temp3,2),temp4),6));
%     output2 = in +((temp1 + temp2*2 + temp3*2 + temp4)/6);
end

