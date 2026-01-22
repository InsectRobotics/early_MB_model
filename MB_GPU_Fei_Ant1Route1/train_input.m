function [I_PN, BA] = train_input( t, input, reward, null_input, C_I_PN_sen, C_I_PN_var)
%GL_INPUT According to the inputs and time, set GL and BA
% t is the current time; total training time is 1000ms

BA = 0;

if t < 15.01           % [ms]
   I_PN = input*C_I_PN_var + C_I_PN_sen;
   if t == 15.0 && reward == 1
        BA = 0.5;
   end
elseif t <= 30
    I_PN = null_input;
end
end

