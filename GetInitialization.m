function [ GetBreak, OptimalValue, OptimalValue_preStep, DownlinkRate_PerGroupPerUser1, DownlinkRate_PerGroupPerUser2, W1_current, W2_current, alpha_current, beta_current, vartheta1_current, vartheta2_current ] = GetInitialization( Pbs, Hk, Gm, RateThreshold, MaxIteration, Conven, Fixed_timegroup_assignment )
%GETINITIALIZATION Summary of this function goes here
%   Detailed explanation goes here


N = size(Hk,1);
K = size(Hk,2);
M = size(Gm,2);

OptimalValue_preStep = 0;
OptimalValue = 0;

n=0;
GetBreak = 0;


W1_current = randn(N*K,1)+1i*randn(N*K,1);
W2_current = randn(N*M,1)+1i*randn(N*M,1);

if (~Conven)
    alpha_current = 2;
    beta_current = 2;
    x = 1/alpha_current;
    y = 1/beta_current;
else
    alpha_current = 1;
    beta_current = 1;
    x = 1;
    y = 0;
end

vartheta1_current = rand(K, 1);
vartheta2_current = rand(M, 1);

while (~(OptimalValue>=log(2) || n >= MaxIteration) )%(~(OptimalValue>=log(2) || n >= MaxIteration) )
    
    if (n>0)
        OptimalValue_preStep = OptimalValue;
    end


    
    [ OptimalValue, DownlinkRate_PerGroupPerUser1, DownlinkRate_PerGroupPerUser2, W1_next, W2_next, alpha_next, beta_next, vartheta1_next, vartheta2_next] = ... 
        Get_optSolutionPerIteration(Pbs, Hk, Gm, RateThreshold, Conven, Fixed_timegroup_assignment, W1_current, W2_current, alpha_current, beta_current, vartheta1_current, vartheta2_current, 1);
    
%     OptimalValue
    
    W1_current = W1_next;
    W2_current = W2_next;
    
    alpha_current = alpha_next;
    beta_current = beta_next;
    
    if (~Conven)
        x = 1/alpha_current;
        y = 1/beta_current;
    else
        x = 1;
        y = 0;
    end
    
    
    vartheta1_current = vartheta1_next;
    vartheta2_current = vartheta2_next;
    
    
    disp(['---- For Pbs = ' num2str(10*log10(Pbs)) ' dB & RateThreshold = ' num2str(RateThreshold)]);
    n = n+1
    
%     if (n>MaxIteration || (round(OptimalValue)==0) || (OptimalValue>1000))
%         GetBreak = 1;
%         break;
%     end
    
end

disp(['Iterations: ' num2str(n) ]);
OptimalValue_preStep
OptimalValue 
DownlinkRate_PerGroupPerUser1
DownlinkRate_PerGroupPerUser2


end

