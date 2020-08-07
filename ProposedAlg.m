%%% Written by Van-Dinh Nguyen
%%% Authors: Van-Dinh Nguyen, Hoang Duong Tuan, Trung Q. Duong, Oh-Soon Shin, and H. Vincent Poor
%%% Corresponding to the paper: "Joint Fractional Time Allocation and Beamforming for Downlink Multiuser MISO Systems," 
% % % IEEE Communications Letters, vol. 21, no. 12, pp. 2650-2653, Dec. 2017
%%%  

function [ isBreak, global_OptValue1,global_OptValue2, global_OptValueChain, global_DownlinkRate_PerGroupPerUser1, global_DownlinkRate_PerGroupPerUser2 ] = ProposedAlg( Pbs, Hk, Gm, Rate_Threshold, MaxIteration, Conven, Fixed_timegroup_assignment )
% function [ isBreak, global_OptValue, global_OptValueChain, global_DownlinkRate_PerGroupPerUser1, global_DownlinkRate_PerGroupPerUser2 ] = ProposedAlg( Pbs, Hk, Gm, Rate_Threshold, MaxIteration, Conven, Fixed_timegroup_assignment )
%PROPOSEDALG Summary of this function goes here
%   Detailed explanation goes here


N = size(Hk,1);
K = size(Hk,2)
M = size(Gm,2)
sigma_K = 0.01*ones(1, max(K,M));

global_DownlinkRate_PerGroupPerUser1=0;
global_DownlinkRate_PerGroupPerUser2=0;

strdisp = ['For Pbs = ' num2str(10*log10(Pbs)) ' dB & RateThreshold = ' num2str(Rate_Threshold) ];

if (Fixed_timegroup_assignment)
    strdisp = ['Fixed time -- ' strdisp];
else
    strdisp = ['Opt. time -- ' strdisp];
end


        disp(['***** Getting a feasible point ..... ' strdisp]);

    
        [ GetBreak, OptimalValue, OptimalValue_preStep, DownlinkRate_PerGroupPerUser1, DownlinkRate_PerGroupPerUser2, W1_current, W2_current, alpha_current, beta_current, vartheta1_current, vartheta2_current ] = GetInitialization( Pbs, Hk, Gm, Rate_Threshold, MaxIteration, Conven, Fixed_timegroup_assignment);
        


        %% Loop for optimization

        if (GetBreak)

            isBreak = 1;
            global_OptValue = 0;
            global_OptValueChain = [];

            global_UplinkRate_PerGroupPerUser = 0;
            global_DownlinkRate_PerGroupPerUser = 0;
            
            

        else

            isBreak = 0;

            disp(['---------------------------------- ' strdisp ' - Run iterative algorithm --------------------------------']);

            n = 0;
            OptimalValue_current = OptimalValue;
            if (abs(OptimalValue_preStep)<10)
                OptValueChain = [OptimalValue_preStep OptimalValue];

            else
                OptValueChain = [OptimalValue];

            end

            while (n<MaxIteration)

                disp(['******************* ' strdisp ' --- Iteration: ' num2str(n+1) ' *********************']);

                [ OptimalValue, DownlinkRate_PerGroupPerUser_next1, DownlinkRate_PerGroupPerUser_next2, W1_next, W2_next, alpha_next, beta_next, vartheta1_next, vartheta2_next] = ... 
                Get_optSolutionPerIteration(Pbs, Hk, Gm, Rate_Threshold, Conven, Fixed_timegroup_assignment, W1_current, W2_current, alpha_current, beta_current, vartheta1_current, vartheta2_current, 0);

                % Update

                if (~Conven)
                time_next = [1/alpha_next 1/beta_next]; 
                else
                    time_next = 1;
                end
                
                if (isFeasible(time_next, OptimalValue_current, OptimalValue, sum(DownlinkRate_PerGroupPerUser_next1/log(2),2), sum(DownlinkRate_PerGroupPerUser_next2/log(2),2), Rate_Threshold, Conven))


                    W1_current = W1_next;
                    W2_current = W2_next;

                    
                    alpha_current = alpha_next;
                    beta_current = beta_next;
                    
                    vartheta1_current = vartheta1_next;
                    vartheta2_current = vartheta2_next;
                    
                    
                    DownlinkRate_PerGroupPerUser1 = DownlinkRate_PerGroupPerUser_next1
                    DownlinkRate_PerGroupPerUser2 = DownlinkRate_PerGroupPerUser_next2


                    % Check convergence

 
                    if (checkConvergence(OptValueChain, OptimalValue))
                        break;
                    else
                        OptimalValue_current = OptimalValue
                        OptValueChain = [OptValueChain OptimalValue_current];

                    end

                else

                    disp('Infeasible point --> keeping the latest feasible point');
                    OptValueChain = [OptValueChain OptimalValue_current];

                    break;

                end

                n = n + 1;

            end
            
            disp('-------------------- Converging point ----------------------------');

            OptimalValue_permode = 1/log(2)*OptimalValue_current 
            
            DownlinkRate_PerGroupPerUser1 = DownlinkRate_PerGroupPerUser1/log(2)
            DownlinkRate_PerGroupPerUser2 = DownlinkRate_PerGroupPerUser2/log(2)
%             DownlinkRate_PerUser = sum(DownlinkRate_PerGroupPerUser,2)
                                
            OptValueChain = OptValueChain/log(2)
           

            
            disp('------------------------------------------------------------------');
            

                global_time = [1/alpha_current 1/beta_current];
                global_OptValue = OptimalValue_permode;
                
                               
                global_DownlinkRate_PerGroupPerUser1 = DownlinkRate_PerGroupPerUser1;
                global_DownlinkRate_PerGroupPerUser2 = DownlinkRate_PerGroupPerUser2;

% %%% Best rate decoded by UEs in other zones
                 W1_next_mat_all = reshape(W1_next, N, M);
                 Rzone1_leaked = [];  
                 for k = 1:1:K
                    W_next_mat = W1_next_mat_all(:,[1:(k-1) (k+1):K]);
                    for m = 1:1:M
                      Sum_vec = [(Gm(:,m)'*W_next_mat) sigma_K(m)];


                      Rzone1_leaked =[Rzone1_leaked (1/alpha_current)*log2(1+norm(Gm(:,m)'*W1_next(((k-1)*N+1):(k*N)))^2/norm(Sum_vec)^2)];          
                
                    end

                 end
                global_OptValue1 = max(Rzone1_leaked)
                
                  W2_next_mat_all = reshape(W2_next, N, M);
                 Rzone2_leaked = [];  
                 for m = 1:1:M
                    W_next_mat = W2_next_mat_all(:,[1:(m-1) (m+1):M]);
                    for k = 1:1:K
                      Sum_vec = [(Hk(:,k)'*W_next_mat) sigma_K(k)];


                      Rzone2_leaked =[Rzone2_leaked (1/beta_current)*log2(1+norm(Hk(:,k)'*W2_next(((m-1)*N+1):(m*N)))^2/norm(Sum_vec)^2)];          

                    end

                 end
                global_OptValue2 = max(Rzone2_leaked)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
   
   global_OptValueChain = OptValueChain;

        end



        if (global_OptValue>0)

            disp([' ################### FINAL RESULT -- ' strdisp ' #############################']);
            
            
            global_time 
            global_OptValue 
            

            global_DownlinkRate_PerGroupPerUser1
            global_DownlinkRate_PerGroupPerUser2

            
            global_OptValueChain
            disp(' ##############################################################');
        end

end

function [check] = isFeasible(time_next, OptimalValue_current, OptimalValue, DownlinkRate1, DownlinkRate2, Threshold, Conven)

% check if the current point is feasible and the next point is infeasible

check = 1;

if ((sum(( (DownlinkRate1-Threshold) >= -10^(-1) )) < length(DownlinkRate1)) || (sum(( (DownlinkRate1) < 100 )) < length(DownlinkRate1)))
    check = 0;
    disp('Downlink1 --> Not pass Threshold check');
end

if (~Conven)
if ((sum(( (DownlinkRate2-Threshold) >= -10^(-1) )) < length(DownlinkRate2)) || (sum(( (DownlinkRate2) < 100 )) < length(DownlinkRate2)))
    check = 0;
    disp('Downlink2 --> Not pass Threshold check');
end
end

check_time =1;

if (sum(time_next>1)>0 || sum(time_next<0)>0)
    check_time = 0;
    disp('Time_next --> Not pass range check');
end


check_time_next = 1;



check3 = 1;

if ((OptimalValue_current>0 && OptimalValue<0) )%|| (round(OptimalValue_current,20)==0) || (OptimalValue_current>10000))
    check3 = 0;
    disp('OptimalValue --> Not pass range check');
    OptimalValue_current
    OptimalValue
end

check = (check_time && check_time_next && check3 && check);

end


function [check] = checkConvergence(OptValueChain, OptimalValue)

check = 0;

if (length(OptValueChain)>10)
    if (abs(OptimalValue-OptValueChain(length(OptValueChain))) < 10^-5)
        check = 1;
    end
%     if (abs(OptimalValue-OptValueChain(length(OptValueChain)-5)) < 0.1)
%         check = 1;
%     end
end

end

