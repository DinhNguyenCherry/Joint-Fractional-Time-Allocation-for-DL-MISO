function [ isBreak, global_OptValue, global_OptValueChain, global_DownlinkRate_PerGroupPerUser1, global_DownlinkRate_PerGroupPerUser2 ] = ProposedAlgNOMAperzone( Pbs, Hk, Gm, Rate_Threshold, MaxIteration, Conven, Fixed_timegroup_assignment )
%PROPOSEDALG Summary of this function goes here
%   Detailed explanation goes here


N = size(Hk,1);
K = size(Hk,2)
M = size(Gm,2)

global_DownlinkRate_PerGroupPerUser1=0;
global_DownlinkRate_PerGroupPerUser2=0;

strdisp = ['For Pbs = ' num2str(10*log10(Pbs)) ' dB & RateThreshold = ' num2str(Rate_Threshold) ];

if (Fixed_timegroup_assignment)
    strdisp = ['Fixed time -- ' strdisp];
else
    strdisp = ['Opt. time -- ' strdisp];
end


        disp(['***** Getting a feasible point ..... ' strdisp]);

    
        [ GetBreak, OptimalValue, OptimalValue_preStep, DownlinkRate_PerGroupPerUser1, DownlinkRate_PerGroupPerUser2, W1_current, W2_current, alpha_current, beta_current, vartheta1_current, vartheta2_current ] = GetInitializationNOMAperzone( Pbs, Hk, Gm, Rate_Threshold, MaxIteration, Conven, Fixed_timegroup_assignment);
        


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
                Get_optSolutionPerIterationNOMAperzone(Pbs, Hk, Gm, Rate_Threshold, Conven, Fixed_timegroup_assignment, W1_current, W2_current, alpha_current, beta_current, vartheta1_current, vartheta2_current, 0);

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

%                 if (Conven)
%                     global_Sum_Downlink_power = real(trace(W1_current*W1_current'))
%                 else
%                     global_Sum_Downlink_power = real(trace([W1_current W2_current]*diag(global_time)*[W1_current W2_current]'))
%                 end
                

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

