function [ OptimalValue, DownlinkRate_PerUser1, DownlinkRate_PerUser2, W1, W2, alpha, beta, vartheta1, vartheta2 ] = Get_optSolutionPerIterationNOMA1( Pbs, Hk, Gm, Rate_Threshold, Conven, Fixed_timegroup_assignment, W1_current, W2_current, alpha_current, beta_current, vartheta1_current, vartheta2_current, Init )
%GET_OPTSOLUTIONPERITERATION Summary of this function goes here
%   Detailed explanation goes here

if (nargin<18)
    eta = 1;
end

max_min = 0;

N = size(Hk,1);
K = size(Hk,2);
M = size(Gm,2);


sigma_K = 0.01*ones(1, max(K,M));

W1_next = sdpvar(N*K,1,'full','complex');
W2_next = sdpvar(N*M,1,'full','complex');

alpha_next = ones(1,1);
beta_next = ones(1,1);

Theta1_next = sdpvar(K, 1, 'full', 'real');
Theta2_next = sdpvar(M, 1, 'full', 'real');
Theta21_next = sdpvar(K, 1, 'full', 'real');

vartheta1_next = sdpvar(K, 1, 'full', 'real');
vartheta2_next = sdpvar(M, 1, 'full', 'real');

varrho = sdpvar(1, 1);
phi = sdpvar(1, 1);

obj = 0;
cons = [];

%% add objective function and constraints for DOWNLINK rate in inner zone

PerGroupDownlinkRate1 = sdpvar(K,1);

tempsum1_down = sdpvar(K,1);


for k = 1:1:K
    
    temp_down = 0;
    
          
		
		% get parameters 

        [nu, xi, lambda] = GetParameter_Downlink(vartheta1_current(k), alpha_current);
		
		% add first term 
		
        PerGroupDownlinkRate1(k) = real(nu);

        temp_down = temp_down + real(nu);
		
		
		% add second term
		

        PerGroupDownlinkRate1(k) = PerGroupDownlinkRate1(k) + real(xi)*vartheta1_next(k);

        temp_down = temp_down + real(xi)*vartheta1_next(k);
        
        
        % add third term 
        
        PerGroupDownlinkRate1(k) = PerGroupDownlinkRate1(k) + real(lambda)*alpha_next;

        temp_down = temp_down + real(lambda)*alpha_next;
        
        if (~Init)
        
        cons = [cons, PerGroupDownlinkRate1(k)>=0];
        
        end


    if (Init || max_min==1)
        
        cons = [cons, real(temp_down/Rate_Threshold)>= varrho];
        
    else
        
        if (max_min==0)

            cons = [cons, real(temp_down/log(2))>=Rate_Threshold];
            cons = [cons, temp_down>=tempsum1_down(k)];
            
        else
            
%             cons = [cons, real(temp_down/Rate_Threshold)>= phi_down];
           % cons = [cons, real(temp_down/Rate_Threshold)>= eta*phi];
            
        end
        
    end
    
    
end


obj = obj + sum(tempsum1_down);



%% add objective function and constraints for DOWNLINK rate in outer zone

if (~Conven)

PerGroupDownlinkRate2 = sdpvar(M,1);

tempsum2_down = sdpvar(M,1);

for m = 1:1:M
    
    temp_down = 0;
    
		
		% get parameters 

        [nu, xi, lambda] = GetParameter_Downlink(vartheta2_current(m), beta_current);
		
		% add first term 
		
        PerGroupDownlinkRate2(m) = real(nu);

        temp_down = temp_down + real(nu);
		% add second term 
		

        PerGroupDownlinkRate2(m) = PerGroupDownlinkRate2(m) + real(xi)*vartheta2_next(m);

        temp_down = temp_down + real(xi)*vartheta2_next(m);
        
        
        % add third term 
        
        PerGroupDownlinkRate2(m) = PerGroupDownlinkRate2(m) + real(lambda)*beta_next;

        temp_down = temp_down + real(lambda)*beta_next;
        
        if (~Init)
        
        cons = [cons, PerGroupDownlinkRate2(m)>=0];
        
        end



    if (Init || max_min==1)
        
        cons = [cons, real(temp_down/Rate_Threshold)>= varrho];
        
    else
        
        if (max_min==0)

            cons = [cons, real(temp_down/log(2))>=Rate_Threshold];
            cons = [cons, temp_down>=tempsum2_down(m)];
            
        else
            
%             cons = [cons, real(temp_down/Rate_Threshold)>= phi_down];
            cons = [cons, real(temp_down/Rate_Threshold)>= eta*phi];
            
        end
        
    end
    
    
end


obj = obj + sum(tempsum2_down);

end



%% Add power constraints

   
    W1_next_mat_all = reshape(W1_next, N, K);
    W2_next_mat_all = reshape(W2_next, N, M);
    
    for k = 1:1:K
    
        W1_next_mat = W1_next_mat_all(:,[1:(k-1) (k+1):K]);  
        
        % Constraints
        
        cons = [cons, real(Hk(:,k)'*W1_next(((k-1)*N+1):(k*N)))>=0 ];
        
        cons = [cons, 2*real(Hk(:,k)'*W1_current(((k-1)*N+1):(k*N)))*real(Hk(:,k)'*W1_next(((k-1)*N+1):(k*N)))-(real(Hk(:,k)'*W1_current(((k-1)*N+1):(k*N))))^2>=Theta1_next(k) ];
        
        cons = [cons, Theta1_next(k)>=0];
        
        
        % Constraints 
        
%         cons = [cons, vartheta_next(k,g)>0];
        
        Sum_vec = [(Hk(:,k)'*W1_next_mat) sigma_K(k)];
        Sum_vec = [Sum_vec  (Hk(:,k)'*W2_next_mat_all(:,[1:(k-1) (k+1):K]))];
		
        cons = [cons, cone([Sum_vec 0.5*(Theta1_next(k)-vartheta1_next(k))], 0.5*(Theta1_next(k)+vartheta1_next(k)) )];
        
        
    end


    
    
    for m = 1:1:M
    
        W2_next_mat = W2_next_mat_all(:,[1:(m-1) (m+1):M]);  
        
        % Constraints 
        
        cons = [cons, real(Gm(:,m)'*W2_next(((m-1)*N+1):(m*N)))>=0 ];
        %cons = [cons, real(Hk(:,m)'*W2_next(((m-1)*N+1):(m*N)))>=0 ];
        
        cons = [cons, 2*real(Gm(:,m)'*W2_current(((m-1)*N+1):(m*N)))*real(Gm(:,m)'*W2_next(((m-1)*N+1):(m*N)))-(real(Gm(:,m)'*W2_current(((m-1)*N+1):(m*N))))^2>=Theta2_next(m) ];
        cons = [cons, 2*real(Hk(:,m)'*W2_current(((m-1)*N+1):(m*N))*conj(Hk(:,m)'*W2_next(((m-1)*N+1):(m*N))))-(norm(Hk(:,m)'*W2_current(((m-1)*N+1):(m*N))))^2>=Theta21_next(m) ];
        
        
        cons = [cons, Theta2_next(m)>=0];
        cons = [cons, Theta21_next(m)>=0];
        
        % Constraints 
        
        
        Sum_vec2 = [(Gm(:,m)'*W2_next_mat) sigma_K(m)];
		Sum_vec2 = [Sum_vec2  (Gm(:,m)'*W1_next_mat_all(:,[1:m  m+1:K]))];
        cons = [cons, cone([Sum_vec2 0.5*(Theta2_next(m)-vartheta2_next(m))], 0.5*(Theta2_next(m)+vartheta2_next(m)) )];
        
        Sum_vec21 = [(Hk(:,m)'*W2_next_mat) sigma_K(m)]; 
		Sum_vec21 = [Sum_vec21  (Hk(:,m)'*W1_next_mat_all(:,[1:m  m+1:K]))];
        cons = [cons, cone([Sum_vec21 0.5*(Theta21_next(m)-vartheta2_next(m))], 0.5*(Theta21_next(m)+vartheta2_next(m)) )];
        
    end





% Constraints - power constraint with Pbs

if (~Conven)

%alpha_tilde_next = sdpvar(1, 1, 'full', 'real');
%beta_tilde_next = sdpvar(1, 1, 'full', 'real');

cons = [cons, cone([W1_next(:); W2_next(:)], sqrt(Pbs ))];

%cons = [cons, cone([W2_next(:); 0.5*(beta_tilde_next-beta_next)], 0.5*(beta_tilde_next+beta_next)) ];

%cons = [cons, alpha_tilde_next + beta_tilde_next <= (2*real(W1_current(:)'*W1_next(:))/beta_current - norm(W1_current(:))^2/(beta_current^2)*beta_next + Pbs) ];

else
    
%     cons = [cons, cone(W1_next, sqrt(Pbs) ) ];
    
end




if (Init || max_min==1)
    
%     myops = sdpsettings('solver','mosek');
    myops = sdpsettings('solver','sedumi','verbose',0);
    diagnotics = solvesdp(cons, -varrho, myops);
    
    OptimalValue = double(varrho);
    
else
    
    if (max_min==0)

    %     myops = sdpsettings('solver','mosek');
        myops = sdpsettings('solver','sedumi','verbose',0);
        diagnotics = solvesdp(cons, -obj, myops);

        OptimalValue = double(obj);
%         sumrate_up = double(sum(tempsum_up));
%         sumrate_down = double(sum(tempsum_down));
    else
        
%         phi = phi_up + phi_down;
        
        %     myops = sdpsettings('solver','mosek');
        myops = sdpsettings('solver','sedumi','verbose',0);
        diagnotics = solvesdp(cons, -phi, myops);

        OptimalValue = double(phi);
        
    end
    
end




DownlinkRate_PerUser1 = real(double(PerGroupDownlinkRate1));

if(~Conven)
    
DownlinkRate_PerUser2 = real(double(PerGroupDownlinkRate2));

else
    DownlinkRate_PerUser2 = 0;  
end




W1 = double(W1_next);
W2 = double(W2_next);



alpha = double(alpha_next);
beta = double(beta_next);
vartheta1 = double(vartheta1_next)
vartheta2 = double(vartheta2_next)

norm(W1)^2/alpha + norm(W2)^2/beta;

end

