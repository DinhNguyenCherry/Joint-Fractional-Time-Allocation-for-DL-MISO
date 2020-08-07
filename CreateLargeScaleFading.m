function [ D_Hk, D_Gm, positionDownlinkUsers ] = CreateLargeScaleFading( K, M, Parameters )
%CREATED Summary of this function goes here
%   Detailed explanation goes here

RadiusOfCell = Parameters(1);
InnerZoneRadius = Parameters(2);
RadiusOfNearestUser = Parameters(3);
StandardDeviation = 10^(Parameters(4)/10);
ploss = Parameters(5);

L = K + M;


%% large-fading for K downlink users - inner zone

    Zvector = StandardDeviation*randn(1,K);
    rvector = RadiusOfNearestUser*ones(1,K) + (InnerZoneRadius-RadiusOfNearestUser)*rand(1,K);
    anglevector = 2*pi*rand(1,K);
    
    positionDownlinkUsers = [(rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];
        
    distanceBS_To_DownlinkUsers = sqrt(sum(positionDownlinkUsers'.^2));
        
    betavector_downlink = (10.^(Zvector/10))./((distanceBS_To_DownlinkUsers/RadiusOfNearestUser).^ploss);

        
    D_Hk = (diag(betavector_downlink)).^(0.5)
    
    
%% large-fading for M downlink users  - outer zone

    Zvector = StandardDeviation*randn(1,M);
    rvector = InnerZoneRadius*ones(1,M) + (RadiusOfCell-InnerZoneRadius)*rand(1,M);
    anglevector = 2*pi*rand(1,M);
    
    positionDownlinkUsers = [positionDownlinkUsers; (rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];

        
    distanceBS_To_DownlinkUsers = sqrt(sum(positionDownlinkUsers(K+1:L)'.^2));
        
    betavector_downlink = (10.^(Zvector/10))./((distanceBS_To_DownlinkUsers/RadiusOfNearestUser).^ploss);

        
    D_Gm = (diag(betavector_downlink)).^(0.5)
    
    

end

