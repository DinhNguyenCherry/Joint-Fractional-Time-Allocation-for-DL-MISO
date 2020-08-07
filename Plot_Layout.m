function [ done ] = Plot_Layout( RadiusOfCell, InnerZoneRadius, positionDownlinkUsers, K )
%PLOT_LAYOUT Summary of this function goes here
%   Detailed explanation goes here

hold on

ang=0:0.01:2*pi; 
xp=RadiusOfCell*cos(ang);
yp=RadiusOfCell*sin(ang);
ax = plot(0+xp,0+yp);

xp=InnerZoneRadius*cos(ang);
yp=InnerZoneRadius*sin(ang);
bx = plot(0+xp,0+yp);

scatter(0,0 , 'k^ ')

if (~isempty(positionDownlinkUsers))
    scatter(positionDownlinkUsers(1:K,1), positionDownlinkUsers(1:K,2), 'r* ')
    scatter(positionDownlinkUsers(K+1:end,1), positionDownlinkUsers(K+1:end,2), 'bs ')
    for iUser = 1:1:size(positionDownlinkUsers,1)
        
        Order = iUser;
        
        if (Order>K)
            Order = Order - K;
        end
        
        text(positionDownlinkUsers(iUser,1)-30, positionDownlinkUsers(iUser,2)-10, num2str(Order));
        
    end
end


set(gca,'XTick',[-500 :100: 500])
set(gca,'YTick',[-500 :100: 500])

xlim([-RadiusOfCell-5 RadiusOfCell+5])
ylim([-RadiusOfCell-5 RadiusOfCell+5])
axis('square')

done = 1;


end

