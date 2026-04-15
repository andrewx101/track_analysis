%% Load data first



%% Plot chi contour along with msd
dt=msd(1,1);
lags=chi4_d{1}.*dt;
d=chi4_d{2};
chi4=chi4_d{3};

close all
h=figure;
ax=axes(h);
h.Units="inches";
h.Position(3)=3.37;
h.Position(4)=2.53;
hold(ax,'on');

[X,Y]=meshgrid(lags(:),d(:));
Z=chi4;
% Z(Z<0)=0;
[~,hc]=contourf(ax,X,Y,Z);
hc.LevelList=linspace(min(Z(:)),1.5e-3,250);

ax.XScale='log';
ax.YScale='log';


l=plot(msd(:,1),sqrt(msd(:,2)),'r--','LineWidth',1);

hc.LineStyle='none';
hc.Fill='on';

cb=colorbar(ax);

cb.LineWidth=1;
cb.FontSize=10;
clim(ax,[0,3e-3]);

hold(ax,'off');

ax.Box='on';
ax.LineWidth=1;
ax.FontSize=10;
ax.XScale='log';
ax.YScale='log';
ax.XLabel.String='lag times (s)';
ax.YLabel.String='tolarance (m)';
ax.XLabel.FontSize=12;
ax.YLabel.FontSize=12;

ax.XMinorTick='on';
ax.YMinorTick='on';

cb.Position(2)=0.2;
ax.Position(1)=0.2;
ax.Position(2)=0.2;
ax.Position(3)=(cb.Position(1)-ax.Position(1)).*0.95;
ax.Position(4)=(1-ax.Position(2))*0.95;
cb.Position(4)=ax.Position(4)*0.9;


h.Color='none';
exportgraphics(h,'PICs-chi4-2a05-0.80-2024.emf',"ContentType",'vector');

% Extract chi values on msd
Z_r=interp2(X,Y,Z,msd(:,1),sqrt(msd(:,2)));

%% Plot gel point
% close all
% symbol_list='sod^<>h';
% h=figure;
% ax=axes(h);
% h.Units="inches";
% h.Position(3)=3.37;
% h.Position(4)=2.53;
% 
% hold(ax,'on');
% l=cell(length(freq),1);
% for i=1:length(freq)
%     l{i}=plot(g_stages,g_tan_trans(:,i),['-' symbol_list(i)],...
%         'LineWidth',1,'DisplayName',num2str(freq(i)));
% end
% c=errorbar(cross(1),cross(2),cross(4),cross(4),cross(3),cross(3),'k.',...
%     'LineWidth',1);
% hold(ax,'off');
% 
% ax.YScale='log';
% 
% ax.FontSize=10;
% ax.Box='on';
% ax.LineWidth=1;
% 
% ax.XLabel.String="\itc\rm (mg/ml)";
% ax.YLabel.String="tan\itδ\rm";
% 
% lg=legend(ax);
% lg.Location="southwest";
% lg.Box='off';
% lg.FontSize=8;
% 
% msg={join([num2str(n(1),2) "±" num2str(n(2),1)]),...
%     join([num2str(cross(1),4) "±" num2str(cross(3),1) "mg/ml"])};
% 
% t=annotation(h,'textbox');
% t.FitBoxToText='on';
% t.String=msg;
% t.LineStyle='none';
% t.FontSize=8;
% t.Position(1)=ax.Position(1)+ax.Position(3)-t.Position(3);
% t.Position(2)=ax.Position(2)+ax.Position(4)-t.Position(4);



