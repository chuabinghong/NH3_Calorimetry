size = 40;
div = 2.5;

f = figure();
f.Position = [100 100 1000 800];

T = readtable('../i_data_literature/lit_cp/chan_1964.csv');

% tbl = chan1964;
% plot3(tbl{:,2},tbl{:,4}*100,tbl{:,5}/1000,'.',markersize=size)
% hold on
% tbl = hildenbrand1953;
% plot3(tbl{:,2},tbl{:,4}*100,tbl{:,5}/1000,'.',markersize=size)
% tbl = wrewsky1924;
% plot3(tbl{:,2},tbl{:,4}*100,tbl{:,5}/1000,'.',markersize=size)
% tbl = chernenkaya1971;
% plot3(tbl{:,2},tbl{:,4}*100,tbl{:,5}/1000,'.',markersize=size)
% tbl = allred1981;
% plot3(tbl{:,2},tbl{:,4}*100,tbl{:,5}/1000,'.',markersize=size)
% tbl = fujita2008;
% plot3(tbl{:,2},tbl{:,4}*100,tbl{:,5}/1000,'.',markersize=size)


% tbl = zliquiduscp;
% plot3(tbl{:,1},tbl{:,2}*100,tbl{:,3},Color='k',linewidth=5)

% BaptisteTRF = readtable('TF98_cp_all');
% T = 180:1:340;
% wt = 0:1:100;

TRFcut = BaptisteTRF{:,6:31};
T = 180:1:340;
wt = 5:1:30;

[X,Y] = meshgrid(T,wt);
C = Y;
h = surf(X,Y,(TRFcut(:,:)/1000)',C,'FaceAlpha',.6,EdgeColor = 'k');

hold on

Cdata = h.CData;
colormap(flipud(turbo))
cmap = colormap;

 % make it into a index image.
 cmin = min(Cdata(:));
 cmax = max(Cdata(:));
 m = length(cmap);

c = fix((5.2-cmin)/(cmax-cmin)*m)+1;
c1 = fix((8.2-cmin)/(cmax-cmin)*m)+1;
c2 = fix((8.4-cmin)/(cmax-cmin)*m)+1;
c3 = fix((10.0-cmin)/(cmax-cmin)*m)+1;
c4 = fix((14.3-cmin)/(cmax-cmin)*m)+1;
c5 = fix((20.07-cmin)/(cmax-cmin)*m)+1;
c6 = fix((26.912-cmin)/(cmax-cmin)*m)+1;

dd = '../i_data_processed/';
data = readtable(strcat(dd,'5.2wt%_cp_cut_pure_4.5386g.csv'));
data1 = readtable(strcat(dd,'8.2wt%_cp_cut_pure_4.1943g.csv'));
data2 = readtable(strcat(dd,'8.4wt%_cp_cut_pure_4.5858g.csv'));
data3 = readtable(strcat(dd,'10.0wt%_cp_cut_pure_4.5202g.csv'));
data4 = readtable(strcat(dd,'14.3wt%_cp_cut_pure_3.8153g.csv'));
data5 = readtable(strcat(dd,'20.07wt%_cp_cut_pure_3.7107g.csv'));
data6 = readtable(strcat(dd,'26.912wt%_cp_cut_pure_3.7778g.csv'));


tbl = data;
plot3(tbl{:,1},tbl{:,2}*100,tbl{:,4},'o',markersize=size/div,Color='k',MarkerFaceColor=cmap(c,:))
hold on
tbl = data1;
plot3(tbl{:,1},tbl{:,2}*100,tbl{:,4},'o',markersize=size/div,Color='k',MarkerFaceColor=cmap(c1,:))
tbl = data2;
plot3(tbl{:,1},tbl{:,2}*100,tbl{:,4},'o',markersize=size/div,Color='k',MarkerFaceColor=cmap(c2,:))
tbl = data3;
plot3(tbl{:,1},tbl{:,2}*100,tbl{:,4},'o',markersize=size/div,Color='k',MarkerFaceColor=cmap(c3,:))
tbl = data4;
plot3(tbl{:,1},tbl{:,2}*100,tbl{:,4},'o',markersize=size/div,Color='k',MarkerFaceColor=cmap(c4,:))
tbl = data5;
plot3(tbl{:,1},tbl{:,2}*100,tbl{:,4},'o',markersize=size/div,Color='k',MarkerFaceColor=cmap(c5,:))
tbl = data6;
plot3(tbl{:,1},tbl{:,2}*100,tbl{:,4},'o',markersize=size/div,Color='k',MarkerFaceColor=cmap(c6,:))


xlim([210 315])
zlim([3.4 4.6])
ylim([5 30])

xlabel('Temperature (K)')
ylabel('Mass Fraction (wt%)')
zlabel('Specific Heat  (J g^-^1 K^-^1)')
leg = legend('Equation of State','Data','Location','southeast');
set(findall(gcf,'-property','FontSize'),'FontSize',30)

c= colorbar;
c.Label.String = 'Mass Fraction (wt%)';
title('Liquid Phase Specific Heat')
view(0,0)

saveas(gcf,'../o_supplementaryPlots/3D.png')