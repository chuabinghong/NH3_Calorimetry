%%
% created by Bing Hong CHUA 29Sep22Density

% script objective:
% extract pure phase specific heat from NH3 EOS based on measured wt%

%% SCRIPT ARGUMENTS

wt_str = [5.2, 8.2, 8.4, 10.0, 14.3, 20.07, 26.912]; % based off liquidus alignment
P_MPa = .101325;
T_K = 240:1:330;

%% GET CP

m = (1000.*(wt_str./100)./17.031)./(1-(wt_str./100));
in = {P_MPa,T_K,m};
out = SeaFreeze(in,'NH3');
cp = squeeze(out.Cp/1000);



t = array2table([T_K',cp]);

wt_str = convertStringsToChars(string(wt_str));
t.Properties.VariableNames(1:7) = {'T(K)',wt_str{1},wt_str{2},wt_str{3},wt_str{4},wt_str{5},wt_str{6}};
writetable(t,'SF_cp_liq.csv');



%%
%% GET CP
P_MPa = .101325;
T_K = 0:1:273;
in = {P_MPa,T_K};
out = SeaFreeze(in,'Ih');
cp = squeeze(out.Cp/1000);

figure()
hold on
plot(T_K,cp,LineWidth=5)
plot(iceCp{:,1},iceCp{:,2}/1000,'red--',LineWidth=5)

legend('SeaFreeze','FW06')
xlabel('T (K)')
ylabel('Cp (J/gK)')
%%
t = array2table([T_K',cp]);

wt_str = convertStringsToChars(string(wt_str));
t.Properties.VariableNames(1:7) = {'T(K)',wt_str{1},wt_str{2},wt_str{3},wt_str{4},wt_str{5},wt_str{6}};
writetable(t,'SF_cp_liq.csv');
