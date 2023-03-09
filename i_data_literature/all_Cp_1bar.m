 %NH3=0.01703026;  %kg/molH2O=0.01801528;  %kg/molnw= 55.508435;  %moles of water in 1 kg of waterAllred1981_P_T_m_X_Cp=[  %  MPa K m(mol/kg) molfracX Cpo (J/K/mol)  % Cp_Allred=(m.*Cp0+Cp_water)./(1+m*M2);  M2=17.87371609e-3; %kg/mole  0.1  283.15  0.04755  0.00085589  NaN  0.1  283.15  0.07801  0.0014034  68.65  0.1  283.15  0.10587  0.0019036  68.28  0.1  283.15  0.12893  0.0023173  69.04  0.1  283.15  0.13411  0.0024102  69.99  0.1  283.15  0.17084  0.0030683  69.41  0.1  283.15  0.17115  0.0030738  69.72  0.1  283.15  0.19316  0.0034678  70.19  0.1  283.15  0.19964  0.0035837  69.99  0.1  283.15  0.20821  0.0037369  70.26  0.1  283.15  0.23808  0.0042708  69.7  0.1  283.15  0.25835  0.0046327  70.05  0.1  283.15  0.31079  0.0055678  69.62  0.1  283.15  0.3437  0.0061538  69.47  0.1  298.15  0.04356  0.00078413  NaN  0.1  298.15  0.06649  0.0011964  72.99  0.1  298.15  0.11  0.0019778  73.42  0.1  298.15  0.15601  0.0028027  74.23  0.1  298.15  0.18195  0.0032672  74.15  0.1  298.15  0.21196  0.003804  74.57  0.1  298.15  0.23261  0.0041731  75.03  0.1  298.15  0.28192  0.0050532  75.37  0.1  298.15  0.31392  0.0056236  75.67  0.1  298.15  0.3437  0.0061538  75.18  0.1  313.15  0.03472  0.0006251  NaN  0.1  313.15  0.05971  0.0010745  77.54  0.1  313.15  0.06407  0.0011529  78.83  0.1  313.15  0.08261  0.001486  78.14  0.1  313.15  0.1159  0.0020836  79.23  0.1  313.15  0.12255  0.0022029  79.43  0.1  313.15  0.15636  0.002809  78.98  0.1  313.15  0.17819  0.0031999  79.02  0.1  313.15  0.19692  0.003535  79.06  0.1  313.15  0.20519  0.0036829  NaN  0.1  313.15  0.22806  0.0040918  79.86  0.1  313.15  0.25635  0.004597  79.47  0.1  313.15  0.27855  0.0049931  79.15  0.1  313.15  0.28241  0.0050619  80.16  0.1  313.15  0.31392  0.0056236  79.83];Allred1981_P_T_m_X_Cp=Allred1981_P_T_m_X_Cp(find(~isnan(Allred1981_P_T_m_X_Cp(:,5))),:);outw=SeaFreeze(Allred1981_P_T_m_X_Cp(:,[1 2]),'water1');Cpw=outw.Cp;kgovmol=Allred1981_P_T_m_X_Cp(:,4)*NH3+(1-Allred1981_P_T_m_X_Cp(:,4))*H2O;Cp_allred=(Allred1981_P_T_m_X_Cp(:,3).*Allred1981_P_T_m_X_Cp(:,5)+Cpw)./(1+Allred1981_P_T_m_X_Cp(:,3)*NH3);Allred1981_P_T_m_X_Cp=[Allred1981_P_T_m_X_Cp(:,[1 2 4]) Cp_allred];%%Hildenbrand1953_P_T_X_Cp=[  %  MPa K molfracX (J/K/mole)   0.1  197.12  0.5  59.329  0.1  201.38  0.5  60.124  0.1  206.38  0.5  60.982  0.1  211.61  0.5  61.965  0.1  216.93  0.5  62.969  0.1  222.37  0.5  64.015  0.1  227.89  0.5  65.04  0.1  233.47  0.5  66.107  0.1  239.13  0.5  67.195  0.1  244.84  0.5  68.262  0.1  250.62  0.5  69.392  0.1  256.44  0.5  70.396  0.1  262.35  0.5  71.525  0.1  268.28  0.5  72.571  0.1  274.12  0.5  73.534  0.1  279.81  0.5  74.517  0.1  285.33  0.5  75.438  0.1  290.21  0.5  76.128  0.1  196.82  0.66667  65.842  0.1  201.08  0.66667  66.47  0.1  205.86  0.66667  67.251  0.1  211.08  0.66667  67.92  0.1  217.74  0.66667  68.938  0.1  223.36  0.66667  69.691  0.1  228.55  0.66667  70.514  0.1  234.14  0.66667  71.295  0.1  239.9  0.66667  72.048  0.1  245.65  0.66667  72.816  0.1  251.63  0.66667  73.583  0.1  257.78  0.66667  74.35  0.1  264  0.66667  75.075  0.1  270.16  0.66667  75.758  0.1  201.56  0.65393  66.14  0.1  207.3  0.65393  67.085  0.1  213.4  0.65393  67.957  0.1  219.35  0.65393  68.829  0.1  225.14  0.65393  69.629  0.1  237.26  0.65393  71.301  0.1  243.37  0.65393  72.1  0.1  249.63  0.65393  73.045  0.1  256.22  0.65393  73.844  0.1  201.3  0.62556  65.227  0.1  207.13  0.62556  66.101  0.1  213.35  0.62556  67.047  0.1  219.76  0.62556  67.994  0.1  225.87  0.62556  68.94  0.1  232.04  0.62556  69.887  0.1  238.18  0.62556  70.76  0.1  244.18  0.62556  71.561  0.1  250.31  0.62556  72.507  0.1  256.64  0.62556  73.381  0.1  201.33  0.60184  64.367  0.1  207.3  0.60184  65.388  0.1  213.72  0.60184  66.263  0.1  220.14  0.60184  67.356  0.1  226.42  0.60184  68.304  0.1  232.63  0.60184  69.324  0.1  238.82  0.60184  70.272  0.1  245.09  0.60184  71.22  0.1  251.46  0.60184  72.24  0.1  257.91  0.60184  73.261];kgovmol=Hildenbrand1953_P_T_X_Cp(:,3)*NH3+(1-Hildenbrand1953_P_T_X_Cp(:,3))*H2O;Hildenbrand1953_P_T_X_Cp(:,4)=Hildenbrand1953_P_T_X_Cp(:,4)./kgovmol;%%Chan1964_P_T_X_Cp=[  %  MPa K molfracX (J/K/mole)   0.1  183.72  0.3343  53.039  0.1  192.28  0.3343  54.308  0.1  200.07  0.3343  56.08  0.1  208.11  0.3343  57.809  0.1  215.9  0.3343  59.496  0.1  223.49  0.3343  61.184  0.1  230.93  0.3343  62.899  0.1  238.48  0.3343  64.489  0.1  245.97  0.3343  66.121  0.1  253.07  0.3343  67.502  0.1  259.93  0.3343  68.673  0.1  266.98  0.3343  70.305  0.1  274.34  0.3343  71.407  0.1  281.63  0.3343  72.662  0.1  288.15  0.3343  73.666];kgovmol=Chan1964_P_T_X_Cp(:,3)*NH3+(1-Chan1964_P_T_X_Cp(:,3))*H2O;Chan1964_P_T_X_Cp(:,4)=Chan1964_P_T_X_Cp(:,4)./kgovmol;Wrewsky1924_P_T_X_Cp=[  %  MPa K molfracX (J/K/mole)     0.1  275.55  0.014692  75.82  0.1  275.55  0.030415  75.537  0.1  275.55  0.039899  75.415  0.1  275.55  0.16562  72.991  0.1  275.55  0.24732  71.933  0.1  275.55  0.31316  71.559  0.1  275.55  0.41457  73.292  0.1  293.75  0.015537  74.408  0.1  293.75  0.030309  75.334  0.1  293.75  0.042426  75.299  0.1  293.75  0.08979  75.044  0.1  293.75  0.15804  74.321  0.1  293.75  0.25092  74.253  0.1  293.75  0.31306  74.708  0.1  293.75  0.33541  74.941  0.1  314.15  0.015537  75.259  0.1  314.15  0.042005  75.835  0.1  314.15  0.086123  75.839  0.1  314.15  0.15502  76.337  0.1  314.15  0.21917  76.514  0.1  334.05  0.030309  75.733  0.1  334.05  0.086332  76.341  0.1  334.05  0.12878  76.859];kgovmol=Wrewsky1924_P_T_X_Cp(:,3)*NH3+(1-Wrewsky1924_P_T_X_Cp(:,3))*H2O;Wrewsky1924_P_T_X_Cp(:,4)=Wrewsky1924_P_T_X_Cp(:,4)./kgovmol;  Chernenkaya1971_P_T_X_Cp=[  %  MPa K molfracX (J/K/mole)   0.1  298.15  0.021132  75.214  0.1  298.15  0.042215  75.277  0.1  298.15  0.06325  75.326  0.1  298.15  0.084236  75.329  0.1  298.15  0.10517  75.28  0.1  298.15  0.12606  75.156  0.1  298.15  0.14691  75.002  0.1  298.15  0.1677  74.857  0.1  298.15  0.18845  74.756  0.1  298.15  0.20915  74.656];kgovmol=Chernenkaya1971_P_T_X_Cp(:,3)*NH3+(1-Chernenkaya1971_P_T_X_Cp(:,3))*H2O;Chernenkaya1971_P_T_X_Cp(:,4)=Chernenkaya1971_P_T_X_Cp(:,4)./kgovmol;Fujita2008_P_T_X_Cp=[  %  MPa K molfracX (J/K/mole)   0.1  280  0.2939  68.954  0.1  280  0.2939  69.485  0.1  300  0.1511  73.253];kgovmol=Fujita2008_P_T_X_Cp(:,3)*NH3+(1-Fujita2008_P_T_X_Cp(:,3))*H2O;Fujita2008_P_T_X_Cp(:,4)=Fujita2008_P_T_X_Cp(:,4)./kgovmol;id=find(Fujita2008_P_T_X_Cp(:,1)<5);Fujita2008_P_T_X_Cp=Fujita2008_P_T_X_Cp(id,:);%%t = array2table(Wrewsky1924_P_T_X_Cp);t.Properties.VariableNames(1:4) = {'P(MPa)','T(K)','molefrac','cp(J/kgK)'};sd = 'C:\1_data\OneDrive - Nanyang Technological University\OFYP\CalorimetryAnalysis\i_data_literature\Baptiste\';writetable(t,append(sd,'wrewsky_1924','.csv'));