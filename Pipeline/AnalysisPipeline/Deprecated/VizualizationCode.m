%% Visualization Code


% %% Condition names
% 
% 
% 
% %% Condition names
% 
% CNames = {'Control 1', 'Control 2', 'Control 3', ...
%           'Day 0 - Well 1', 'Day 0 - Well 2','Day 0 - Well 3',...
%           'Day 2 - Well 1', 'Day 2 - Well 2','Day 2 - Well 3'};
% 
% for i = 1:9
%     Ch1R(i).Expression = Results_Ch1{i}.Expression;
%     Ch2R(i).Expression = Results_Ch2{i}.Expression;
%     Ch1R(i).F          = Results_Ch1{i}.SumF;
%     Ch2R(i).F          = Results_Ch2{i}.SumF;
%     Ch1R(i).Fcell      = Results_Ch1{i}.FPerCellA;
%     Ch2R(i).Fcell      = Results_Ch2{i}.FPerCellA;
%     Ch1R(i).Fvol       = Results_Ch1{i}.FPerCellV;
%     Ch2R(i).Fvol       = Results_Ch2{i}.FPerCellV;
%     Ch1R(i).Fvox       = Results_Ch1{i}.FPerVox;
%     Ch2R(i).Fvox       = Results_Ch2{i}.FPerVox;
% end
% 
% 
% %% All
% ALLRESULTS = [];
% 
% for i = 1:9
%     Ch1Class           = logical([Ch1R(i).Expression{:}]');
%     Ch2Class           = logical([Ch2R(i).Expression{:}]');
%     ALLRESULTS(i).CH1  = Ch1Class & ~Ch2Class;
%     ALLRESULTS(i).CH2  = Ch2Class & ~Ch1Class;
%     ALLRESULTS(i).Both = Ch1Class & Ch2Class;
%     ALLRESULTS(i).None = ~Ch1Class & ~Ch2Class;
% end
% 
% 
% % Creating plots
% 
% for i = 1:9
%     ALLRESULTS2(i).Ch1 = sum(ALLRESULTS(i).CH1); 
%     ALLRESULTS2(i).Ch2 = sum(ALLRESULTS(i).CH2); 
%     ALLRESULTS2(i).Both = sum(ALLRESULTS(i).Both); 
%     ALLRESULTS2(i).None = sum(ALLRESULTS(i).None);
% end
% 
% CTypes = {'Col2a Only', 'ColX Only', 'Both', 'Neither'};
% 
% 
% for i = 1:9
%     NucCount(i) = numel(ALLRESULTS(i).CH1);
% end
% 
% 
% 
% CH1COUNTS  = [ALLRESULTS2(:).Ch1];
% CH2COUNTS  = [ALLRESULTS2(:).Ch2];
% BOTHCOUNTS = [ALLRESULTS2(:).Both];
% NoneCOUNTS = [ALLRESULTS2(:).None];
% 
% Totals           = [CH1COUNTS; CH2COUNTS; BOTHCOUNTS; NoneCOUNTS]';
% NormalizedTotals = Totals./NucCount';
% NormalizedTotals = NormalizedTotals*100;
% 
% 
% figure
% ax = subplot(1, 2, 1);
% bar(NormalizedTotals, 'stacked')
% title('Protein Expression per Nucleus')
% ylabel('Fraction of Nuclei (%)')
% xlabel('Group')
% legend(CTypes)
% xticklabels(CNames)
% ax.FontSize = 14;
% 
% ax = subplot(1, 2, 2);
% bar(Totals, 'stacked')
% title('Protein Expression Across Image Volume')
% ylabel('Number of Nuclei (Counts)')
% xlabel('Group')
% legend(CTypes)
% xticklabels(CNames)
% ax.FontSize = 14;
% 
% %% Bar Plots with error bars
% 
% 
% for i = 1:9
%     % Getting average fluor
%     AVGEXP(i, 1) = mean([Ch1R(i).Fvox{:}], 'omitmissing');
% 
%     AVGEXP(i, 2) = mean([Ch2R(i).Fvox{:}], 'omitmissing');
%     AVGVAR(i, 1) = var([Ch1R(i).Fvox{:}], 'omitmissing');
%     AVGVAR(i, 2) = var([Ch2R(i).Fvox{:}], 'omitmissing');
% 
% 
%     AVGSTD(i, 1) = std([Ch1R(i).Fvox{:}], 'omitmissing');
%     AVGSTD(i, 2) = std([Ch2R(i).Fvox{:}], 'omitmissing');
% end
% 
% nPerGroup = 3;
% for i = 1:9
%     PerGroup(1)      = mean(AVGEXP(1:3, 1));
%     PerGroup(2)      = mean(AVGEXP(4:6, 1));
%     PerGroup(3)      = mean(AVGEXP(7:9, 1));
%     PerGroup(2,1)    = mean(AVGEXP(1:3, 2));
%     PerGroup(2,2)    = mean(AVGEXP(4:6, 2));
%     PerGroup(2,3)    = mean(AVGEXP(7:9, 2));
%     PerGroupSTD(1,1) = mean(AVGSTD(1:3, 1));
%     PerGroupSTD(1,2) = mean(AVGSTD(4:6, 1));
%     PerGroupSTD(1,3) = mean(AVGSTD(7:9, 1));
%     PerGroupSTD(2,1) = mean(AVGSTD(1:3, 2));
%     PerGroupSTD(2,2) = mean(AVGSTD(4:6, 2));
%     PerGroupSTD(2,3) = mean(AVGSTD(7:9, 2));
% end
% 
% % Plot
% BarErrPlot(PerGroup, PerGroupSTD)