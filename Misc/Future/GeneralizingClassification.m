CNames = {'Control 1', 'Control 2', 'Control 3', ...
          'Day 0 - Well 1', 'Day 0 - Well 2','Day 0 - Well 3',...
          'Day 2 - Well 1', 'Day 2 - Well 2','Day 2 - Well 3'};

% Init
nFiles    = numel(CNames);
nChannels = 2;
RESULTS   = cell(nFiles, nChannels);


for j = 1:nChannels
    RESULTS{:, j} = ClassifyExpression(Nuclei, Ch1, r, dx, dz, cropregion, []);
    for i = 1:nFiles
        ExpResults(i,j).Expression = RESULTS{i, j}.Expression;
        ExpResults(i,j).F          = RESULTS{i, j}.SumF;
        ExpResults(i,j).Fcell      = RESULTS{i, j}.FPerCellA;
        ExpResults(i,j).Fvol       = RESULTS{i, j}.FPerCellV;
        ExpResults(i,j).Fvox       = RESULTS{i, j}.FPerVox;

        AVGEXP(i, j) = mean([ExpResults(i, j).Expression{:}]);
        AVGVAR(i, j) = var([ExpResults(i, j).Expression{:}]);
    end
end

%% All


ClassifyExpression(Nuclei, Ch1, r, dx, dz, cropregion, []);

ALLRESULTS  = [];
nChannels   = size(classification, 2);
allchannels = 1:nChannels;

Classification = cell(nFiles, nChannels);

for j = 1:nChannels
    for i = 1:9
        Classification{i,j}    = logical([ExpResults(i, j).Expression{:}]');
        chidx                  = allchannels(allchannels~=j);
        chname                 = ['Ch' num2str(j)];
        ALLRESULTS(i).(chname) = Classification{i,j} & ~Classification{i,chidx};
        ALLRESULTS(i).CH2   = Ch2Class & ~Classification;
        ALLRESULTS(i).Both  = Classification & Ch2Class;
        ALLRESULTS(i).None  = ~Classification & ~Ch2Class;
    end
end

% Creating plots

for i = 1:9
    ALLRESULTS2(i).Ch1 = sum(ALLRESULTS(i).CH1); 
    ALLRESULTS2(i).Ch2 = sum(ALLRESULTS(i).CH2); 
    ALLRESULTS2(i).Both = sum(ALLRESULTS(i).Both); 
    ALLRESULTS2(i).None = sum(ALLRESULTS(i).None);
end

CTypes = {'Col2a Only', 'ColX Only', 'Both', 'Neither'};


for i = 1:9
    NucCount(i) = numel(ALLRESULTS(i).CH1);
end


figure
CH1COUNTS  = [ALLRESULTS2(:).Ch1];
CH2COUNTS  = [ALLRESULTS2(:).Ch2];
BOTHCOUNTS = [ALLRESULTS2(:).Both];
NoneCOUNTS = [ALLRESULTS2(:).None];

Totals = [CH1COUNTS; CH2COUNTS; BOTHCOUNTS; NoneCOUNTS]';
Totals = Totals./NucCount';
Totals = Totals*100;
bar(Totals, 'stacked')
title('Protein Expression')
ylabel('Fraction of Nuclei (%)')

xlabel('Group')
ax = gca;
ax.FontSize = 14;

legend(CTypes)
xticklabels(CNames)



figure
Totals = [CH1COUNTS; CH2COUNTS; BOTHCOUNTS; NoneCOUNTS]';

bar(Totals, 'stacked')
title('Nuclei')
ylabel('Number of Nuclei (Counts)')

xlabel('Group')
ax = gca;
ax.FontSize = 14;

legend(CTypes)
xticklabels(CNames)