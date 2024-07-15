%% Add the path with results

addpath('D:\MichaelPalmer_Gillis\2024_04_26 - ATDC5 Cells\IlluminationProfileCorrected\Results')

%% Select nuclear channel detections - POST DBSCAN 
load NuclearDetections.mat

%% Select channel of interest segmentation
load Segmentation_488Info.mat SegmentationResults
Ch1 = cat(2, SegmentationResults.SegmentationInfo);

load Segmentation_640Info.mat SegmentationResults
Ch2 = cat(2, SegmentationResults.SegmentationInfo);

% Clean up memory
clear SegmentationResults


%% 

FIDs = [{Ch1(1,:).FilePath} {Ch2(1,:).FilePath}];

% Finding file with z shift (missing z slices)
imds      = imageDatastore(FIDs);
FilesInfo = cellfun(@dir, imds.Files, 'UniformOutput', false);
Filesz    = cellfun(@(x) x.bytes, FilesInfo);

% Largest file - greatest indices
[~, fmx]  = max(Filesz);
fileInfo  = imfinfo(imds.Files{fmx});
Zmx       = numel(fileInfo);
Mmx       = fileInfo(1).Height;
Nmx       = fileInfo(1).Width;

% Smallest file - lesser indices
[~, fmn]  = min(Filesz);
fileInfo  = imfinfo(imds.Files{fmn});
Zmn       = numel(fileInfo);
Mmn       = fileInfo(1).Height;
Nmn       = fileInfo(1).Width;

% Indices to use for region cropping
y = [50 Mmn-50];  % The line scan has some artifacts along z at the ends
x = [1 Nmn];
z = [1 Zmn];

% y = [1 Mmn];




%% Classification of Ch1

% Perform analysis on images
cropregion = [x' y' z'];
Results_Ch1    = ClassifyExpression(Nuclei, Ch1, r, dx, dz, cropregion, []);

%% Classification of Ch2

% Perform analysis on images
Results_Ch2    = ClassifyExpression(Nuclei, Ch2, r, dx, dz, cropregion, []);

%% Condition names

CNames = {'Control 1', 'Control 2', 'Control 3', ...
          'Day 0 - Well 1', 'Day 0 - Well 2','Day 0 - Well 3',...
          'Day 2 - Well 1', 'Day 2 - Well 2','Day 2 - Well 3'};

for i = 1:9
    Ch1R(i).Expression = Results_Ch1{i}.Expression;
    Ch2R(i).Expression = Results_Ch2{i}.Expression;
    Ch1R(i).F          = Results_Ch1{i}.SumF;
    Ch2R(i).F          = Results_Ch2{i}.SumF;
    Ch1R(i).Fcell      = Results_Ch1{i}.FPerCellA;
    Ch2R(i).Fcell      = Results_Ch2{i}.FPerCellA;
    Ch1R(i).Fvol       = Results_Ch1{i}.FPerCellV;
    Ch2R(i).Fvol       = Results_Ch2{i}.FPerCellV;
    Ch1R(i).Fvox       = Results_Ch1{i}.FPerVox;
    Ch2R(i).Fvox       = Results_Ch2{i}.FPerVox;
end

%%

for i = 1:9
    AVGEXP(i, 1) = mean([Ch1R(i).Expression{:}]);
    AVGVAR(i, 1) = var([Ch1R(i).Expression{:}]);
    AVGEXP(i, 2) = mean([Ch2R(i).Expression{:}]);
    AVGVAR(i, 2) = var([Ch2R(i).Expression{:}]);
end

%% All
ALLRESULTS = [];

for i = 1:9
    Ch1Class           = logical([Ch1R(i).Expression{:}]');
    Ch2Class           = logical([Ch2R(i).Expression{:}]');
    ALLRESULTS(i).CH1  = Ch1Class & ~Ch2Class;
    ALLRESULTS(i).CH2  = Ch2Class & ~Ch1Class;
    ALLRESULTS(i).Both = Ch1Class & Ch2Class;
    ALLRESULTS(i).None = ~Ch1Class & ~Ch2Class;
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