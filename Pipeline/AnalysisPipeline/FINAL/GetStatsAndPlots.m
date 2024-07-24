function [Stats, ALLRESULTS] = GetStatsAndPlots(ExpressionResults, titlemod, ax1, ax2, xlab, ylab)

    if nargin < 6 || isempty(ylab)
        ylab = false;
    end

    if nargin < 5 || isempty(xlab)
        xlab = false;
    end

    if nargin < 2
        titlemod = [];
    end

    TrialNumber = ExpressionResults.Trial;
    Param1      = ExpressionResults.Param1;
    Param2      = ExpressionResults.Param2;
    
    CNames = {'Control 1', 'Control 2', 'Control 3', ...
              'Day 0 - Well 1', 'Day 0 - Well 2','Day 0 - Well 3',...
              'Day 2 - Well 1', 'Day 2 - Well 2','Day 2 - Well 3'};
    
    nfiles    = size(ExpressionResults.Results, 1);
    nchannels = size(ExpressionResults.Results, 2);
    Results   = ExpressionResults.Results;
   
    % Getting the fields
    fields  = fieldnames(Results);
    nfields = numel(fields);

    for f_idx = 1:nfields
        field = fields{f_idx};
        for ch = 1:nchannels
            for f = 1:nfiles
                % Grabbing the field and checking for valid data
                Val                     = [Results(f, ch).(field)];
                nodata                  = all(cellfun(@isempty, Val));
                if nodata
                    % Continue to next field when the field has no data
                    continue
                end

                % New field includes the name of the stat
                newfield                = [field 'Mean'];
                Stats(f, ch).(newfield) = GetMeasurementStat(Val, 'mean');
                newfield                = [field 'Var'];
                Stats(f, ch).(newfield) = GetMeasurementStat(Val, 'var');
            end
        end
    end
    % There might be an issue where some values are somehow going to INF
    % in the fluorescence quantification... or go to NaN ... figure out
    % whats up

%% Everything below this has to be adjusted to work on the above stats struct
    ALLRESULTS = [];
    
    for f = 1:9
        Ch1Class           = logical([Results(f, 1).Expression{:}]');
        Ch2Class           = logical([Results(f, 2).Expression{:}]');
        ALLRESULTS(f).CH1  = Ch1Class & ~Ch2Class;
        ALLRESULTS(f).CH2  = Ch2Class & ~Ch1Class;
        ALLRESULTS(f).Both = Ch1Class & Ch2Class;
        ALLRESULTS(f).None = ~Ch1Class & ~Ch2Class;
    end
    
    
    % Creating plots
    
    for f = 1:9
        ALLRESULTS2(f).Ch1 = sum(ALLRESULTS(f).CH1); 
        ALLRESULTS2(f).Ch2 = sum(ALLRESULTS(f).CH2); 
        ALLRESULTS2(f).Both = sum(ALLRESULTS(f).Both); 
        ALLRESULTS2(f).None = sum(ALLRESULTS(f).None);
    end
    
    CTypes = {'Col2a Only', 'ColX Only', 'Both', 'Neither'};
    
    
    for f = 1:9
        NucCount(f) = numel(ALLRESULTS(f).CH1);
    end
    

    % Modifying the title to have more information about the plots
    if ~isempty(titlemod)
        titlemod = [titlemod ': '];
    end

    titlemod = [titlemod newline...
                'P1: ' num2str(Param1) 'P2: ' num2str(Param2) newline];

    % First set of stats to show in bar plot
    CH1COUNTS  = [ALLRESULTS2(:).Ch1];
    CH2COUNTS  = [ALLRESULTS2(:).Ch2];
    BOTHCOUNTS = [ALLRESULTS2(:).Both];
    NoneCOUNTS = [ALLRESULTS2(:).None];
    Totals = [CH1COUNTS; CH2COUNTS; BOTHCOUNTS; NoneCOUNTS]';
    Totals = Totals./NucCount';
    Totals = Totals*100;


    % Setting up figure window
    if isempty(ax1)
        fig1 = figure;
        ax1  = axes(fig1);
    end

    % Plotting and labels
    bar(ax1, Totals, 'stacked')
    % ftitle = [titlemod  'Cell Count Normalized Protein Expression'];

    ftitle = titlemod;

    title(ax1, ftitle)

    if ylab
        ylabel(ax1, 'Fraction of Nuclei (%)')
    else
        yticks(ax2, [])
    end


    if xlab
        xlabel(ax1, 'Group')
        xticklabels(ax1, CNames)
    else
        xticks(ax1, [])
    end

    % Setting up legends and bar labels
    % legend(ax1, CTypes)
    

    % Adjusting the font size
    ax1.FontSize = 12;
    
    % Non normalized totals
    Totals = [CH1COUNTS; CH2COUNTS; BOTHCOUNTS; NoneCOUNTS]';

    % Setting up figure window
    if isempty(ax1)
        fig2 = figure;
        ax2  = axes(fig2);
    end

    % Plotting and labels
    bar(ax2, Totals, 'stacked')
    % ftitle = [titlemod  'Protein Expression'];


    title(ax2, ftitle)

    if ylab
        ylabel(ax2, 'Number of Nuclei (Counts)')
    else
        yticks(ax2, [])
    end

    if xlab
        xlabel(ax2, 'Group')
        xticklabels(ax2, CNames)
    else
        xticks(ax2, [])
    end

    % Adjusting the font size
    ax2.FontSize = 12;
    
    % Setting up legends and bar labels
    % legend(ax2, CTypes)
    
end

function Stat = GetMeasurementStat(V, stype)
% GETMEASUREMENTSTAT will get a statistical measurement of a field entry in
% a data struct such that it accounts for cell arrays as well as numerical
% arrays.

    switch stype
        case 'mean'
            if iscell(V)
                [dimchoice, V] = AdaptiveStat(V);
                Stat           = cellfun(@(x) mean(x, dimchoice, 'omitmissing'), V, 'UniformOutput', false);
            else
                Stat           = mean(V, 1, 'omitmissing');
            end
        case 'var'
            if iscell(V)
                [dimchoice, V] = AdaptiveStat(V);
                Stat           = cellfun(@(x) var(x, 1, dimchoice, 'omitmissing'), V, 'UniformOutput', false);
            else
                Stat           = var(V, 1, 'omitmissing');
            end
    end

    function [dimchoice, v] = AdaptiveStat(v)
        % ADAPTIVESTAT will adapt how the mean is computed depending on
        % whether the field entry has a numerical array, cell array, or a
        % nested cell array with collections of indices/points
        v_empty   = cellfun(@isempty, v);
        v_check   = v(~v_empty);
        sz        = cellfun(@(x) size(x), v_check, 'UniformOutput', false);
        hascoords = cellfun(@(x) x(1)==1 & x(2)==3, sz);
        if all(hascoords)
            dimchoice = 'all';
            v         = cat(1, v{:});
        else
            dimchoice = 1;
        end
    end
end