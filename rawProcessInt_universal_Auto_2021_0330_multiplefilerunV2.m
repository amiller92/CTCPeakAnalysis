% Last edit on 3/30/2021 by Alex Miller

close all
clear all
clc
clf

tic
% leave '\' at the end of the directory line and ',' between each
% directory
directory = {'D:\Data\2021\2021-03-30 BE NSCLC\TBM MB10985\BE\',...      
     '\\rowley.mit.edu\manalis\abmiller\_SORTER 2 2021\2021-03-30 BE NSCLC\HM MB11389 BE\'};

% paste the .bin file between '  ' marks, separated by ','
filename1 = {'2021-03-30_1108_PMT',...
      '2021-03-30_1047_PMT'};
   

% confirm that directory and filename have the same lengths
if length(directory) ~= length(filename1)
    error('check directories and filenames-- unequal lengths')
end

for direc = 1:length(directory)
    
    filename = strcat(directory(direc),filename1(direc));
    filename = cellstr(filename);
    
    figDir = char(strcat(directory(direc),'PeakFigures_Auto',strrep(datestr(now),':','_'),'\'));
    mkdir(figDir)
    
    % format
    % parameters
    %mainThreshold    = 0.05;
    
    filtering    = 1;           % low pass filtering
    medianFilter = 1;           % subtract median (window by window)
    
    chunkSize    = 2000;       % size (in datapoints) for analysis window
    dataRate     = 30000;      % acquisition data rate
    
    limitMaxPeakWidth = 35;
    limitMinPeakWidth = 3;
    limitMinPeakSeparation = 35;
    limitMaxPeakSeparation = 600; % 230 equals minimum speed of 30 mm/sec between
    % two laser lines 460 um apart
    limitConsecPeaks = 4;
    
    minPeakRatio = 0.01;
    maxPeakRatio = 1/minPeakRatio;
    
    SNR = 6; %for line 125 to set the theshold
    a_ratio = 5; %the smaller, the narrower the peak is
    
    % output
    output = 1;
    
    %%%% get PMT data
    data = [];
    for i=1:1:length(filename)
        fileID = fopen([filename{i} '.bin']);
        data_temp = fread(fileID,Inf,'uint16',0,'b');
        data = [data;data_temp];
    end
    
    data = data*5/65536;
    
    fileSize = length(data);
    disp(['File size: ' num2str(fileSize) ' data points'])
    
    flength_sec = fileSize/dataRate;
    x_start_of_analysis = flength_sec-3600;
    
    %%%% get file creation timestamp
    listing = dir([filename{1} '.bin']);
    assert(numel(listing) == 1, 'No such file: %s', [filename '.bin']);
    startTime = listing.datenum - fileSize/dataRate/86400;
    finishTime = listing.datenum;
    disp(['File start time: ' datestr(startTime, 'HH:MM:SS')])
    disp(['File end time: ' datestr(finishTime, 'HH:MM:SS')])
    disp(['File length: ' datestr(fileSize/dataRate/86400, 'HH:MM:SS')])
    
    % keep original data
    originalData = data;
    
    
    %%%% setup auxiliary variables
    peaks = zeros(2000000,2);
    peaksCounter = 0;
    plotCount = 0;
    TotalPeaks = 0;
    
    %%%% To view the figures in a specific corner of your screen, type
    % figure
    %%then move plot to the designated area then type to know its position
    %p = get(gcf, 'Position')
    %%finally set that position declaring
    %p = [-1326, 620, 1063, 594];
    %set(0,'DefaultFigurePosition',p);
    
    p = [-1378, 309, 1063, 594];
    set(0,'DefaultFigurePosition',p);
    peak_time = [];
    
    %%%% SPLIT data into chunks in order to not overload memory.
    %   run through all of the sections
    last_n = floor(length(data)/chunkSize)-1;
    
    for n = 0:1:floor(length(data)/chunkSize)-1
        clf
        %loop_percent = n/last_n*100
        peaks = [];
        peakCounter = 0;
        congregate = [];
        cong_ind = [];
        betweenPeak = [];
        peakIndex = [];
        peaksLocs = [];
        peakROI = [];
        percentAnalyzed = floor((chunkSize*n+1)/length(data)*100);
        %disp(['Percent Competed: ' num2str(percentAnalyzed) '%'])
        dataROI = originalData(chunkSize*n+1:chunkSize*(n+1)); %original data
        dataROI_f = originalData(chunkSize*n+1:chunkSize*(n+1)); %filtered data
        
        % apply a median filter to remove AC component
        if medianFilter == 1
            med = median(dataROI);
            dataROI_f = dataROI_f - median(dataROI_f);
            dataROI = dataROI - median(dataROI);
        end
        
        [bb_filt,aa_filt] = cheby1(2,0.05,.07,'low');
        dataROI_f = filter(bb_filt,aa_filt,dataROI_f); %filtered
        format long;
        %threshold = mean(dataROI)+ mainThreshold;
        mainThreshold = mean(dataROI)+0.04;
        noise = std(dataROI_f(dataROI_f<mainThreshold));
        threshold = noise*SNR;
        %-------------------- find points above threshold--------------------------
        locs = find(dataROI_f>threshold);
        
        % % -------- find data behavior below threshold (if it's bumpy)--------------
        %     % the idea here is that if the baseline is too bumpy, then exclude
        %     % the peaks within this window and move to the next iteration
        %     locsBelow = find(dataROI<threshold);
        %     maxBelow = max(dataROI(locsBelow));
        %     minBelow = min(dataROI(locsBelow));
        
        % ----- check how bumpy this window is (removal of bubbles or debris)------
        belowZ = length(find(dataROI<0));
        aboveZ = length(find(dataROI>0 & dataROI<threshold));
        
        if abs(belowZ-aboveZ)>300
            continue
        end
        %     if (maxBelow-minBelow) > 2.5*threshold
        %         continue
        %     end
        
        % -------------------------------------------------------------------------
        if ~isempty(locs)
            plotCount = plotCount + 1;
            
            % monitor the data
            xlabel('Time (s)')
            ylabel('PMT Arbitrary Intensity');
            ylim([-0.2 0.5]);
            
            plot((n*chunkSize+linspace(0,chunkSize-1,chunkSize))/dataRate,dataROI_f,'-.b')
            hold on
            grey = [0.8, 0.8, 0.8];
            plot((n*chunkSize+linspace(0,chunkSize-1,chunkSize))/dataRate,dataROI,'Color',grey)
            plot(([1 chunkSize]+chunkSize*n)/dataRate,[threshold,threshold], '-r')
            plot(([chunkSize chunkSize].*n)/dataRate,[0,.1], ':g')
            
            % identify continuous (within limitNoise) regions of points above threshold
            locsOI = diff(locs,1);
            locsOI(locsOI<2) = 1; %so any peak that goes over and under thershold line
            % within few data points is disregarded
            
            % but before remember singlets
            singleLocs = find(locsOI>limitMaxPeakWidth);
            singleLocs = singleLocs(diff(singleLocs)==1)+1;
            
            locsOI(locsOI>1) = 0;
            
            % get edges of those regions
            edgesLocsOI = diff([0;locsOI;0]);
            edgesLocsOI = [find(edgesLocsOI==1) find(edgesLocsOI==-1)];
            
            % put back singlets
            edgesLocsOI = sort([edgesLocsOI;singleLocs singleLocs]);
            
            % check the width of each peak to make sure it is not too big
            if ~isempty(edgesLocsOI)
                peakWidthCheck = edgesLocsOI(:,2)-edgesLocsOI(:,1);
                tooBig = find (peakWidthCheck > limitMaxPeakWidth);
                tooSmall = find (peakWidthCheck < limitMinPeakWidth);
                tooBig_tooSmall = [tooBig;tooSmall];
                edgesLocsOI (tooBig_tooSmall,:)=[];
            end
            
            % go trough each region and determine maximum value and its index
            if ~isempty(edgesLocsOI)
                edgesSize = size(edgesLocsOI);
                for i=1:1:edgesSize(1)
                    peakLocs = locs(edgesLocsOI(i,1):edgesLocsOI(i,2));
                    % check three extra data points before and after to make
                    % sure that the maximum of the raw data is included
                    if peakLocs(1)>3 && peakLocs(end)<chunkSize-3
                        peakLocs = [peakLocs(1)-3;peakLocs(1)-2;peakLocs(1)-1;...
                            peakLocs;peakLocs(end)+1;peakLocs(end)+2;peakLocs(3)+3];
                    end
                    peakROI = dataROI(peakLocs);
                    peakROI_f = dataROI_f(peakLocs);
                    [peakVal,peakIndex] = max(peakROI);
                    [peakVal_f,peakIndex_f] = max(peakROI_f);
                    peakWidth = length(find(peakROI_f>threshold));
                    %determine peak height to width
                    peakIndex = peakLocs(peakIndex)+n*chunkSize;
                    peakIndex_f = peakLocs(peakIndex_f)+n*chunkSize;
                    %----------------------------------------------------------
                    %--------------------CHECK PEAK HEIGHT TO WIDTH------------
                    if (peakVal_f/peakWidth)>(threshold/a_ratio)
                        peaksCounter = peaksCounter + 1;
                        peaks(peaksCounter,:) = [peakIndex peakVal];
                        
                        plot((peakIndex-1)/dataRate,peakVal,'.m')
                        plot((peakIndex_f-1)/dataRate,peakVal_f,'.m')
                        %---------------------------------------------
                        %---------------------------------------------
                    end
                end
            end
            
            
            % remove empty slots from pre-allocation
            peaks(peaksCounter+1:end,:) = [];
            %         peaks
            %         size(peaks,1)
            %
            %-------------------- more than 3 congregating peaks ----------------------
            if  ~isempty(peaks) && size(peaks,1)>=3
                congregate = peaks(3:end,1)-peaks(1:end-2,1);
                cong_ind = find(congregate<150);
                cong_ind = union(cong_ind,cong_ind+2);
                if ~isempty(cong_ind)
                    peaks(cong_ind,:) = [];
                    
                end
            end
            
            if ~isempty(peaks)
                % --------------- remove peaks that are too far apart ---------------------
                peaks(:,3) = [diff(peaks(:,1));10000];
                pairedPeaks = find(peaks(:,3)<limitMaxPeakSeparation);
                totalpairedPeaks = length(pairedPeaks);
                
                %         % PLOTTING TIMES AND OVERALL IMAGE SETUP
                %         time = (1:1:length(data))./dataRate;
                %         originalTime = (1:1:length(originalData))./dataRate;
                
                selectedPeaks = zeros(totalpairedPeaks,6);
                
                %-- eliminate improperly paired peaks | based on user-defined peak ratios--
                i2 = 1;
                for i=1:1:totalpairedPeaks
                    n2 = pairedPeaks(i);
                    if (i==1) || (peaks(n2,1) ~= selectedPeaks(max(i2-1,1),3))
                        peakRatio = peaks(n2,2)/peaks(n2+1,2);
                        if (peakRatio > minPeakRatio) && (peakRatio < maxPeakRatio)
                            selectedPeaks(i2,:) = [peaks(n2,1) peaks(n2,2) peaks(n2+1,1) peaks(n2+1,2) peaks(n2+1,1)-peaks(n2,1) peakRatio];
                            i2 = i2 + 1;
                        end
                    end
                end
                
                % trim empty elements, leftovers from pre-allocated variable
                selectedPeaks(i2:end,:) = [];
                
                %------------- remove peaks that are too close together--------------------
                selectedPeaks(selectedPeaks(:,5)<limitMinPeakSeparation,:) = [];
                
                %------- final check is the data stability between the two peaks ----------
                finalPeakCount = size(selectedPeaks);
                finalPeakCount = finalPeakCount(1);
                betweenPeaks = zeros(finalPeakCount,1);
                if ~isempty(selectedPeaks)
                    for i=1:finalPeakCount
                        betweenPeaks(i) = mean(data(selectedPeaks(i,1)+5:selectedPeaks(i,3)-5)-med);
                    end
                    
                    % ------- Select only data with mean intensity between the two peaks ------
                    % below some specific value - this varies between
                    % experiments
                    finalLocs = find(betweenPeaks<threshold);
                    
                    if ~isempty(finalLocs)
                        selectedPeaks = selectedPeaks(finalLocs,:);
                        text(selectedPeaks(:,1)./dataRate,selectedPeaks(:,2),'\rightarrow')
                        text(selectedPeaks(:,3)./dataRate,selectedPeaks(:,4),'\leftarrow')
                        x=1;
                    else
                        continue
                    end
                else
                    x = 0;
                end
                
                if x==1
                    disp(['Percent Competed: ' num2str(percentAnalyzed) '%'])
                    n
                    peaksWithin = size(finalLocs,1);
                    if peaksWithin > 1
                        TotalPeaks = (TotalPeaks+1:TotalPeaks+peaksWithin)'
                    else
                        TotalPeaks = TotalPeaks+1
                    end
                    
                    peak_time(TotalPeaks,1) = TotalPeaks;
                    peak_time(TotalPeaks,2) = selectedPeaks(:,1)./dataRate;
                    peak_time(TotalPeaks,3) = (selectedPeaks(:,1)./dataRate)./60;
                    if peaksWithin > 1
                        p = zeros(peaksWithin,1);
                        p(1)=peaksWithin;
                        peak_time(TotalPeaks,4) = p;
                    else
                        peak_time(TotalPeaks,4) = peaksWithin;
                    end
                    peak_time(TotalPeaks,5) = n;
                    selectedPeaksSize = size(selectedPeaks,1);
                    peak_time(TotalPeaks,6:11) = selectedPeaks;
                    TotalPeaks = TotalPeaks(end);
                    fname = ['Fig',num2str(TotalPeaks)];
                    saveas(gcf, strcat(figDir,fname),'fig');
                    fname = [];
                end
                hold off
            end
            
        end
        
    end
    
    if(output == 1)
        output_filename = 'peak_timesAuto';
        %output_filename = char(strcat(directory(direc),output_filename));
        output_filename = char(strcat(directory(direc),output_filename,'_',strrep(datestr(now),':','_'),'.csv'));
        dlmwrite(output_filename,  'count, time (sec), time (min), No. of Peaks Within, n, P1_Loc, P1_Height, P2_Loc, P2_Height, Distance, P2_P1 Ratio','');
        dlmwrite(output_filename, peak_time, '-append', 'precision', 14);
    end
end
toc