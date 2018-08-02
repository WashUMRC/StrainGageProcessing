clear all;
% function strainGagingAnalysis()

answer = inputdlg('How many mice do you have to analyze for this strain gage group?');
num = str2num(answer{1});

answer = inputdlg('Please enter the gage factor for going from voltage to strain');
factor = str2num(answer{1});

answer = inputdlg('Do you need to apply any voltage offsets to any signals? y or n');
if strcmpi(answer{1},'y') == 1
    answer = inputdlg('Please enter signals to which an offset will be applied in the order of mouse,load,offset,mouse,load,offset,etc');
    offsets = num2str(answer{1});
else
    offsets = [];
end

answer = inputdlg('Do you think you will need to crop your signals? y or n');
cropYN = answer{1};

strainDiv       = 1000;     % Mystery division factor supplied by Michael.
excVoltage      = 3.33;  

for out = 1:num
    
    %User select folder first file, then it goes back for each additional
    %mouse
    if out == 1
        [forceFiles, forceDir ] = uigetfile([pwd '\*.*'],['Please select the force files of interest for mouse ' num2str(out)],'MultiSelect','on');
        [strainFiles, strainDir] = uigetfile([forceDir '\*.*'],['Please select the strain files of interest for mouse ' num2str(out)],'MultiSelect','on');
    else
        [forceFiles, forceDir ] = uigetfile([forceDir '\*.*'],['Please select the force files of interest for mouse ' num2str(out)],'MultiSelect','on');
        [strainFiles, strainDir] = uigetfile([strainDir '\*.*'],['Please select the strain files of interest for mouse ' num2str(out)],'MultiSelect','on');
    end
    
    %requires equal number of force and strain/voltage files
    if length(strainFiles) ~= length(forceFiles)
        error('Unequal number of files picked! Try again goober.')
    end
    
    %load in force and strain data. Strain files behave poorly, so used a
    %more robust reading technique
    for i = 1:length(forceFiles)
        forceData(i) = importdata([forceDir '\' forceFiles{i}],'\t',23);
        fid = fopen([strainDir '\' strainFiles{i}],'r');
        c=0;
        while ~feof(fid)
            c=c+1;
            in{c} = fgets(fid);
        end
        in = in(2:end);
        for j = 1:length(in)
            line = str2num(in{j});
            strainData(i).data(j,1) = line(1);
            strainData(i).data(j,2) = line(2);
        end
    end
    
    %Read in data and find mins
    for i = 1:length(forceData)
        %generate discreet vectors for readability
        force = forceData(i).data(:,3);
        if ~isempty(strcmpi(cropYN,'y'))
            plot(force)
            [x y] = ginput(2);
            force = force(x(1):x(2));
            close all;
        end
        forceTime = forceData(i).data(:,1);
        strain = strainData(i).data(:,2);
        strainTime = strainData(i).data(:,1);
        if ~isempty(strcmpi(cropYN,'y'))
            plot(strain)
            [x y] = ginput(2);
            strain = strain(x(1):x(2));
            close all;
        end
       
        %identify signal baseline between each force peak
%         [alarmForce] = envelop_hilbert_v2(force,1,1,30,1);
%         [alarmStrain] = envelop_hilbert_v2(strain,1,1,30,1);
        
        %find force range
        forceMin = min(force);
        forceMax = max(force);
        forceThresh = ((forceMin - forceMax)*0.3) + forceMax;%set floor for identifying peaks
        
        %find force mins
        [forceMins1 forceMinInds] = lmin(force,60);
        [forceMinsMinsInd] = find(forceMins1 < forceThresh);
        forceMins{out,i} = force(forceMinInds(forceMinsMinsInd));
        %discard first force min
        forceMins{out,i} = forceMins{out,i}(2:end);
        
        
        %find strain range
        strainMin = min(strain);
        strainMax = max(strain);
        strainThresh = ((strainMin - strainMax)*0.3) + strainMax;%set floor for identifying peaks
        
        %find strain mins
        [strainMins1 strainMinInds] = lmin(strain,60);
        [strainMinsMinsInd] = find(strainMins1 < strainThresh);
        strainMins{out,i} = strain(strainMinInds(strainMinsMinsInd));
        %discard first strain min
        strainMins{out,i} = strainMins{out,i}(2:end);
        
        %account for possible unequal lengths
        slen = length(strainMins{out,i});
        flen = length(forceMins{out,i});
        lenDif = slen-flen;
        if length(offsets) == 0;
            offset1 = 0;
        end
        for k = 1:3:length(offsets)
            if offsets(k) == i && offsets(k+1) == (j)
                offset1 = offsets(k+2);
            else
                offset1 = 0;
            end
        end
        smins = strainMins{out,i}(end-(end-2-lenDif):end-offset1); 
        fmins = forceMins{out,i}(end-(end-2):end); 
        
        sIndsAll = strainMinInds(strainMinsMinsInd);
        sInds = sIndsAll(end-(end-2-lenDif):end-offset1); 
        
        fIndsAll = forceMinInds(forceMinsMinsInd);
        fInds = fIndsAll(end-(end-2):end); 
        
        
        
        
        
        %correct force mins with current signal floor
        
        indsDiff = diff(fIndsAll);
        forceFloorInds = fIndsAll(1:end-1) + indsDiff ./ 2;
        clear forceCorrectionValues
        for j = 1:length(forceFloorInds)
            forceCorrectionValues(j) = mean(force(forceFloorInds(j)-50:forceFloorInds(j)+50));
        end
        for j = 1:length(forceCorrectionValues)
            forceMins{out,i}(j) = forceMins{out,i}(j) - forceCorrectionValues(j);
        end
        
        %correct strain mins with current signal floor
        
        indsDiff = diff(sIndsAll);
        strainFloorInds = sIndsAll(1:end-1) + indsDiff ./ 2;
        clear strainCorrectionValues
        for j = 1:length(strainFloorInds)
            strainCorrectionValues(j) = mean(strain(strainFloorInds(j)-50:strainFloorInds(j)+50));
        end
        for j = 1:length(strainCorrectionValues)
            strainMins{out,i}(j) = strainMins{out,i}(j) - strainCorrectionValues(j);
        end
        
    end
end

%check for equal numbers of forces and strains
for out = 1:length(strainMins(:,1))
    for i = 1:length(strainMins(out,:))
        strainLength = length(strainMins{out,i});
        forceLength = length(forceMins{out,i});
        if forceLength ~= strainLength
            answer = inputdlg(['Warning! Force and strain peak numbers are not equal for mouse ' num2str(out) ' load ' num2str(i) '. Proceed? y or n']);
            if strcmp(answer{1},'y') == 1
            elseif strcmp(answer{1},'n') == 1
                error('Check your original data for equal lengths');
            else
                error('Choose y or n like instructed next time');
            end
        else
            answer{1} = 'y';
        end

        if forceLength ~= 10 || strainLength ~= 10
            answer1 = inputdlg(['Warning! Force and strain peak numbers are not equal to 10 for mouse ' num2str(out) ' load ' num2str(i) '. Proceed? y or n']);
            if strcmp(answer1{1},'y') == 1
            elseif strcmp(answer1{1},'n') == 1
                error('Check your original data for number of peaks');
            else
                error('Choose y or n like instructed next time');
            end
        else
            answer1{1} = 'y';
        end
    end
end


smins = [];
fmins = [];
if strcmpi(answer1{1},'y') == 1 && strcmpi(answer{1},'y') == 1
    for i = 1:length(strainMins(:,1))
        for j = 1:length(strainMins(i,:))
            slen = length(strainMins{i,j});
            flen = length(forceMins{i,j});
            lenDif = slen-flen;
            if length(offsets) == 0;
                offset1 = 0;
            end
            for k = 1:3:length(offsets)
                if offsets(k) == i && offsets(k+1) == (j)
                    offset1 = offsets(k+2);
                else
                    offset1 = 0;
                end
            end
            smins = [smins;(strainMins{i,j}(end-(end-2-lenDif):end)- offset1)] ; 
            fmins = [fmins;forceMins{i,j}(end-(end-2):end)]; 
        end
    end
end


for i = 1:length(smins)
   smin1(i) = -4.*(((smins(i)./strainDiv))./excVoltage)./...
        (factor.*(1+2.*((smins(i)./strainDiv)-(smins(i)/strainDiv))./excVoltage));
end
smin1 = abs(smin1) .* 1000000;
fmins = abs(fmins);

fid = fopen([forceDir '\ForceStrainOutput.tab'],'w');
fprintf(fid,'%s\n','Remember, the output is in the order in which you selected your files');
for i = 1:length(smin1)
    if i == 1
        fprintf(fid,'%s\t','Microstrain');
        fprintf(fid,'%s\n','Force (N)');
    end
    fprintf(fid,'%s\t',num2str(smin1(i)));
    fprintf(fid,'%s\n',num2str(fmins(i)));
end
fclose(fid);
    
% function [lmval,indd]=lmin(xx,filt)
% %LMIN 	function [lmval,indd]=lmin(x,filt)
% %	Find local minima in vector X, where LMVAL is the output
% %	vector with minima values, INDD is the corresponding indeces 
% %	FILT is the number of passes of the small running average filter
% %	in order to get rid of small peaks.  Default value FILT =0 (no
% %	filtering). FILT in the range from 1 to 3 is usially sufficient to 
% %	remove most of a small peaks
% %	Examples:
% %	xx=0:0.01:35; y=sin(xx) + cos(xx ./3); 
% %	plot(xx,y); grid; hold on;
% %	[a b]=lmin(y,2)
% %	 plot(xx(a),y(a),'r+')
% %	see also LMAX, MAX, MIN
% 	
% %
% %**************************************************|
% % 	Serge Koptenko, Guigne International Ltd., |
% %	phone (709)895-3819, fax (709)895-3822     |
% %--------------06/03/97----------------------------|
% 
% x=xx;
% len_x = length(x);
% 	fltr=[1 1 1]/3;
%   if nargin <2, filt=0; 
% 	else
% x1=x(1); x2=x(len_x); 
% 
% 	for jj=1:filt,
% 	c=conv(fltr,x);
% 	x=c(2:len_x+1);
% 	x(1)=x1;  
%         x(len_x)=x2; 
% 	end
%   end
% 
% lmval=[];
% indd=[];
% i=2;		% start at second data point in time series
% 
%     while i < len_x-1,
% 	if x(i) < x(i-1)
% 	   if x(i) < x(i+1)	% definite min
% lmval =[lmval x(i)];
% indd = [ indd i];
% 
% 	   elseif x(i)==x(i+1)&x(i)==x(i+2)	% 'long' flat spot
% %lmval =[lmval x(i)];	%1   comment these two lines for strict case 
% %indd = [ indd i];	%2 when only  definite min included
% i = i + 2;  		% skip 2 points
% 
% 	   elseif x(i)==x(i+1)	% 'short' flat spot
% %lmval =[lmval x(i)];	%1   comment these two lines for strict case
% %indd = [ indd i];	%2 when only  definite min included
% i = i + 1;		% skip one point
% 	   end
% 	end
% 	i = i + 1;
%     end
% 
% if filt>0 & ~isempty(indd),
% 	if (indd(1)<= 3)|(indd(length(indd))+2>length(xx)), 
% 	   rng=1;	%check if index too close to the edge
% 	else rng=2;
% 	end
% 
% 	   for ii=1:length(indd), 
% 		[val(ii) iind(ii)] = min(xx(indd(ii) -rng:indd(ii) +rng));
% 		iind(ii)=indd(ii) + iind(ii)  -rng-1;
% 	   end
%   indd=iind; lmval=val;
% else
% end
% 
