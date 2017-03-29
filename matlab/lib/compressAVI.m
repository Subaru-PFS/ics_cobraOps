function output=compressAVI(filename, threshold)
% detect movement in frames and eliminate no-movement frames
% quick look at a few frames indicates that min/max for noise is
% ~+-40 DN.

    if ~exist('threshold','var'), threshold = 50; end;
    if ~exist('compressed','dir'), mkdir('compressed'); end;

    % allow this many frames with no motion detection after last one.
    max_inactive_frames = 20;
    n_inactive_frames = 0;
    
    % crop image if necessary
    rowRange = 250:1124;
    colRange = 725:1524;
    mvpix = []; % defined to handle no-data case.
    maxpx = [];
    
    vidIN = VideoReader(filename);
    % if the images is already cropped, take all of the data.
    if (vidIN.Height == length(rowRange)) & (vidIN.Width == length(colRange))
        rowRange = 1:vidIN.Height;
        colRange = 1:vidIN.Width;
        outfilename = 'deleteme.avi';
    else
        outfilename = ['compressed/' filename(1:end-8) ...
                       filename(end-4) '.avi'];
        % skip existing files
        if exist(outfilename,'file')
           fprintf(1,'%s already exists... exiting...\n', ...
                   outfilename);
           return;
        end
    end
    fprintf(1,'%s -> %s\n',filename,outfilename);
    
    vidOUT= VideoWriter(outfilename,'Grayscale AVI');
    open(vidOUT);
    
    firstFrame = readFrame(vidIN);
    firstFrame = firstFrame(rowRange,colRange,1);
   
    lastFrame = firstFrame;
    
    NothingWrittenYet = true;
    WriteVideoFrame   = false;
    
    indx = 1;
    while hasFrame(vidIN);
        thisFrame = readFrame(vidIN);
        thisFrame = thisFrame(rowRange,colRange,1);

        delta1 = thisFrame - lastFrame;
        delta2 = lastFrame - thisFrame;
        delta = max(delta1,delta2);
        
        mvpix(indx) = length(find(delta > threshold));
        maxpx(indx) = max(reshape(delta,1,[]));
        
        %% write all frames once motion has been detected.
        if maxpx(indx) > threshold
            if NothingWrittenYet
                % 10 initial frames and 10 next-to last frames to
                % pause video at the initial condition.
                for kk =1:10
                    writeVideo(vidOUT,firstFrame);
                end
                for kk =1:10
                    writeVideo(vidOUT,lastFrame);
                end
                NothingWrittenYet = false; 
            end
            WriteVideoFrame = true;
            n_inactive_frames = 0;
        else
            n_inactive_frames = n_inactive_frames + 1;
            if n_inactive_frames > max_inactive_frames
                WriteVideoFrame = false;
            end
        end
        
        if WriteVideoFrame
            writeVideo(vidOUT,thisFrame);
        end
        
        indx = indx+1;
        lastFrame = thisFrame;
    end
    close(vidOUT)
    
    DELTA = abs(int16(lastFrame) - int16(firstFrame));
    if length(find(DELTA > threshold)) == 0
        fprintf(1,'No motion detected in %s\n',filename);
    end

    
    output = packstruct(DELTA,mvpix,maxpx);
    
% $$$     if ~WriteVideoFrame
% $$$         move_detected = find(maxpx > threshold);
% $$$         figure;
% $$$         plot(maxpx); hold on;
% $$$         plot(move_detected,maxpx(move_detected),'ro');
% $$$         hold off;
% $$$         refline(0,threshold,0,'r--');
% $$$         xlabel('frame number');
% $$$         ylabel('max delta signal between frames');
% $$$         title(strrep(filename,'_','\_'));
% $$$         grid on;
% $$$     end