function output = placeArrows(filename, threshold)
% Place white arrows next to the Cobras that had movements
% filename is .avi
% Size of the avi should be 850 x 875
% Author: Mitsuko Roberts
    
    vidIN = VideoReader(filename);
    
    vidOUT= VideoWriter([filename(1:end-4) '_arrow.avi'],'Grayscale AVI');
    open(vidOUT);
    
    vidFrames = read(vidIN);
    vidFrames = reshape (vidFrames(:,:,1,:),vidIN.Height,vidIN.Width,vidIN.NumberOfFrames);
    delta = zeros(vidIN.Height,vidIN.Width,vidIN.NumberOfFrames);

    
    % Print warning if video size is different
    
    if vidIN.Width~=850
        fprintf('Video width not 850\n');
    end
    
    if vidIN.Height~=875
        fprintf('Video Height not 875\n');
    end

    
    % Calculate difference of first and last frames (delta3)
    
    delta1 = vidFrames(:,:,1) - vidFrames(:,:,end);
    delta2 = vidFrames(:,:,end) - vidFrames(:,:,1);
    delta3 = max(delta1, delta2);
    
    [row,col] = find(delta3>threshold);
    
    % how many pixels are over threshold per region of interest - this
    % region might need to get adjusted 
    
    for i=1:13
    pid(i)=length(find(delta3(100+52*(i-1):150+52*(i-1),...
        50+52*(i-1):100+52*(i-1))>threshold));
    end
    
    
    stillFrame = vidFrames(:,:,end);
    
    % insert white arrows if more than 2 pixels had difference over
    % threshold
    for j=1:13
        if pid(j)>2
            stillFrame(125+52*(j-1):130+52*(j-1),10+52*(j-1):40+52*(j-1))=255;
            for i=1:5
                stillFrame(130+i+52*(j-1),35-i+52*(j-1):40-i+52*(j-1))=255;
                stillFrame(125-i+52*(j-1),35-i+52*(j-1):40-i+52*(j-1))=255;
            end
        end
    end
    
    figure;
    imshow(insertMarker(vidFrames(:,:,end),[col,row],'x'))
    figure;
    imshow(stillFrame)
    % insert command to write video

    close(vidOUT)
    
