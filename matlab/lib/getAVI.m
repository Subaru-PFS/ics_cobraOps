function output=getAVI(filename)
%% import video and difference frames

    vid = VideoReader(filename);
    if vid.Duration > 1
        disp('Video length > 1 sec');
        return
    end
    
    jj = 1;
    while hasFrame(vid)
        f(:,:,jj) = int16(readFrame(vid));
        jj=jj+1;
    end
    
    df = diff(f,1,3);
    keyboard;
    output=packstruct(f,df);
    
