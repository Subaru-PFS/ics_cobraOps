function FIDM = readFidFile(fidfile)

    FH = fopen(fidfile);
    A = fscanf(FH, 'raw[%f,%f], horned[%f,%f] ',[4,inf]);
    A = A.';
    
    % Approx X coord of fiducials to identify them by
    FIDX = [1644,1257,94];
    
    fid{1} = [];
    fid{2} = [];
    fid{3} = [];
    
    for ii=1:length(A(:,3))
        
        % Identify which Fidicual centroid belongs to
        [minD, fidI] = min(abs(FIDX - A(ii,3)));
        
        fid{fidI} = [fid{fidI} A(ii,3)+A(ii,4)*i];
    end
    
    try
        minVL = min([length(fid{1}),length(fid{2}),length(fid{3})]);
        fid{1} = fid{1}(1:minVL);
        fid{2} = fid{2}(1:minVL);
        fid{3} = fid{3}(1:minVL);
        FIDM = [fid{1}.',fid{2}.',fid{3}.'];
    catch err
        disp(err)
        keyboard;
    end
    
    figure
    plot(FIDM,'o')
    axis equal
        
    
end
        
        
        
        
        