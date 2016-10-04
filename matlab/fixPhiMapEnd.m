function output = fixPhiMapEnd(phiMotorMap,nsigma)

flds = fields(phiMotorMap);

for ii=1:length(flds)
    tmap = phiMotorMap.(flds{ii});
    goodI = [];
    
    for kk = 2:length(tmap(2,:));
        tstd = std(tmap(2,[1:kk]));
        if abs(tmap(2,kk)-tmap(2,kk-1)) > nsigma*tstd
            tmap2 = tmap(:,[1:kk-1]);
            break
        end
    end
    
    tmap2 = [tmap2 [200;tmap2(2,end)]];
    
    figure(100+ii)
    subplot(1,2,1)
    plot(tmap(1,:),tmap(2,:),'bx')
    hold on
    plot(tmap2(1,:),tmap2(2,:),'r')
    subplot(1,2,2)
    plot(tmap2(1,:),tmap2(2,:),'r')
    
    phiMotorMap.(flds{ii}) = tmap2;
    
    
end



output = phiMotorMap;

end
