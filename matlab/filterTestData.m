function output = filterTestData(pos, minSteps, minAngle, maxAngle, maxSpeed, filterNegMoves)
% minSteps -  minimum Number of commanded steps in data
% minAngle - minimum move angle [deg]
% maxAngle - maximum move angle [deg]
% maxSpeed - maximum motor map value [deg/step]


for id=1:numel(pos) % 1:the number of positioners used.

  fwdminSteps = pos(id).s1.steps > minSteps;
  rvsminSteps = pos(id).s1.steps < -minSteps;

  fwdmaxAngle = abs(pos(id).s1.moveSizes) < maxAngle;
  rvsmaxAngle = abs(pos(id).s1.moveSizes) < maxAngle;

  fwdminAngle = abs(pos(id).s1.moveSizes) > minAngle;
  rvsminAngle = abs(pos(id).s1.moveSizes) > minAngle;
  
  if filterNegMoves
  fwdmaxSpeed = pos(id).s1.mmap < maxSpeed  &  pos(id).s1.mmap > 0;
  rvsmaxSpeed = abs(pos(id).s1.mmap) < maxSpeed  &  pos(id).s1.mmap > 0;
  else 
  fwdmaxSpeed = pos(id).s1.mmap < maxSpeed; %  &  pos(id).s1.mmap > 0;
  rvsmaxSpeed = abs(pos(id).s1.mmap) < maxSpeed; %  &  pos(id).s1.mmap > 0;    
  end
  fwdmaxIsNNaN = ~isnan(pos(id).s1.moveSizes);
  rvsmaxIsNNaN = ~isnan(pos(id).s1.moveSizes);
  
  fwdmoves = fwdminSteps & fwdmaxAngle & fwdmaxSpeed & fwdmaxIsNNaN & fwdminAngle;
  rvsmoves = rvsminSteps & rvsmaxAngle & rvsmaxSpeed & rvsmaxIsNNaN & rvsminAngle; 
  
  % Filter forward Moves
  pos(id).s1f.moveSizes        = pos(id).s1.moveSizes(  fwdmoves);
  pos(id).s1f.steps            = pos(id).s1.steps(      fwdmoves);
  pos(id).s1f.startAngle       = pos(id).s1.startAngle( fwdmoves);
  pos(id).s1f.finishAngle      = pos(id).s1.finishAngle(fwdmoves);
  pos(id).s1f.mmap             = pos(id).s1.mmap(       fwdmoves);
  pos(id).s1f.iteration        = pos(id).s1.iteration(  fwdmoves);
  pos(id).s1f.targetNo         = pos(id).s1.targetNo(  fwdmoves);

  % Filter reverse Moves       
  pos(id).s1r.moveSizes        = pos(id).s1.moveSizes(  rvsmoves);
  pos(id).s1r.steps            = pos(id).s1.steps(      rvsmoves);
  pos(id).s1r.startAngle       = pos(id).s1.startAngle( rvsmoves);
  pos(id).s1r.finishAngle      = pos(id).s1.finishAngle(rvsmoves);
  pos(id).s1r.mmap             = pos(id).s1.mmap(       rvsmoves);
  pos(id).s1r.iteration        = pos(id).s1.iteration(  rvsmoves);
  pos(id).s1r.targetNo        = pos(id).s1.targetNo(  rvsmoves);


  fwdminSteps = pos(id).s2.steps > minSteps;
  rvsminSteps = pos(id).s2.steps < -minSteps;
  
  fwdmaxAngle = abs(pos(id).s2.moveSizes) < maxAngle;
  rvsmaxAngle = abs(pos(id).s2.moveSizes) < maxAngle;

  fwdminAngle = abs(pos(id).s2.moveSizes) > minAngle;
  rvsminAngle = abs(pos(id).s2.moveSizes) > minAngle;
  if filterNegMoves
  fwdmaxSpeed = pos(id).s2.mmap < maxSpeed &  pos(id).s2.mmap> 0; 
  rvsmaxSpeed = abs(pos(id).s2.mmap) < maxSpeed &  pos(id).s2.mmap > 0;
  else 
  fwdmaxSpeed = pos(id).s2.mmap < maxSpeed; % &  pos(id).s2.mmap> 0; 
  rvsmaxSpeed = abs(pos(id).s2.mmap) < maxSpeed; % &  pos(id).s2.mmap > 0;
  end 
  
  fwdmaxIsNNaN = ~isnan(pos(id).s2.moveSizes);
  rvsmaxIsNNaN = ~isnan(pos(id).s2.moveSizes);

  fwdmoves = fwdminSteps & fwdmaxAngle & fwdmaxSpeed & fwdmaxIsNNaN & fwdminAngle;
  rvsmoves = rvsminSteps & rvsmaxAngle & rvsmaxSpeed & rvsmaxIsNNaN & rvsminAngle;
  
  

  % Filter forward Moves
  pos(id).s2f.moveSizes        = pos(id).s2.moveSizes(fwdmoves);
  pos(id).s2f.steps            = pos(id).s2.steps(fwdmoves);
  pos(id).s2f.startAngle       = pos(id).s2.startAngle(fwdmoves);
  pos(id).s2f.finishAngle      = pos(id).s2.finishAngle(fwdmoves);
  pos(id).s2f.mmap             = pos(id).s2.mmap(fwdmoves);
  pos(id).s2f.iteration        = pos(id).s2.iteration(fwdmoves);
  pos(id).s2f.targetNo        = pos(id).s2.targetNo(fwdmoves);

  % Filter reverse Moves       
  pos(id).s2r.moveSizes        = pos(id).s2.moveSizes(rvsmoves);
  pos(id).s2r.steps            = pos(id).s2.steps(rvsmoves);
  pos(id).s2r.startAngle       = pos(id).s2.startAngle(rvsmoves);
  pos(id).s2r.finishAngle      = pos(id).s2.finishAngle(rvsmoves);
  pos(id).s2r.mmap             = pos(id).s2.mmap(rvsmoves);
  pos(id).s2r.iteration        = pos(id).s2.iteration(rvsmoves);
  pos(id).s2r.targetNo         = pos(id).s2.targetNo(rvsmoves);

  
end

output = pos;

end