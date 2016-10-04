function Jn = modify_damping_small(Jn)
  
  pp = polyfit(Jn.reqsml,Jn.err0sml,1); % line fit to get starting point
  % fit slope with y-int == 0;
  fitslope = fit(Jn.reqsml, Jn.err0sml, @(a,x) a*x, 'Startpoint', pp(1));
  Jn.errsml = Jn.err0sml - feval(fitslope,Jn.reqsml);
  Jn.slopesml = fitslope;
  
  %% pos and neg slopes
  pos = find(Jn.reqsml > 0);
  neg = find(Jn.reqsml < 0);
  posslope = fit(Jn.reqsml(pos), Jn.err0sml(pos), @(a,x) a*x, 'Startpoint', pp(1));
  if length(neg)>1
    negslope = fit(Jn.reqsml(neg), Jn.err0sml(neg), @(a,x) a*x, 'Startpoint', pp(1));
  else
    negslope.a = NaN;
  end
  
  Jn.posslopesml = posslope;
  Jn.negslopesml = negslope;
  
end