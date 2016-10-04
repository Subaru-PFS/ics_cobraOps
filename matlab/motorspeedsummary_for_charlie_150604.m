digits_precision = 2;
Fast_OR_Slow = 'FAST';


cfg = loadCfgXml;



disp(cfg.cfgFile);
disp('RP   F_1_fwd   F_1_rev   F_2_fwd   F_2_rev');
disp('------------------------------------------');
for railID = 1:27
  try
    [F1f x] = getMMap(cfg,railID,Fast_OR_Slow,1,'fwd');
    [F1r x] = getMMap(cfg,railID,Fast_OR_Slow,1,'rev');
    [F2f x] = getMMap(cfg,railID,Fast_OR_Slow,2,'fwd');
    [F2r x] = getMMap(cfg,railID,Fast_OR_Slow,2,'rev');
    fprintf(1,'%2d:  %s %s %s %s\n',railID,...
            format_data(mean(F1f(3:end)),std(F1f(3:end)),digits_precision),...
            format_data(mean(F1r(3:end)),std(F1r(3:end)),digits_precision),...
            format_data(mean(F2f(3:end)),std(F2f(3:end)),digits_precision),...
            format_data(mean(F2r(3:end)),std(F2r(3:end)),digits_precision) ...
            );
  end
end

clear all;