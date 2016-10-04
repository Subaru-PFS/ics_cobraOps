function label_out = autolabel(struct)
% create x/y/z label based on "label" or "info" and "unit" from structure
% USAGE: [xyz]label(autolabel(data.field))

  if isfield(struct,'label')
    label_out = struct.label;
  elseif isfield(struct,'info')
    label_out = struct.info;
  else
    warning('AUTOLABEL: no appropriate field available');
    label_out = '';
  end
  
  if isfield(struct,'unit')
    label_out = [label_out ' [' struct.unit ']'];
  else
    warning('AUTOLABEL: no unit field');
  end
  