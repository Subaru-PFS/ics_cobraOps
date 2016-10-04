function q=charge_state(data)
% calculate charge state of a given data set

  V_terminal = data.globals.accel.tps.GvmVR.data;
  V_ESA      = data.globals.he.esa.VC.data;
  q = round(30*(V_terminal + .01) / (V_ESA - 30*V_terminal));
