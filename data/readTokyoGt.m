function [pos, t] = readTokyoGt(dataPath)

tb = readtable(dataPath);

t = tb.GPSTOW_s_;
x = tb.ECEFX_m_;
y = tb.ECEFY_m_;
z = tb.ECEFZ_m_;
pos = [x'; y'; z'];

end