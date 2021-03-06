function plotme(refsol_path)
addpath(genpath(fullfile(pwd,'../../refsol/ode15s_v3/')))
addpath(genpath(fullfile(pwd,'../../../createfigures/')))
if nargin < 1
  refsol_path='';
end
[sr_err, sr_calls] = read_multirate_datapoints('BDF2_DEF.BDF2_DEF.3.GSCELLFirst_v3', 'slow', refsol_path);
[mr_err, mr_calls] = read_multirate_datapoints('BDF2_DEF.BDF2_DEF.3.GSFastFirst_v3', 'slow', refsol_path);
  
N = sr_calls(1,:) + sr_calls(2,:);
X1 = [ N; N; N]';
Y1 = sr_err';
N = mr_calls(1,:) + mr_calls(2,:);
X2 = [ N; N; N]';
Y2 = mr_err';
N = 10.^[5.1:0.1:5.9];
X3 = [ N; N]';
Y3 = [10^7*(10.^[5.1:0.1:5.9]).^(-2); 10^2*(10.^[5.1:0.1:5.9]).^(-1)]';
Createfigure(X1, Y1, X2, Y2, X3, Y3);
end
