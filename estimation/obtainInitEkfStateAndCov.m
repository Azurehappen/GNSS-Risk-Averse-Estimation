function [state, cov] = obtainInitEkfStateAndCov(p, estState)

% position, velocity, acceleration,
% clock biases (m), clock drift (m/s)
if p.state_mode == p.pva_mode
    state = [estState.pos;estState.vel;zeros(3, 1)];
    cov_diag = p.ekf_para.q_pos * ones(1,9);
elseif p.state_mode == p.pos_mode
    state = [estState.pos;estState.clock_bias;0];
    cov_diag = p.ekf_para.q_pos * ones(1,3);
end

if ~isnan(estState.clock_sys(p.gps.sys_num))
    state = [state; estState.clock_sys(p.gps.sys_num)];
    cov_diag = [cov_diag, 40^2];
end
if ~isnan(estState.clock_sys(p.glo.sys_num))
    state = [state; estState.clock_sys(p.glo.sys_num)];
    cov_diag = [cov_diag, 40^2];
end
if ~isnan(estState.clock_sys(p.gal.sys_num))
    state = [state; estState.clock_sys(p.gal.sys_num)];
    cov_diag = [cov_diag, 40^2];
end
if ~isnan(estState.clock_sys(p.bds.sys_num))
    state = [state; estState.clock_sys(p.bds.sys_num)];
    cov_diag = [cov_diag, 40^2];
end

state = [state; estState.clock_drift];
cov_diag = [cov_diag, 5^2];
cov = diag(cov_diag);

% if length(state) == size(estState.state_cov, 1)+3
%     cov(1:6,1:6) = estState.state_cov(1:6,1:6);
%     cov(7:9,7:9) = p.ekf_para.q_accHor*eye(3);
%     cov(1:6,10:end) = estState.state_cov(1:6,7:end);
%     cov(10:end,1:6) = estState.state_cov(7:end,1:6);
%     cov(10:end,10:end) = estState.state_cov(7:end,7:end);
% end