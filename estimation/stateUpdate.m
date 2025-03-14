function [p,estState,res] = stateUpdate(p, cpt, dt)

%-------------------%
% Initialize
estState.clock_sys = dictionary;
estState.clock_sys(p.gps.sys_num) = NaN;
estState.clock_sys(p.glo.sys_num) = NaN;
estState.clock_sys(p.gal.sys_num) = NaN;
estState.clock_sys(p.bds.sys_num) = NaN;
x_minus = p.state0;
num_user_states = p.modeToNumUserStates(p.state_mode);
[H_clk,x_clk] = formClkStatesAndH(cpt.num_sv);
if length(x_clk) + num_user_states + 1 ~= length(p.state0)
    error('current No. of sys does not match the previous epoch');
end
%------------------%
y_rho = cpt.corr_range;
p.num_sats_window = [p.num_sats_window(2:length(p.num_sats_window)), length(y_rho)];
y_dop = [];
if p.state_mode == p.pva_mode
    y_dop = cpt.doppler;
    H_clk = [H_clk; zeros(length(y_dop),size(H_clk, 2))];
end
num = length(y_rho) + length(y_dop); % The number of measurement
H = zeros(num,num_user_states+length(x_clk)+1);
H(:,num_user_states+1:num_user_states+length(x_clk)) = H_clk;
if p.state_mode == p.pva_mode
    H(length(y_rho)+1:end,end) = ones(length(y_dop),1);
    res_v = zeros(length(y_dop),1);
end
Range = zeros(length(y_rho),1);
r = zeros(length(y_rho),1);

s_pos_ecef = cpt.s_pos_ecef;
if p.post_mode == 1 && p.IGS_enable == 1
    s_pos_ecef = cpt.s_pos_prc;
end

s_v_ecef = cpt.sat_v_Rcorr;
for j=1:length(y_rho)
    Range(j)=norm(s_pos_ecef(:,j)-x_minus(1:3));
    los = (x_minus(1:3)-s_pos_ecef(:,j))'/Range(j)+...
        [-s_pos_ecef(2,j)*p.omge/p.c s_pos_ecef(1,j)*p.omge/p.c 0];
    H(j,1:3)=los;
    range_r = norm(cpt.sat_pos_Rcorr(:,j)-x_minus(1:3));
    los_r = (x_minus(1:3)-cpt.sat_pos_Rcorr(:,j))'/range_r;
    if p.state_mode == p.pva_mode
        H(j+length(y_rho),4:6) = los_r;
        res_v(j) = y_dop(j) - los_r*(x_minus(4:6) - s_v_ecef(:,j));
    end
    r(j) = Range(j)+sagnac(p,s_pos_ecef(:,j),x_minus(1:3));
end
H_os = H;
[R, ~] = constructMeasNoise(p, cpt, dt); %cpt.elev, cpt.svprn_mark
% measurement residual
res = y_rho - r - H_os(1:length(y_rho),4:end)*x_minus(4:end);
[x_minus, p.state_cov, flag] = checkClockReset(p, x_minus, p.state_cov, ...
    num_user_states, res, cpt); 
if flag == true
    res = y_rho - r - H_os(1:length(y_rho),4:end)*x_minus(4:end);
end
if p.state_mode == p.pva_mode
    res_v = res_v - x_minus(end);
    Rdop = p.sigmaSquare_dop*eye(length(y_dop));
    R = [R,zeros(length(y_rho),length(y_dop));
        zeros(length(y_dop),length(y_rho)),Rdop];
    res_all=[res;res_v];
else
    res_all=res;
end

% y - f(x0) = H (x - x0);
zk = res_all + H_os * x_minus;
switch p.est_mode
    case p.ekf_est
        tic;
        [x_plus, cov_plus] = ekfUpdate(x_minus, p.state_cov, res_all, H_os, R);
        comp_t = toc;
    case p.map_est
        lla = ecef2lla(x_minus(1:3)', 'WGS84');
        R_e2g=computeRotForEcefToNed(lla');
        R_pva = [R_e2g, zeros(3,6);
            zeros(3,3), R_e2g, zeros(3,3);
            zeros(3,6), R_e2g];
        Rot_e2g = [R_pva, zeros(9,length(x_minus)-9);
            zeros(length(x_minus)-9, 9), eye(length(x_minus)-9)];
        tic;
        H_os = H_os * Rot_e2g';
        cov_prior = Rot_e2g' * p.state_cov * Rot_e2g;
        [x_plus,cov_ned,p.infor_ned,p.augcost] = ...
            mapUpdate(ones(num,1),x_minus,cov_prior,res_all,H_os,R,Rot_e2g);
        cov_plus = Rot_e2g' * cov_ned * Rot_e2g;
        comp_t = toc;
        p.num_meas_used = num;
    case p.td_est
        [flag_rapid,p.num_sats_window] = checkRapidNumSatChange(p.num_sats_window, sum(cpt.num_sv~=0));
        if p.state_mode == p.pva_mode && flag_rapid == true
            p.state_cov = zeros(size(p.state_cov));
            p.state_cov(1:end-1,1:end-1) = diag(150^2*ones(1,length(x_minus)-1));
            p.state_cov(end,end) = 20^2;
        end
        lla = ecef2lla(x_minus(1:3)', 'WGS84');
        R_e2g=computeRotForEcefToNed(lla');
        R_pva = [R_e2g, zeros(3,6);
            zeros(3,3), R_e2g, zeros(3,3);
            zeros(3,6), R_e2g];
        Rot_e2g = [R_pva, zeros(9,length(x_minus)-9); 
            zeros(length(x_minus)-9, 9), eye(length(x_minus)-9)];

        tic;
        b = thresholdTest(p.td_lambda,p.state_cov, res_all, H_os, R);
        H_os = H_os * Rot_e2g';
        cov_prior = Rot_e2g' * p.state_cov * Rot_e2g;
        [x_plus,cov_ned,p.infor_ned,p.augcost] = ...
            mapUpdate(b, x_minus, cov_prior, res_all, H_os, R, Rot_e2g);
        cov_plus = Rot_e2g' * cov_ned * Rot_e2g;
        comp_t = toc;
        p.num_meas_used = sum(b);
    case p.raps_ned_est
        % Solve in NED frame
        % tic
        lla_deg = ecef2lla(x_minus(1:3)', 'WGS84');
        R_eg=computeRotForEcefToNed(lla_deg);
        if p.state_mode == p.pva_mode
            R_pva = [R_eg, zeros(3,6);
                zeros(3,3), R_eg, zeros(3,3);
                zeros(3,6), R_eg];
            [flag_rapid,p.num_sats_window] = checkRapidNumSatChange(p.num_sats_window, sum(cpt.num_sv~=0));
            % Rapid num of sat detected, may entering an open sky area,
            % reset prior covariance
            if flag_rapid == true
                p.state_cov = zeros(size(p.state_cov));
                p.state_cov(1:end-1,1:end-1) = diag(150^2*ones(1,length(x_minus)-1));
                p.state_cov(end,end) = 20^2;
            end
        else
            R_pva = R_eg;
        end
        Rot_e2g = [R_pva, zeros(num_user_states,length(x_minus)-num_user_states);
            zeros(length(x_minus)-num_user_states, num_user_states), eye(length(x_minus)-num_user_states)];
        Ht = H_os * Rot_e2g';
        xt_minus = Rot_e2g*(x_minus - x_minus);
        Pt_minus = Rot_e2g*p.state_cov*Rot_e2g';
        if p.state_mode == p.pva_mode
            num_constrain = 6;
            cov_spec_ecef = diag([p.raps.poshor_cov_spec; ...
                p.raps.poshor_cov_spec; p.raps.posver_cov_spec;...
                p.raps.velhor_cov_spec; p.raps.velhor_cov_spec;...
                p.raps.velver_cov_spec]);
            p_clk = diag([p.raps.va_cov_spec*ones(3,1);...
                p.raps.clk_cov_spec*ones(length(x_clk),1);...
                p.raps.dclk_cov_spec]);
            p_u = [cov_spec_ecef, zeros(6, length(x_minus)-6);
                zeros(length(x_minus)-6, 6), p_clk];
        elseif p.state_mode == p.pos_mode
            num_constrain = 3;
            cov_spec_ecef = diag([p.raps.poshor_cov_spec; ...
                p.raps.poshor_cov_spec; p.raps.posver_cov_spec]);
            p_clk = diag([p.raps.clk_cov_spec*ones(length(x_clk),1);...
                p.raps.dclk_cov_spec]);
            p_u = [cov_spec_ecef, zeros(3, length(x_minus)-3);
                zeros(length(x_minus)-3, 3), p_clk];
        end
        J_l = p_u^(-1);
        if p.raps_mode == p.nb_diag
            [flag,x_ned,cov_ned,b,J_out,p.augcost,num_iter,constraint,p.pos_risk,p.raps_penalty,comp_t] = ...
                mapRiskAverseNonBiSlackMaxJ(num_constrain,res_all,Ht,Pt_minus,R,...
                diag(diag(J_l)),xt_minus);
        elseif p.raps_mode == p.bi_diag
            [flag,x_ned,cov_ned,b,J_out,p.augcost,num_iter,constraint,p.pos_risk,p.raps_penalty,comp_t] = ...
                mapRiskAverseSlack(num_constrain,res_all,Ht,Pt_minus,R,...
                diag(diag(J_l)),xt_minus);
        elseif p.raps_mode == p.bi_diag_cvx
            [flag,x_ned,cov_ned,b,J_out,p.augcost,num_iter,constraint,p.pos_risk,p.raps_penalty,comp_t] = ...
                mapRiskAverseCvx(num_constrain,res_all,Ht,Pt_minus,R,...
                diag(diag(J_l)),xt_minus);
        end
        cov_plus = Rot_e2g' * cov_ned * Rot_e2g;
        p.num_meas_used = sum(b>0.001);
        b(b>0.01) = 1;
        b(b<=0.01) = 0;
        if p.state_mode == p.pva_mode
            p.raps_num_sat = sum(b(1:length(y_rho)));
        else
            p.raps_num_sat = sum(b);
        end
        p.infor_ned = J_out;
        p.raps_num_iter = num_iter;
        p.constraint = constraint;
        p.raps_flag = flag;
        x_plus = Rot_e2g'*x_ned + x_minus;
    otherwise
        error('Incorrect state estimation mode configuration');
end

p.comp_t = comp_t;

p.state0 = x_plus;
p.state_cov = cov_plus;

HH = [H_os(1:length(y_rho),1:3),H_clk(1:length(y_rho),:)];
if p.est_mode ~= p.ekf_est && p.est_mode ~= p.map_est
    b_rho = b(1:length(y_rho));
    HH = diag(b_rho)*HH;
end
hSqrtInv = (HH'*HH)^(-1);
p.GDOP = sqrt(trace(hSqrtInv));

estState.pos = x_plus(1:3);
if p.state_mode == p.pva_mode
    estState.vel = x_plus(4:6);
end
estState.clock_bias = x_plus(num_user_states+1);
estState.clock_drift = x_plus(end);

clk_est = x_plus(num_user_states+1:end-1);
j = 1;
for i = 1:length(cpt.num_sv)
    if cpt.num_sv(i) == 0
        continue;
    end
    estState.clock_sys(i) = clk_est(j);
    j=j+1;
end

end


