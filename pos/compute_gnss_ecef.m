function  log = compute_gnss_ecef(p,eph,obs)
% This function is to implement GNSS positioning with
% standard mode (without Iono, Trop, Es correction)
% or PPP mode.
% Output: log is a data struct that store the results.
%----------------------------------%
N = length(obs.tr_sow); % The number of positioning points
log = initOutputLog(p, obs);
satlog = initConstellationLog(p,log);

for i = 1:p.inval:N
    % Find the groud truth position
    if ~isnan(p.Grdpos.t(1))
        index = abs(p.Grdpos.t - obs.tr_sow(i)) < 0.01;
        grd.pos = p.Grdpos.pos(:,index);
        if isempty(grd.pos)
            continue;
        end
        if isfield(p.Grdpos, 'vel')
            grd.vel = p.Grdpos.vel(:,index);
        end
    else
        grdpos = p.Grdpos.pos;
    end

    p.i = i; % To debug
    cpt = obtainSatInfoStruct(p,satlog,eph,obs,i,log); 
    if max(cpt.num_sv(cpt.num_sv>0)) < p.min_sv
        continue;
    end
    if isempty(log.epoch_t) || (seconds(obs.datetime(i) - log.epoch_t(end)) > 1.5 && p.post_mode ~= p.mode_sps)
        % Rotate the sat pos to common reference frame
        [estState,~] = initialLsSolver(p,cpt);
        p.state0(1:3) = estState.pos;
        cpt = earth_rotation_corr(p,cpt);
        % Check elevation
        cpt = elevaz_check(p,cpt,estState.pos);
        % Open sky condition check
        cpt.is_open_sky = checkOpenSky(cpt.gps_range, cpt.gps_sat_pos, p.state0(1:3));
        tic
        if (p.post_mode == p.mode_dgnss || p.post_mode == p.mode_rtkfloat)...
                && ~isempty(p.eph_b) && ~isempty(p.obs_b)
            [cpt,n] = diff_corr_compute(p,cpt,obs.tr_posix(i));
            if max(cpt.num_sv(cpt.num_sv>0)) < p.min_sv
                continue;
            end 
            cpt = cpt_clear(cpt); % Clear the data where don't have diff correction
            cpt.corr_range = cpt.corr_range - cpt.diff_corr;
            cpt.phase_m = cpt.phase_m - cpt.diff_phase;
            cpt.doppler = cpt.doppler - cpt.diff_doppler;
            [estState,res] = weightLsSolver(p,cpt,true);
            if p.post_mode == p.mode_rtkfloat
                p.state0(1:3) = estState.pos;
                [estState,res] = weightLsSolver(p,cpt,false);
            end
        elseif p.post_mode == p.mode_ppp
            tdoy = doy(obs.tr_prime(1:3,i)); % Day of year
            [rt.week, rt.dow, rt.sow] = date2gnsst(obs.tr_prime(:,i)');
            cpt.IoFac = zeros(length(cpt.corr_range),1);
            cpt = correctBeidouCodeError(p,cpt);
            cpt = trop_iono_compute(p,eph,cpt,obs,p.state0(1:3),tdoy,rt,obs.tr_posix(i));
            if ~isempty(find(cpt.iono_delay~=0, 1))
                cpt.corr_range = cpt.corr_range - cpt.trop_delay - cpt.iono_delay;
                [estState,res] = weightLsSolver(p,cpt,true);
            end
        else
            [estState,res] = weightLsSolver(p,cpt,true);
        end
        % p.comp_t = toc;
        p.comp_t = NaN;
        p.comp_t_bcd = NaN;
        if ~isempty(estState.pos)
            % Save the initial state but not the first epoch.
            [p.state0, p.state_cov] = obtainInitEkfStateAndCov(p, estState);
            p.infor_ned = (p.state_cov)^(-1);
            p.augcost = NaN;
            p.augcost_bcd = NaN;
            p.pos_risk = NaN;
            p.pos_risk_bcd = NaN;
            p.num_meas_used = sum(cpt.num_sv)*2;
            if size(log.state_cov,1) == size(p.state_cov,1)
                log.epoch_t = [log.epoch_t, obs.datetime(i)];
                log = save_result(p,cpt,log,i,estState,res,grd,obs.datetime(i),true);
            end
        end
        dt = 0;
        continue;
    elseif p.post_mode ~= p.mode_sps
        dt = seconds(obs.datetime(i) - log.epoch_t(end));
        % EKF predict
        [p.state0, p.state_cov] = ekfPredict(p, p.state0, p.state_cov, dt);
    end
    if max(cpt.num_sv(cpt.num_sv>0)) < p.min_sv
        continue;
    end
    % Rotate the sat pos to common reference frame
    % [estState,~] = initialLsSolver(p,cpt);
    cpt = earth_rotation_corr(p,cpt);
    % Check elevation
    cpt = elevaz_check(p,cpt,p.state0(1:3));
    % Open sky condition check
    cpt.is_open_sky = checkOpenSky(cpt.gps_range, cpt.gps_sat_pos, p.state0(1:3));
    switch p.post_mode
        case p.mode_sps % Standard GNSS
            % if p.elev_mark ==0
            % log = save_result(p,cpt,log,i,re_pos,clock_bias,res);
            % else
            % Compute the final position
            % [re_pos,clock_bias,res] = userpos_Rcorr(p,cpt);
            tdoy = doy(obs.tr_prime(1:3,i)); % Day of year
            [rt.week, rt.dow, rt.sow] = date2gnsst(obs.tr_prime(:,i)');
            cpt.IoFac = zeros(length(cpt.corr_range),1);
            cpt = trop_iono_compute(p,eph,cpt,obs,p.state0(1:3),tdoy,rt,obs.tr_posix(i));
            cpt.corr_range = cpt.corr_range - cpt.trop_delay - cpt.iono_delay;
            [estState,res] = weightLsSolver(p,cpt);
            log = save_result(p,cpt,log,i,estState,res,grd,obs.datetime(i),false);
            % end
        case p.mode_ppp % PPP
            tdoy = doy(obs.tr_prime(1:3,i)); % Day of year
            [rt.week, rt.dow, rt.sow] = date2gnsst(obs.tr_prime(:,i)');
            %%-------------%%
            cpt.IoFac = zeros(length(cpt.corr_range),1);
            cpt = correctBeidouCodeError(p,cpt);
            cpt = trop_iono_compute(p,eph,cpt,obs,p.state0(1:3),tdoy,rt,obs.tr_posix(i));
            % Using the correction to the measurements
            if ~isempty(find(cpt.iono_delay~=0, 1))
                cpt.corr_range = cpt.corr_range - cpt.trop_delay - cpt.iono_delay;
                %%-------------%%
                [p,estState,res] = stateUpdate(p,cpt,dt);
                if ~isempty(estState.pos)
                    log = save_result(p,cpt,log,i,estState,res,grd,obs.datetime(i),false);
                end
            end
        case p.mode_dgnss % DGNSS
            if isempty(p.eph_b) || isempty(p.obs_b)
                warning('No differential source given');
                continue;
            end
            [cpt,n] = diff_corr_compute(p,cpt,obs.tr_posix(i));
            if isempty(n) || max(cpt.num_sv(cpt.num_sv>0)) < p.min_sv
                continue;
            end
            cpt.corr_range = cpt.corr_range - cpt.diff_corr;
            cpt.doppler = cpt.doppler - cpt.diff_doppler;
            [p,estState,res] = stateUpdate(p,cpt,dt);
            % [estState,res,~] = weightLsSolver(p,cpt,false);
            if ~isempty(estState.pos)
                log = save_result(p,cpt,log,i,estState,res,grd,obs.datetime(i),false);
            end
        otherwise
            warning('Unsupport positioning option');
    end
end






end