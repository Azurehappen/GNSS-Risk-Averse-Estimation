[pos_ecef, sow, Q, dtime] = readRtklibSol('E:\Github\Multi-GNSS-RAPS\data\univOfTexas\asterx4_rover.log');

% lever_arm = [0.5169;0.3668;0.0930];
lever_arm = [0;0;0];
[pos_gt, vel_ned, sow_gt, dtime_gt] = readTexasTxtGt('E:\Github\Multi-GNSS-RAPS\data\univOfTexas\ground_truth.log', lever_arm);

error_hor = NaN(1,length(sow));
error_ver = NaN(1,length(sow));
for i = 1:length(sow)
    if Q(i) ~= 1
        continue;
    end
    index = abs(sow_gt - sow(i)) < 0.01;
    pos_gt_i = pos_gt(:,index);
    if isempty(pos_gt_i)
        continue;
    end
    lla_gt_deg = ecef2lla(pos_gt_i', 'WGS84');
    wgs84 = wgs84Ellipsoid('meter');
    pos = pos_ecef(:,i);
    [xNorth,yEast,zDown] = ecef2ned(pos(1),pos(2),pos(3),lla_gt_deg(1),lla_gt_deg(2),lla_gt_deg(3),wgs84);
    error_hor(i) = norm(xNorth^2 + yEast^2);
    error_ver(i) = abs(zDown);
end

figure
plot(dtime, error_hor,'.')
grid on
title('Horizontal Error (RTKLIB to GT)')

% figure
% plot(dtime, error_ver,'.')
% grid on
% title('Horizontal Error (RTKLIB to GT)')