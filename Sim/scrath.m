%         if ctrl.loopselect == 0
%             figure('Name', ['Flight Path Angle and Heading' num2str(gyroNoiseDensity(j))]);
%             subplot(2,1,1); plot(tvec(2:end-1),180/pi * cmd_g(1,2:end),tvec,flightPathAngle * 180/pi);
%             ylabel('Flight Path Angle (deg)')
%             legend('Guidance Command','Response')
%             subplot(2,1,2); plot(tvec(2:end-1),180/pi * cmd_g(2,2:end),tvec,heading * 180/pi)
%             ylabel('headinging Angle (deg)'); xlabel('Time(sec)'); grid on
%
%             figure('Name', ['Altitude and Velocity' num2str(gyroNoiseDensity(j))])
%             subplot(2,1,1); plot(tvec,V)
%             title(['Flight Conditions with Gyro noise of ' num2str(gyroNoiseDensity(j))]); ylabel('Velocity (m/s)')
%             subplot(2,1,2); plot(tvec,z-constants.RADIUS_EQ)
%             ylabel('Altitude (m)');xlabel('Time(sec)'); grid on
%         end
%
%         if ctrl.loopselect <= 1
%             figure('Name', ['Angle Outputs' num2str(gyroNoiseDensity(j))]);
%             subplot(3,1,1); plot(tvec,180/pi * cmd_a(1,1:length(aoa)),tvec,bank * 180/pi)
%             title(['Angle Output with Gyro noise of ' num2str(gyroNoiseDensity(j))]); ylabel('bank (deg)')
%             legend('Angle Command','Response')
%             subplot(3,1,2); plot(tvec,180/pi * cmd_a(2,1:length(aoa)),tvec,aoa * 180/pi);
%             ylabel('aoa (deg)')
%             subplot(3,1,3); plot(tvec,180/pi * cmd_a(3,1:length(aoa)),tvec,ss * 180/pi)
%             ylabel('ss (deg)'); xlabel('Time(sec)'); grid on
%         end
%
%         figure('Name', ['Angular Rate Outputs' num2str(gyroNoiseDensity(j))]);
%         subplot(3,1,1); plot(tvec(1:end-1),cmd_r(1,:) * 180/pi,tvec,p * 180/pi)
%         title(['Angular Rate Output with Gyro noise of ' num2str(gyroNoiseDensity(j))]); ylabel('p (deg/sec)')
%         legend('Command','Response')
%         subplot(3,1,2); plot(tvec(1:end-1),cmd_r(2,:) * 180/pi,tvec,q * 180/pi);
%         ylabel('q (deg/sec)')
%         subplot(3,1,3); plot(tvec(1:end-1),cmd_r(3,:) * 180/pi,tvec,r * 180/pi)
%         ylabel('r (deg/sec)'); xlabel('Time(sec)'); grid on
%
%         % figure;
%         % subplot(3,1,1); plot(tvec,p * 180/pi)
%         % title('Angular Rate Output'); ylabel('p (deg/sec)')
%         % subplot(3,1,2); plot(tvec,q * 180/pi);
%         % ylabel('q (deg/sec)')
%         % subplot(3,1,3); plot(tvec,r * 180/pi)
%         % ylabel('r (deg/sec)'); xlabel('Time(sec)'); grid on
%
%         % figure;
%         % subplot(3,1,1); plot(tvec,bank * 180/pi)
%         % title('Angle Displacements'); ylabel('Bank Angle (deg)')
%         % subplot(3,1,2); plot(tvec,aoa * 180/pi);
%         % ylabel('Angle of Attack (deg)')
%         % subplot(3,1,3); plot(tvec,ss * 180/pi)
%         % ylabel('Sideslip Angle (deg)'); xlabel('Time(sec)'); grid on
%
%         figure('Name', ['Angular Accel Output' num2str(gyroNoiseDensity(j))]);
%         subplot(3,1,1); plot(tvec(1:end-1),flin(1,:) * 180/pi,tvec(1:end-1),pDot(1:end-1) * 180/pi)
%         title(['Angular Accel Output with Gyro noise of ' num2str(gyroNoiseDensity(j))]); ylabel('pDot (deg/s^2)')
%         % ylim([-2 2])
%         legend('Accel Command','Response')
%         subplot(3,1,2); plot(tvec(1:end-1),flin(2,:) * 180/pi,tvec(1:end-1),qDot(1:end-1) * 180/pi);
%         ylabel('qDot (deg/s^2)')
%         % ylim([-2 2])
%         subplot(3,1,3); plot(tvec(1:end-1),flin(3,:) * 180/pi, tvec(1:end-1),rDot(1:end-1) * 180/pi)
%         ylabel('rDot (deg/s^2)'); xlabel('Time(sec)'); grid on
%         % ylim([-2 2])
%
%         % figure;
%         % subplot(3,1,1); plot(tvec(1:end-1),pDot(1:end-1) * 180/pi)
%         % title('Angular Accel Output'); ylabel('pDot (deg/s^2)')
%         % subplot(3,1,2); plot(tvec(1:end-1),qDot(1:end-1) * 180/pi);
%         % ylabel('qDot (deg/s^2)')
%         % subplot(3,1,3); plot(tvec(1:end-1),rDot(1:end-1) * 180/pi)
%         % ylabel('rDot (deg/s^2)'); xlabel('Time(sec)'); grid on
%
%         % reorder flaps to match report numbering
%         fdplot = [flapDeflection(2,:);flapDeflection(3,:);flapDeflection(4,:);flapDeflection(5,:);flapDeflection(6,:);flapDeflection(7,:);flapDeflection(8,:);flapDeflection(1,:)];
%         dataSave.finDeflectSave(:,:,j) = fdplot;
%
%         figure('Name', ['Flap Deflections ' num2str(gyroNoiseDensity(j))])
%         plot(tvec(1:end-1),fdplot(:,1:end-1))
%         title(['Flap Deflect Angles with Gyro noise of ' num2str(gyroNoiseDensity(j))]); ylabel('deg')
%         legend('F1','F2','F3','F4','F5','F6','F7','F8')
%         xlabel('Time (Sec)')
%         grid on

end % end of gyro noise loop


% plotGyroNoise(dataSave)

% figure
% plot(tvec,gammaDot,tvec,gammaDotFilt)
% figure
% plot(tvec,headingDot,tvec,headingDotFilt)
%
% figure
% plot(tvec,pDot,tvec,pDotFilt)
% figure
% plot(tvec,qDot,tvec,qDotFilt)
% figure
% plot(tvec,rDot,tvec,rDotFilt)

% figure
% subplot(3,1,1)
% plot(tvec(2:end),vlin(1,:),tvec(2:end),flin(1,:))
% legend('raw','filtered')
% subplot(3,1,2)
% plot(tvec(2:end),vlin(1,:),tvec(2:end),flin(1,:))
% subplot(3,1,3)
% plot(tvec(2:end),vlin(1,:),tvec(2:end),flin(1,:))
% V = zeros(npoints,1); Vdot = V;% vehicle velocity (m/s)
        % flightPathAngle = zeros(npoints,1); gammaDot = flightPathAngle; % flight path angle (rad)
        % heading = zeros(npoints,1); headingDot = heading; % headinging angle (rad)
        % chi = zeros(npoints,1); chidot = chi; % track angle (rad)
        % p = zeros(npoints,1); pDot = p; % angular velocity about x axis (rad/s)
        % q = zeros(npoints,1); qDot = q; % angular velocity about y axis (rad/s)
        % r = zeros(npoints,1); rDot = r; % angular velocity about z axis (rad/s)
        % pDotFilt = r; qDotFilt = r; rDotFilt = r;
        % gammaDotFilt = r; headingDotFilt = r;
        % rotBodyFromECEF = zeros(3,3,npoints);
        % rotNedFromBody  = zeros(3,3,npoints);
        % rotNedFromECEF  = zeros(3,3,npoints);

        % % Kinematics
        % z     = zeros(npoints,1); zdot = z; % orbital radius (m)
        % lat   = zeros(npoints,1); latdot = lat; % latitude (rad)
        % lon   = zeros(npoints,1); londot = lon; % longitude (rad)
        % bank  = zeros(npoints,1); bdot = bank; % bank angle (rad)
        % aoa   = zeros(npoints,1); adot = aoa; % angle of attack (rad)
        % ss    = zeros(npoints,1); ssdot = ss; % sideslip angle (rad)
        % pitch = zeros(npoints,1);
        % roll  = zeros(npoints,1);

        % % Telemetry
        % vlin = zeros(3,npoints-1);
        % flin = zeros(3,npoints-1);
        % if ~exist('cmd_a'); cmd_a = zeros(3,npoints); end
        % cmd_r = zeros(3,npoints-1); % [p, q, r]
        % cmd_g = zeros(2,npoints-1); % [flightPathAngle; heading]

        % %% Initial Conditions
        % M = 30; % Mach number (dless) and speed of sound (m/s)
        % V(1) = M * c;
        % heading(1) = 0;
        % p(1) = 0; q(1) = 0; r(1) = 0;
        % alt = 100000; %118000;
        % alt_init = alt; % save initial altitude for filename
        % z(1) = alt + constants.RADIUS_EQ;
        % lat(1) = deg2rad(37.3352); lon(1) = deg2rad(121.8811);
        % bank(1) = 0; aoa(1) = 0; ss(1) = 0;

        % flapDeflection = zeros(8,npoints); % flap deflection angles (deg)
        % fc = zeros(8,1);