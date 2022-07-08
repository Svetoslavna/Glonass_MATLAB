function [output1] = SmallStepMove(in,step)
    Nt = 251; %current day
    N4 = 5;%current 4 year
    %C1 = 44195; C2 = 45290;
    G = 398600.4418 * 1e9;   %Earth`s gravitation field const, m^3/sec^2
    JulDate = 2451545; % Джулианская дата в течение 12 часов, 1 января 2000 года (UTC (SU))
    N_efemerids = 36525;% число эфемерид в день
    ti = 12300;
    R = sqrt(in.xyz(1)*in.xyz(1)+in.xyz(2)*in.xyz(2)+in.xyz(3)*in.xyz(3));
    j20 = 1082625.75 * 1e-9;        %koefficient pri vtoroi zonalnoi garmoniki 
    %j20 =  0.002420834274481/R;
    ae = 6378136.0;    % ekvatorialnii radius zemli, m
    w = 7.2921151467 * 1e-5;   %uglovaya skorost vrasheniya zemli
    ae2 = ae^2;
    J20_G_ae2 = 1.5 * j20 * G * ae2; %3/2*J20*G*ae^2
    r  = R; r2 = r^2; r3 = r * r2; r5 = r2 * r3;
    mur3 = G / r3;
    jmar = J20_G_ae2 / r5;  %3/2*J20*G*ae^2/(r * r^4)
    zr2 = (5.0 * in.xyz(3)^2) / r2;     %5z^2/r^2
    JDO = 1461*(N4-1) + Nt + 2450082.5 - (Nt - 3)/25;
    Tdelta = (JDO - JulDate)/N_efemerids;
    ERA = 2*pi*(0.779057273264+1.00273781191135448*(JDO-JulDate));
    GMST = ERA + 0.703270726e-7 + 0.0223603658710194*Tdelta + 0.67465784654e-5*Tdelta^2 - 0.21332e-11*Tdelta^3 - 0.1452308e-9*Tdelta^4 - 0.1784e-12*Tdelta^5;
    MT_UTC = 10800; %разница между MT и UTC(SU), sek
%     JDO = 2456177.5; GMST = 29191.44283;JDN = JDO + 0.5; Tdelta = (JDO - JulDate)/N_efemerids;
%     ERA = 2*pi*(0.779057273264+1.00273781191135448*(JDO-JulDate)); 
    S = GMST + w*(ti - MT_UTC);

        %уравнения дифференциала движения
    out.x = step * in.v_xyz(1);
    out.y = step * in.v_xyz(2);
    out.z = step * in.v_xyz(3);
 
    out.vx = step * (-mur3 * in.xyz(1) - jmar * in.xyz(1) * (1 - zr2))+in.a_xyz(1)*cos(S) - in.a_xyz(2)*sin(S);
    out.vy = step * (-mur3 * in.xyz(2) - jmar * in.xyz(2) * (1 - zr2))+in.a_xyz(1)*sin(S) + in.a_xyz(2)*cos(S);
    out.vz = step * (-mur3 * in.xyz(3) - jmar * in.xyz(3) * (3 - zr2))+in.a_xyz(3);

    
    out.ax = 0;
    out.ay = 0;
    out.az = 0;

    output1 = out;
end

