function [output] = ko_ti(out,step)
    Nt = 251; %current day
    N4 = 5;%current 4 year
    JDO = 1461*(N4-1) + Nt + 2450082.5 - (Nt - 3)/25;% + C1 + C2;
    JulDate = 2451545; % Джулианская дата в течение 12 часов, 1 января 2000 года (UTC (SU))
    N_efemerids = 36525;% число эфемерид в день
    JDN = JDO + 0.5; % число Джулиансикх дней для текущей даты
    Tdelta = (JDO - JulDate)/N_efemerids;
    ERA = 2*pi*(0.779057273264+1.00273781191135448*(JDO-JulDate));
    GMST = ERA + 0.703270726e-7 + 0.0223603658710194*Tdelta + 0.67465784654e-5*Tdelta^2 - 0.21332e-11*Tdelta^3 - 0.1452308e-9*Tdelta^4 - 0.1784e-12*Tdelta^5;
    MT_UTC = 10800;ti = 12300;w = 7.2921151467 * 1e-5;   %угловая скорость вращения Земли

    Sti = GMST + w*(ti - MT_UTC);
    one.xyz(1) = step * (out.xyz(1)*cos(Sti) + out.xyz(2)*sin(Sti));
    one.xyz(2) = step * (-out.xyz(1)*sin(Sti) + out.xyz(2)*cos(Sti));
    one.xyz(3) = step * (out.xyz(3));
    one.v_xyz(1) = step * (out.v_xyz(1)*cos(Sti) + out.v_xyz(2)*sin(Sti)+ w*one.xyz(2));
    one.v_xyz(2) = step * (-out.v_xyz(1)*sin(Sti) + out.v_xyz(2)*cos(Sti) - w*one.xyz(1));
    one.v_xyz(3) = step * (out.v_xyz(3));
    one.a_xyz(1) = 0;
    one.a_xyz(2) = 0;
    one.a_xyz(3) = 0;
    output = one;
end


