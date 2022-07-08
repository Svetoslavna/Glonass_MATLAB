%% Исходные данные
clc, clear, close all
GOD = 2012;% ближайший к N предшествубщий високосныйы год

Nt = 251; %current day
N4 = 5;%current 4 year, P1 = 01 - величина интервала между значениями tb 30 мин
% if Nt>=1 && Nt<= 366
%     J = 1;
% elseif Nt>= 367 && Nt<= 731
%     J = 2;
% elseif Nt>= 732 && Nt <=1096 
%     J = 3;
% elseif Nt>= 1097 && Nt <= 1461
%     J = 4;
% end
% YY = 1996+4*(N4 - 1)+(J-1);
G = 398600.4418 * 1e9;   %Earth`s gravitation field const, m^3/sec^2
%Earth`s gravitation field const for the moon and the sun
Gm = 4902.799e9; Gs = 1.3271244e20;
ti = 12300; tb = 11700;dt = ti - tb;step_size = 60; %tb - номер интервала 
if dt < 0 % эффект неравномерной шкалы времени: при переходе через начало суток возникает отрицательная разность, если последнее прохождение восходящего узла было в предыдущую дату, а Ti измеряется в секундах текущей даты.
    dt = dt + 86400;
end
steps = round(dt/step_size);
c = 299792458; %speed of light in vacuum m/s
j20 = 10827 * 1e-7;        %коэф. при 2й зональной гармонике j20n = -484.16495e-6;% нормализованное 
ae = 6378136.0;     % экваториальный радиус Земли, м
w = 7.2921151467 * 1e-5;   %угловая скорость вращения Земли
ae2 = ae^2;
w2 = w^2;
J20_G_ae2 = 1.5 * j20 * G * ae2; %3/2*J20*G*ae^2
as = 1.49598e11;    % осн.орбита земли вокруг солнца, м
es = 0.016719;      %эксцентриситет орбиты Земли вокруг Солнца
am = 3.84385243e8; % осн.орбита земли вокруг луны, м
em = 0.054900489; %эксцентриситет орбиты Земли вокруг луны
im = 0.0898041080; %среднее наклонение орбиты Луны к эклиптике

xyz = [7003.008789e3, -12206.626953e3, 21280.765625e3 ];
v_xyz = [0.7835417e3, 2.8042530e3, 1.3525150e3];
a_xyz = [0,1.7e-9, -5.41e-9];
test_point0 = struct ('xyz', xyz, 'v_xyz', v_xyz, 'a_xyz', a_xyz);
Rteor = sqrt(test_point0.xyz(1)*test_point0.xyz(1)+test_point0.xyz(2)*test_point0.xyz(2)+test_point0.xyz(3)*test_point0.xyz(3));

JD0teor = 2456177.5; GMSTteor = 29191.44283;
%J_moon_teor = [-5.035590e-7, 7.379024e-7, -1.64803e-6];
%J_sun_teor = [4.428827e-7, 3.541631e-7, -8.911601e-7];
C1 = 44195; C2 = 45290;
% _____________from App.S_3
MT_UTC = 10800; %разница между MT и UTC(SU), sek = 3 часа
JulDate = 2451545;%Джулианская дата в течение 12 часов, 1 января 2000 года (UTC (SU)) - в методичке - 2451545
N_efemerids = 36525;% число эфемерид в день

JDO = round(1461*(N4-1) + Nt + 2450092.5 - (Nt - 3)/25);%+ C2 + C1;2450082.5
JDN = JDO - 0.5;

        %%для проверки
        %JulDate = juliandate(datetime(("2012-01-01 12:00:00")));%дата, часы, минуты, секунды
        % JDN = juliandate(datetime(("2012-09-07 12:00:00")));
        % JDO = JDN - 0.5;
            % %проверка високосности года и определение ближайшего прошедшего висок.
            % %года.
            % for i = 1
            %     visok = leapyear(GOD);
            %     if visok == 1
            %         K = 1;
            %     else
            %         i = 0;GODo = GOD;
            %         while i > 1 || i < 1
            %             i = leapyear(GODo);
            %             GODo = GODo - 1; 
            %         end
            %         GODo = GODo + 1;
            %         if GOD == GODo + 1
            %             K = 1.75;
            %         elseif GOD == GODo + 2 
            %             K = 1.5;
            %         elseif GOD == GODo + 3
            %             K = 1.25;
            %         end
            %     end
            % end
            % JD = (4712 + GOD)*365.25 + GOD/400 - GOD/100 + K;

    Tdelta = (JDN - JulDate)/N_efemerids;
    tu = JDN - JulDate;
    ERA = 2*pi*(0.779057273264+1.00273781191135448*tu);%Earth rotation angle,замена tu <--(JDO-JulDate)
            
    
    %считают верно, но нужен перевод из градусов в радианы:
    
        %GMSTT = 6.697374558 + 0.06570982441908 - JulDate + 1.00273790935 +  0.000026 *((JD - JulDate)/36525)^2;
        %     GMST11 = ERA + 0.0223603658710194*Tdelta + 6.697374558*Tdelta^2 + 0.06570982441908*Tdelta^ + 1.00273790935*Tdelta^4 +  0.000026 *((JD - JulDate)/36525)^2*Tdelta^5
%             gst0=1.753368560+628.3319706889*Tdelta+6.7707e-6*Tdelta.^2-4.5e-10*Tdelta.^3
%             we=7.29211585530e-5;
%             gst=rad2deg(wrapTo2Pi(gst0+we*seconds(timeofday(datetime(JDO,'ConvertFrom','juliandate','TimeZone','UTC')))))
%              GMST0 = JD2GMST(JDO)
%              GMST = deltalongeci2ecef(JDO)
    GMST = ERA + 0.703270726e-7 + 0.0223603658710194*Tdelta + 0.67465784654e-5*Tdelta^2 - 0.21332e-11*Tdelta^3 - 0.1452308e-9*Tdelta^4 - 0.1784e-12*Tdelta^5
          %JDO = 2455927.5; GMST = 29191.44283;JDN = JDO + 0.5; Tdelta = (JDO - JulDate)/N_efemerids;
          %ERA = 2*pi*(0.779057273264+1.00273781191135448*(JDO-JulDate));

%промеждуточные коэффициенты
a =fix( JDN + 32044); b = fix((4*a+3)/146097); c1 = fix(a - (146097*b)/4); d = fix((4*c1 + 3)/1461); e = fix(c1 - (1461*d)/4); m = fix((5*e+2)/153);
%Расч. Грегорианского дня, месяца, года
day = e -fix((153*m+2)/5)+1; month = m+3-12*fix((m/10)); year = 100*b+d-4800+fix((m/10));
Stb = GMST + w*(tb - MT_UTC);
    test_point.xyz(1) = test_point0.xyz(1)*cos(Stb) - test_point0.xyz(2)*sin(Stb);
    test_point.xyz(2) = test_point0.xyz(1)*sin(Stb) + test_point0.xyz(2)*cos(Stb);
    test_point.xyz(3) = test_point0.xyz(3);
    test_point.v_xyz(1) = test_point0.v_xyz(1)*cos(Stb) - test_point0.v_xyz(2)*sin(Stb)- w*test_point.xyz(2);
    test_point.v_xyz(2) = test_point0.v_xyz(1)*sin(Stb) + test_point0.v_xyz(2)*cos(Stb) + w*test_point.xyz(1);
    test_point.v_xyz(3) = test_point0.v_xyz(3);
    test_point.a_xyz(1) = test_point0.a_xyz(1);
    test_point.a_xyz(2) = test_point0.a_xyz(2);
    test_point.a_xyz(3) = test_point0.a_xyz(3);
%    

%Расч. коорд. кажущегося солнца для момента (ti)
TkoordTi = (JDO + (ti - MT_UTC)/86400 - JulDate)/N_efemerids;
qs_ti = 6.2400601269 + 628.3019551714*TkoordTi - 0.0000026820*TkoordTi^2; %средняя аномалия Солнца, рад
ws_ti = -7.6281824375 + 0.0300101976*TkoordTi + 0.0000079741*TkoordTi^2; % средняя тропическая долгота солнечной орбиты перигея, рад
eps_ti = 0.4090926006 - 0.0002270711*TkoordTi; % средний наклон земного экватора к эклиптике, рад


% eps_ti = 0.46;
% %средняя тропическая долгота солнечной орбиты перигея, рад
% ws_ti = 4.88 + 628.33*TkoordTi + 0.0053e-3*TkoordTi^2;
% TT = 2*pi*sqrt(as^3/Gs);n = 2*pi/TT;%средняя угловая скорость %n = sqrt(Gs/as^3)
% h=sqrt(n*as^2*(1-es^2));
% M =  n*(ti - TkoordTi);%средняя аномалия
for i = 0
    Es_ti = qs_ti;
    while abs((qs_ti - (Es_ti - es*sin(Es_ti)))/(1-es*cos(Es_ti)))>1e-8 && i<100
        
        Es_ti = Es_ti + ((qs_ti - (Es_ti - es*sin(Es_ti)))/(1-es*cos(Es_ti)));
        i = i + 1;  
    end
    Es_ti = Es_ti + ((qs_ti - (Es_ti - es*sin(Es_ti)))/(1-es*cos(Es_ti)));
end

    sinOmegs1 = ((sqrt((1-es^2))/(1-es*cos(Es_ti))));%sqrt((1-es^2)*sin(Es_ti))/(1-es*cos(Es_ti))
    cosOmegs1 = (cos(Es_ti)*es)/sqrt(1-es^2*cos(Es_ti)^3);%(cos(Es_ti)-es)/(1-es*cos(Es_ti));
    ksis_ti = cosOmegs1*cos(ws_ti) - sinOmegs1*sin(ws_ti);
    etas_ti = (sinOmegs1*cos(ws_ti)+ cosOmegs1*sin(ws_ti))*cos(eps_ti);
    zetas_ti = (sinOmegs1*cos(ws_ti)+cosOmegs1*sin(ws_ti))*sin(eps_ti);
    rs_ti = as*(1-es*cos(Es_ti));%расстояние до Солнца в ti

    % Расч. коорд. кажущегося солнца для момента (ti-delta_t)
    delta_t = rs_ti/c; %время прохождения света от солнца на Землю
    T = (JDO + (ti-delta_t - MT_UTC)/86400 - JulDate)/N_efemerids;
    qs = 6.2400601269 + 628.3019551714*T - 0.0000026820*T^2; %средняя аномалия Солнца, рад
    ws = -7.6281824375 + 0.0300101976*T + 0.0000079741*T^2; % средняя тропическая долгота солнечной орбиты перигея, рад
    eps = 0.4090926006 - 0.2270711e-3*T; % средний наклон земного экватора к эклиптике, рад
for i = 0
    Es = qs;
    while abs((qs - (Es - es*sin(Es)))/(1-es*cos(Es)))>1e-8 && i<100
        Es = Es + ((qs - (Es - es*sin(Es)))/(1-es*cos(Es)));
        i = i + 1;
    end
    Es = Es + ((qs - (Es - es*sin(Es)))/(1-es*cos(Es)));
end
                %из методички     sinOmegs = ((sqrt(1-es^2))*sin(Es))/(1-es*cos(Es));
                %     cosOmegs = (cos(Es)-es)/(1-es*cos(Es));
    sinOmegs = ((sqrt((1-es^2))/(1-es*cos(Es))));
    cosOmegs = (cos(Es)*es)/sqrt(1-es^2*cos(Es)^3);
    ksis = cosOmegs*cos(ws) - sinOmegs*sin(ws);
    etas = (sinOmegs*cos(ws)+ cosOmegs*sin(ws))*cos(eps);
    zetas = (sinOmegs*cos(ws)+cosOmegs*sin(ws))*sin(eps);
    rs = as*(1-es*cos(Es));%расстояние до Солнца в (ti-delta_t)

    %Координаты истинного Солнца в любом произвольном мгновении ti 
    Xso = rs*ksis; Yso = rs*etas; Zso = rs*zetas; 

    %Пересчет координат очевидного и истинного Солнца
    S = GMST + w*(ti - MT_UTC);
    Xs = Xso*cos(S)+Yso*sin(S); Ys = -Xso*sin(S) + Yso*cos(S); Zs = Zso;
 
    Tm = (JDO + (tb - MT_UTC)/86400 - JulDate)/N_efemerids;
    qm = 2.3555557435 + 8328.691425719*Tm + 0.1545547e-3*Tm^2;%средняя аномалия луны, рад
    wm = -2.1824391966 - 33.7570459536*Tm + 0.362262e-4*Tm^2; % средняя долгота восходящего узла Луны, рад
    epsm = 0.4090926006 - 0.2270711e-3*Tm; % средний наклон земного экватора к эклиптике, рад
    g = 1.4547885346 + 71.0176852437*Tm - 0.1801481e-3*Tm^2;% Г' - средняя долгота перигея орбиты Луны, рад 

    
%     %средняя долгота восход. узла орбиты Луны, рад
%     wm = 4.52 - 33.76*Tm + 0.0036e-2*Tm^2 + 0.0039e-5*T^3;
%     %средняя долгота перигея Солнца
%     g = 4.91 + 0.03*Tm + 0.0079e-3*T^2 + 0.0058e-5*Tm^3;
%     TT = 2*pi*sqrt(am^3/Gm);n = 2*pi/TT;%средняя угловая скорость %n = sqrt(Gs/as^3)
%     h=sqrt(n*am^2*(1-em^2));qm =  n*(ti - Tm);%средняя аномалия
%     %TT = 18.5996;%лет, период лунной процессии
%     epsm = 0.027; %g = 1.45;

for i = 0
    Em = qm;
    while abs((qm - (Em - em*sin(Em)))/(1-em*cos(Em)))>1e-8 && i<100
        Em = Em + ((qm - (Em - em*sin(Em)))/(1-em*cos(Em)));
        i = i + 1;
    end
    Em = Em + ((qm - (Em - em*sin(Em)))/(1-em*cos(Em)));
end
            %из методички Sin_m = ((sqrt(1-em^2))*sin(Em))/(1-em*cos(Em)); Cos_m = (cos(Em) - em)/(1 - em*cos(Em));
% https://scask.ru/r_book_mor.php?id=53   стр. 99, (4.46):       
Sin_m = ((sqrt((1-em^2))/(1-em*cos(Em))));%sqrt((1-es^2)*sin(Es_ti))/(1-es*cos(Es_ti)) 
Cos_m = (cos(Em)*em)/sqrt(1-es^2*cos(Em)^3);%(cos(Es_ti)-es)/(1-es*cos(Es_ti));
    
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ksi*, eta*, zeta*
    ksim_zv = 1 - cos(wm)^2*(1 - cos(im)); etam_zv = sin(wm)*sin(im); zetam_zv = cos(wm)*sin(im);

    ksim11 = sin(wm)*cos(wm)*(1-cos(im)); ksim12 = 1 - sin(wm)^2*(1 - cos(im));
    etam11 = ksim_zv*cos(epsm) - zetam_zv*sin(epsm); etam12 = ksim11*cos(epsm) + etam_zv*sin(epsm);
    zetam11 = ksim_zv*sin(epsm) + zetam_zv*cos(epsm); zetam12 = ksim11*sin(epsm) - etam_zv*cos(epsm);

    ksim = (Sin_m*cos(g) + Cos_m*sin(g))*ksim11 + (Cos_m*cos(g) - Sin_m*sin(g))*ksim12;
    etam = (Sin_m*cos(g) + Cos_m*sin(g))*etam11 + (Cos_m*cos(g) - Sin_m*sin(g))*etam12;
    zetam = (Sin_m*cos(g) + Cos_m*sin(g))*zetam11 + (Cos_m*cos(g) - Sin_m*sin(g))*zetam12;
    rm = am*(1 - em*cos(Em));
  
%%%%%Расчет ускорений, вызванных гравитационными возмущениями Солнца и Луны
    Gsn = Gs/rs^2; Xsn = Xs/rs; Ysn = Ys/rs; Zsn = Zs/rs;
    delta_k = ((ksis-Xsn)^2 + (etas-Ysn)^2 + (zetas-Zsn)^2)^(3/2);
    
    Gmn = Gm/rm^2; Xmn = test_point.xyz(1)/rm; Ymn = test_point.xyz(2)/rm; Zmn = test_point.xyz(3)/rm;
    delta_km = ((ksim-Xmn)^2 + (etam-Ymn)^2 + (zetam-Zmn)^2)^(3/2);
    
    Joxs = Gsn*((ksis - Xsn)/delta_k - ksis);Joys = Gsn*((etas - Ysn)/delta_k - etas);Jozs = Gsn*((zetas - Zsn)/delta_k - zetas);
    Joxm = Gmn*((ksim - Xmn)/delta_km - ksim);Joym = Gmn*((etam - Ymn)/delta_km - etam);Jozm = Gmn*((zetam - Zmn)/delta_km - zetam);
  
for step = 0:1:(steps-1)
    test_point = struct(StepMove(test_point, step_size));
%     fprintf("\nStep = %d\n", step);
%     test_point.xyz
end
% 
Sti = GMST + w*(ti - MT_UTC);
    test_point2.xyz(1) = test_point.xyz(1)*cos(Sti) + test_point.xyz(2)*sin(Sti);
    test_point2.xyz(2) = -test_point.xyz(1)*sin(Sti) + test_point.xyz(2)*cos(Sti);
    test_point2.xyz(3) = test_point.xyz(3);
    test_point2.v_xyz(1) = test_point.v_xyz(1)*cos(Sti) + test_point.v_xyz(2)*sin(Sti)+ w*test_point2.xyz(2);
    test_point2.v_xyz(2) = -test_point.v_xyz(1)*sin(Sti) + test_point.v_xyz(2)*cos(Sti) - w*test_point2.xyz(1);
    test_point2.v_xyz(3) = test_point.v_xyz(3);
    test_point2.a_xyz(1) = test_point.a_xyz(1);
    test_point2.a_xyz(2) = test_point.a_xyz(2);
    test_point2.a_xyz(3) = test_point.a_xyz(3);
    
 %данные из методички
%  Joxm = -5.035590e-7; Joym=7.379024e-7;Jozm = -1.64803e-6;
% Joxs = 4.428827e-7;Joys = 3.541631e-7;Jozs = -8.911601e-7;

    test_point3.xyz(1) = test_point2.xyz(1)+(Joxs+Joxm);%*dt^2/2;
    test_point3.xyz(2) = test_point2.xyz(2)+(Joys+Joym);%*dt^2/2;
    test_point3.xyz(3) = test_point2.xyz(3)+(Jozs+Jozm);%*dt^2/2;
    test_point3.v_xyz(1) = test_point2.v_xyz(1)+(Joxs+Joxm);%*dt;
    test_point3.v_xyz(2) = test_point2.v_xyz(2)+(Joys+Joym);%*dt;
    test_point3.v_xyz(3) = test_point2.v_xyz(3)+(Jozs+Jozm);%*dt;
    test_point3.a_xyz(1) = test_point2.a_xyz(1);
    test_point3.a_xyz(2) = test_point2.a_xyz(2);
    test_point3.a_xyz(3) = test_point2.a_xyz(3);
    
% %     
%         test_point3.xyz(1) = test_point2.xyz(1)+(Joxs+Joxm)*dt^2/2;
%     test_point3.xyz(2) = test_point2.xyz(2)+(Joys+Joym)*dt^2/2;
%     test_point3.xyz(3) = test_point2.xyz(3)+(Jozs+Jozm)*dt^2/2;
%     test_point3.v_xyz(1) = test_point2.v_xyz(1)+(Joxs+Joxm)*dt;
%     test_point3.v_xyz(2) = test_point2.v_xyz(2)+(Joys+Joym)*dt;
%     test_point3.v_xyz(3) = test_point2.v_xyz(3)+(Jozs+Jozm)*dt;
%     test_point3.a_xyz(1) = 0;
%     test_point3.a_xyz(2) = 0;
%     test_point3.a_xyz(3) = 0;

respoint = struct ('x',{7523.174819e3},'y',{-10506.961965e3}, 'z',{21999.239413e3}, ...
                    'vx',{0.950126007e3}, 'vy',{2.855687825e3},'vz',{1.040679862e3},'ax',{0},'ay',{0},'az',{0});resPoint = struct2array(respoint);
ResPoint.xyz = resPoint(1:3);ResPoint.v_xyz = resPoint(4:6);ResPoint.a_xyz = resPoint(7:9);

error = mmin(ResPoint,test_point3);
Error = sqrt(error.xyz(1)*error.xyz(1)+error.xyz(2)*error.xyz(2)+error.xyz(3)*error.xyz(3));

fprintf('\nError r = %d',Error);disp(' meters');fprintf("\nerror x = %d\n",error.xyz(1));fprintf("\nerror y = %d\n",error.xyz(2));fprintf("\nerror z = %d\n",error.xyz(3));fprintf("\nerror vx = %d\n",error.v_xyz(1));fprintf("\nerror vy = %d\n",error.v_xyz(2));fprintf("\nerror vz = %d\n",error.v_xyz(3));fprintf("\nerror ax = %d\n",error.a_xyz(1));fprintf("\nerror ay = %d\n",error.a_xyz(2));fprintf("\nerror az = %d\n",error.a_xyz(3));
