function [PHiDE,ThetaDE] = AccDE01_7_21_2021(Dataace2,T,x01,x901,x_901,y01,y901,y_901,z01,z901,z_901,xi01,xi901,xi_901,yi01,yi901,yi_901,zi01,zi901,zi_901,calibT,data_flag)
%output: PHiDE is the spinning anlge of the tool
% 	 ThetaDE is the inclination of angle 
% Inputs: Dataace2  accelerometer input 6x1 matrix 
%	T CCB temperature sensor 
%	xi01,xi901,xi_901,yi01,yi901,yi_901,zi01,zi901,zi_901 calibration data of accelermeter 2 vs temperature  Matrix 18x1 of each variable
%	x01,x901,x_901,y01,y901,y_901,z01,z901,z_901 calibration data of accelermeter 1 vs temperature Matrix 18x1 of each variable
%       calibT 18X1 Temperature at each calibration happned
    FFx=Dataace2(1,4);
    FFy=Dataace2(1,5);
    FFz=Dataace2(1,6);
    FFx1=Dataace2(1,1);
    FFy1=Dataace2(1,2);
    FFz1=Dataace2(1,3);
    hhhd=1;
    jj=1;
    i=1;
    while (hhhd)&&(jj<length(calibT))
        jj=jj+1;
        if T(i)>calibT(jj)
            hhhd=1;
        else 
            hhhd=0;
        end 
    end
    %%
    %%
    %%Temperature calibration for each Accelerometer axis
    %%
    xx0=x01(jj);
    xx02=x01(jj-1);
    xx90=x901(jj);
    xx902=x901(jj-1);
    xx_90=x_901(jj);
    xx_902=x_901(jj-1);
    zz0=z01(jj);
    zz02=z01(jj-1);
    zz90=z901(jj);
    zz902=z901(jj-1);
    zz_90=z_901(jj);
    zz_902=z_901(jj-1);
    yy0=y01(jj);
    yy02=y01(jj-1);
    yy90=y901(jj);
    yy902=y901(jj-1);
    yy_90=y_901(jj);
    yy_902=y_901(jj-1);
    xxi0=xi01(jj);
    xxi02=xi01(jj-1);
    xxi90=xi901(jj);
    xxi902=xi901(jj-1);
    xxi_90=xi_901(jj);
    xxi_902=xi_901(jj-1);
    zzi0=zi01(jj);
    zzi02=zi01(jj-1);
    zzi90=zi901(jj);
    zzi902=zi901(jj-1);
    zzi_90=zi_901(jj);
    zzi_902=zi_901(jj-1);
    yyi0=yi01(jj);
    yyi02=yi01(jj-1);
    yyi90=yi901(jj);
    yyi902=yi901(jj-1);
    yyi_90=yi_901(jj);
    yyi_902=yi_901(jj-1);
    T2=calibT(jj);
    T1=calibT(jj-1);
    x0=(xx0-xx02)*(T(i)-T1)/(T2-T1)+xx02;
    y0=(yy0-yy02)*(T(i)-T1)/(T2-T1)+yy02;
    z0=(zz0-zz02)*(T(i)-T1)/(T2-T1)+zz02;

    %%
    x90=(xx90-xx902)*(T(i)-T1)/(T2-T1)+xx902;
    y90=(yy90-yy902)*(T(i)-T1)/(T2-T1)+yy902;
    z90=(zz90-zz902)*(T(i)-T1)/(T2-T1)+zz902;
    %%
    x_90=(xx_90-xx_902)*(T(i)-T1)/(T2-T1)+xx_902;
    y_90=(yy_90-yy_902)*(T(i)-T1)/(T2-T1)+yy_902;
    z_90=(zz_90-zz_902)*(T(i)-T1)/(T2-T1)+zz_902;
    %%
    xi0=(xxi0-xxi02)*(T(i)-T1)/(T2-T1)+xxi02;
    yi0=(yyi0-yyi02)*(T(i)-T1)/(T2-T1)+yyi02;
    zi0=(zzi0-zzi02)*(T(i)-T1)/(T2-T1)+zzi02;
    %%
    xi90=(xxi90-xxi902)*(T(i)-T1)/(T2-T1)+xxi902;
    yi90=(yyi90-yyi902)*(T(i)-T1)/(T2-T1)+yyi902;
    zi90=(zzi90-zzi902)*(T(i)-T1)/(T2-T1)+zzi902;
    %%
    xi_90=(xxi_90-xxi_902)*(T(i)-T1)/(T2-T1)+xxi_902;
    yi_90=(yyi_90-yyi_902)*(T(i)-T1)/(T2-T1)+yyi_902;
    zi_90=(zzi_90-zzi_902)*(T(i)-T1)/(T2-T1)+zzi_902;
    %%%Same PArt as last time 
    %%Nothing changed in this part 
    xin=FFx(i);
    yin=FFy(i);
    zin=FFz(i);
    if FFx(i)>xi0
        xaccel=(xin- xi0)/(xi90-  xi0);
    else
       xaccel=-(xin-  xi0)/(xi_90-  xi0);
    end 
    if FFy(i)>yi0
        yaccel=(yin-yi0)/(yi90-yi0);
    else
        yaccel=-(yin-yi0)/(yi_90-yi0);
    end 
    if FFz(i)>zi0
        zaccel=(zin- zi0)/(zi90-zi0);
    else
        zaccel=-(zin- zi0)/(zi_90-zi0);
    end 
    xaccel1=sqrt(2)*(xaccel+yaccel)/2;
    yaccel1=sqrt(2)*(xaccel-yaccel)/2;
    zaccel1=zaccel;
    theta(i)=atan(((yaccel1^2+zaccel^2))^0.5/xaccel1)*180/pi;
    if xaccel1<0
        theta(i)=theta(i)+180;
    end 
    if (zaccel>=0)
        phi(i)=atan(yaccel1/zaccel)*180/pi;
    elseif(yaccel1>=0)
        phi(i)=atan(yaccel1/zaccel)*180/pi+180;

    else
        phi(i)=atan(yaccel1/zaccel)*180/pi-180; 
    end 
     %%%
    xin=FFx1(i);
    yin=FFy1(i);
    zin=FFz1(i);
    if FFx1(i)>x0
        xaccel=(xin- x0)/(x90- x0);
    else
        xaccel=-(xin- x0)/(x_90- x0);
    end 
    if FFy1(i)>y0
        yaccel=(yin-y0)/(y90-y0);
    else
        yaccel=-(yin-y0)/(y_90-y0);
    end 
    if FFz1(i)>z0
        zaccel=(zin- z0)/(z90- z0);
        TT(i)=zaccel;
    else
        zaccel=-(zin- z0)/(z_90- z0);

    end 
    if (xaccel>=0)
        phi1(i)=atan(yaccel/xaccel)*180/pi;
    elseif(yaccel>=0)
        phi1(i)=atan(yaccel/xaccel)*180/pi+180;
    else
        phi1(i)=atan(yaccel/xaccel)*180/pi-180;   
    end
     theta1(i)=atan(((yaccel^2+xaccel^2))^0.5/zaccel)*180/pi;
    if zaccel<0
        theta1(i)=theta1(i)+180;
    end 
    %%
    PHiDE=(phi+phi)/2;
    if data_flag == 1
        PHiDE = -PHiDE;
    end
    ThetaDE=(theta+theta)/2;
    %%Kalman filter to be implemented 
end

