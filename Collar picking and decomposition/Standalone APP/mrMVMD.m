function MVMDdata = decomp_MVMD(vdl, CH_select)

clc; clear all; close all, clf
B=xlsread("CAV_well_data.xlsx");
t=B(:,1);
%xx=B(:,86); %Channel 85
signal=(B(:,2:97))';

windowSize = 8; 
 b = (1/windowSize)*ones(1,windowSize);
a = 1;
signal_filter = filter(b,a,signal);
a=1:6001;  
b=[6001];
c=setdiff(a,b);
signal_filter=signal_filter(:,c);
%signal=signal(:,c);



[u, u_hat, omega] = MVMD_ver1(signal_filter, 2000, 0, 5, 0, 1, 1e-7); 
   figure(1)
    subplot(5,1,1)
    plot(t(c),  signal_filter(40,:))
    title('Filtered signal CH85', 'FontSize', 20)
    hold on
    subplot(5,1,2)
    plot(t(c), u(1,:,40))
     title('IMF1-MVMD', 'FontSize', 20)
    hold on
    subplot(5,1,3)
    plot(t(c), u(2,:,40))
    title('IMF2-MVMD', 'FontSize', 20)
    hold on
     subplot(5,1,4)
    plot(t(c), u(3,:,40))
     title('IMF3-MVMD', 'FontSize', 20)
    hold on
    subplot(5,1,5)
     plot(t(c), u(4,:,40))
      title('IMF4-MVMD', 'FontSize', 20)
    hold on





 

%% Extract baseline
u1 = (squeeze(u(1,:, :)+u(2,:, :)))';  
[uu1, u_hat, omega] = MVMD_ver1(u1, 2000, 0, 4, 0, 1, 1e-7);

%mru1=cat(2,uu1,uu2);

figure(2)
    subplot(5,1,1)
    plot(t(c), u1(40,:))
    title('IMF1 of MVMD, CH85', 'FontSize', 20)
    hold on
    subplot(5,1,2)
    plot(t(c), uu1(1,:,40))
     title('IMF1-MVMD of IMF1', 'FontSize', 20)
    hold on
    subplot(5,1,3)
    plot(t(c), uu1(2,:,40))
    title('IMF2-MVMD of IMF1', 'FontSize', 20)
    hold on
     subplot(5,1,4)
    plot(t(c), uu1(3,:,40))
     title('IMF3-MVMD of IMF1', 'FontSize', 20)
    hold on
    subplot(5,1,5)
     plot(t(c), uu1(4,:,40))
      title('IMF4-MVMD of IMF1', 'FontSize', 20)
    hold on

newmap=newmap();

x = [10:96];

IMF2_1=squeeze(uu1(1,:,:));
yy=t(c);
figure(3);
imagesc(x, yy, IMF2_1(:, 10:96));
colormap(newmap);
caxis([-0.02, 0.02]);
colorbar('northoutside');
title('baseline', 'FontSize', 20)
hold on;


%% Extract initial colar
u2 = signal_filter-(squeeze(uu1(1,:,:)+u(end,:,:)))';  
[uu2, u_hat, omega] = MVMD_ver1(u2, 2000, 0, 4, 0, 1, 1e-7);

cu2=u2;
for i=1:96
    index=find(cu2(i,:)<0);
    cu2(i,index)=0;
end
[cuu2, u_hat, omega] = MVMD_ver1(cu2, 2000, 0, 4, 0, 1, 1e-7);


figure(4)
    subplot(5,2,1)
    plot(t(c), u2(40,:))
    title('IMF1 of MVMD, CH85', 'FontSize', 20)
    hold on
     subplot(5,2,2)
    plot(t(c), cu2(40,:))
    title('IMF1 of MVMD, CH85', 'FontSize', 20)
    hold on

    subplot(5,2,3)
    plot(t(c), uu2(1,:,40))
     title('IMF1-MVMD of IMF1', 'FontSize', 20)
    hold on
     subplot(5,2,4)
    plot(t(c), cuu2(1,:,40))
     title('IMF1-MVMD of IMF1', 'FontSize', 20)
    hold on
    subplot(5,2,5)
    plot(t(c), uu2(2,:,40))
    title('IMF2-MVMD of IMF1', 'FontSize', 20)
    hold on
     subplot(5,2,6)
    plot(t(c), cuu2(2,:,40))
    title('IMF2-MVMD of IMF1', 'FontSize', 20)
    hold on
   subplot(5,2,7)
    plot(t(c), uu2(3,:,40))
     title('IMF3-MVMD of IMF1', 'FontSize', 20)
    hold on
     subplot(5,2,8)
    plot(t(c), cuu2(3,:,40))
     title('IMF3-MVMD of IMF1', 'FontSize', 20)
    hold on
    subplot(5,2,9)
     plot(t(c), uu2(4,:,40))
      title('IMF4-MVMD of IMF1', 'FontSize', 20)
    hold on
     subplot(5,2,10)
     plot(t(c), cuu2(4,:,40))
      title('IMF4-MVMD of IMF1', 'FontSize', 20)
    hold on



IMF2_1=squeeze(uu2(1,:,:));
yy=t(c);
figure(5);
imagesc(x, yy, IMF2_1(:, 10:96));
colormap(newmap);
caxis([-0.02, 0.02]);
colorbar('northoutside');
title('colar', 'FontSize', 20)
hold on;

IMF2_1=squeeze(cuu2(1,:,:));
yy=t(c);
figure(6);
imagesc(x, yy, IMF2_1(:, 10:96));
colormap(newmap);
caxis([-0.02, 0.02]);
colorbar('northoutside');
title('colar', 'FontSize', 20)
hold on;



%% Refine colar
u3=u2-(squeeze(uu2(1,:,:)+uu2(end,:,:)))';
[uu3, u_hat, omega] = MVMD_ver1(u3, 2000, 0, 4, 0, 1, 1e-7);

cu3=cu2-(squeeze(cuu2(1,:,:)+cuu2(end,:,:)))';
for i=1:96
    index=find(cu3(i,:)<0);
    cu3(i,index)=0;
end
%cu3=(filter(b,a,cu3'))';


[cuu3, u_hat, omega] = MVMD_ver1(cu3, 2000, 0, 4, 0, 1, 1e-7);


figure(7)
subplot(2,1,1)
    plot(t(c), uu3(1,:,40))
     title('IMF1', 'FontSize', 20)
    hold on
subplot(2,1,2)
    plot(t(c), cuu3(1,:,40))
     title('IMF1', 'FontSize', 20)
    hold on







IMF2_1=squeeze(uu3(1,:,:));
figure(8);
imagesc(x, yy, IMF2_1(:, 10:96));
colormap(newmap);
caxis([-0.02, 0.02]);
colorbar('northoutside');
title('colar', 'FontSize', 20)
hold on;


IMF2_1=squeeze(cuu3(1,:,:));
figure(9);
imagesc(x, yy, IMF2_1(:, 10:96));
colormap(newmap);
caxis([-0.02, 0.02]);
colorbar('northoutside');
title('colar', 'FontSize', 20)
hold on;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


