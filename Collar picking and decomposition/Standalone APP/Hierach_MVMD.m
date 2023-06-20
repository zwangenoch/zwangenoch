function [u, u_hat, omega] = Hierach_MVMD(vdl, depth)
B=xlsread("CAV_well_data.xlsx");
t=B(:,1);
%xx=B(:,86); %Channel 85
signal=(B(:,2:97))';

windowSize = 8; 
 b = (1/windowSize)*ones(1,windowSize);
a = 1;
signal_filter = filter(b,a,signal);
a=1:6001; 

% If the left side has a strong end-effect, cut-off the end-effect
b=[1:100, 6001];   
c=setdiff(a,b);
signal_filter=signal_filter(:,c);


newmap=newmap();
 % plot the cutoff signal
x=[1:96];
yy=t(c);
figure(1);
imagesc(x, yy, (signal_filter(10:96,:))');
colormap(newmap);
caxis([-0.02, 0.02]);
colorbar('northoutside');
title('Cut-off end and filtered signal', 'FontSize', 20)
hold on;


% Perform intial MVMD for the signal with cut-off end
[u, u_hat, omega] = MVMD_ver1(signal_filter, 2000, 0, 5, 0, 1, 1e-7); 
   figure(2)
    subplot(5,1,1)
    plot(t(c),  signal_filter(60,:))
    title('Filtered signal CH60', 'FontSize', 20)
    hold on
    subplot(5,1,2)
    plot(t(c), u(1,:,60))
     title('IMF1-MVMD', 'FontSize', 20)
    hold on
    subplot(5,1,3)
    plot(t(c), u(2,:,60))
    title('IMF2-MVMD', 'FontSize', 20)
    hold on
     subplot(5,1,4)
    plot(t(c), u(3,:,60))
     title('IMF3-MVMD', 'FontSize', 20)
    hold on
    subplot(5,1,5)
     plot(t(c), u(4,:,60))
      title('IMF4-MVMD', 'FontSize', 20)
    hold on





%% Apply MVMD for u1 to indentify baseline (corrosion signal)
u1 = (squeeze(u(1,:, :)))';  
[uu1, u_hat, omega] = MVMD_ver1(u1, 2000, 0, 3, 0, 1, 1e-7);


  figure(3)
    subplot(4,1,1)
    plot(t(c), u1(60,:))
    title(' MVMD of IMF1, CH60', 'FontSize', 20)
    hold on
    subplot(4,1,2)
    plot(t(c), uu1(1,:,60))
     title('IMF1-MVMD of IMF1', 'FontSize', 20)
    hold on
    subplot(4,1,3)
    plot(t(c), uu1(2,:,60))
    title('IMF2-MVMD of IMF1', 'FontSize', 20)
    hold on
     subplot(4,1,4)
    plot(t(c), uu1(3,:,60))
     title('IMF3-MVMD of IMF1', 'FontSize', 20)
    hold on
   


% Plot the figure of baseline
x = [10:96];
IMF2_1=squeeze(uu1(1,:,:));
yy=t(c);
figure(4);
imagesc(x, yy, IMF2_1(:, 10:96));
colormap(newmap);
caxis([-0.02, 0.02]);
colorbar('northoutside');
title('baseline', 'FontSize', 20)
hold on;


%Apply MVMD for u2 (defined below) to identify the colar signal
u2 = signal_filter-(squeeze(uu1(1,:,:)+u(end,:,:)))';  
cu2=u2;
for i=1:96
    index=find(cu2(i,:)<0);
    cu2(i,index)=0;
end
[cuu2, u_hat, omega] = MVMD_ver1(cu2, 2000, 0, 4, 0, 1, 1e-7);


figure(5)
     subplot(5,1,1)
    plot(t(c), cu2(60,:))
    title('cut off  signal', 'FontSize', 20)
    hold on
     subplot(5,1,2)
    plot(t(c), cuu2(1,:,60))
     title('IMF1', 'FontSize', 20)
    hold on
     subplot(5,1,3)
    plot(t(c), cuu2(2,:,60))
    title('IMF2', 'FontSize', 20)
    hold on
     subplot(5,1,4)
    plot(t(c), cuu2(3,:,60))
     title('IMF3', 'FontSize', 20)
    hold on
     subplot(5,1,5)
     plot(t(c), cuu2(4,:,60))
      title('IMF4', 'FontSize', 20)
    hold on

% cut off cuu2(1,;,:) if it implies the initial colar signal.  In some cases, cuu(2,:,;) may represent the initial colar
% signal. 
colar0=squeeze(cuu2(1,:,:));
for i=1:96
    index=find(colar0(:,i)<0);
    colar0(index, i)=0;
end

yy=t(c);
figure(6);
imagesc(x, yy, colar0(:, 10:96));
colormap(newmap);
caxis([-0.02, 0.02]);
colorbar('northoutside');
title('colar0', 'FontSize', 20)
hold on;






%Hierarchical MVMD to identify colar1 (more clear colar)

[cuu3, u_hat, omega] = MVMD_ver1((colar0)', 2000, 0, 2, 0, 1, 1e-7);

figure(7)
    subplot(3,1,1)
    plot(t(c), colar0(:,60))
    title(' MVMD , CH60', 'FontSize', 20)
    hold on
    subplot(3,1,2)
    plot(t(c), cuu3(1,:,60))
     title('IMF1', 'FontSize', 20)
    hold on
    subplot(3,1,3)
    plot(t(c), cuu3(2,:,60))
    title('IMF2', 'FontSize', 20)
    
%% We think cuu3(2,:,:) represents colar1 signal 
colar1=squeeze(cuu3(2,:,:));
for i=1:96
    index=find(colar1(:,i)<0);
    colar1(index, i)=0;
end
figure(8);
imagesc(x, yy, colar1(:, 10:96));
colormap(newmap);
caxis([-0.02, 0.02]);
colorbar('northoutside');
title(' colar1', 'FontSize', 20)
hold on;


%%% Use one of the following colar signal to take numerical differentiation
figure (9)
subplot(2,1,1)
plot(t(c), colar0(:,60));
 title('Colar0, CH60', 'FontSize', 20);
 hold on
 subplot(2,1,2)
 plot(t(c), colar1(:,60));
 title('Colar1,CH60', 'FontSize', 20);



%% We may choose either colar0 or colar1 as the final colar signal. This choice depends on the specific case. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


