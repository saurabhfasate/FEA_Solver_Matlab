clc;clear
close all;
res = load('Disp.txt');
% figure(1)
    img = imread('book.PNG');
    image('CData',img,'XData',[0 1 ],'YData',[2975 0 ])
hold on
plot(-res(:,3),2000/100*res(:,1),'DisplayName','2000 N/m');
% plot(res0(:,1),res0(:,2),'*','DisplayName','AKChopra Book');
%     axis([0 1 0 2000])
    legend
    title('Crisfield Von mises truss')
    xlabel('Disp (m)') 
    ylabel('Force (N)')  
    grid on
    hold off





Stop
% ==============================
% All Plots Compare
% ==============================
res1 = load('Disp_2000.txt');
res2 = load('Disp_1000.txt');
res3 = load('Disp_500.txt');
res4 = load('Disp_250.txt');
% res0 = load('dydisp_AKChopra_Result.txt');
% figure(1)
    img = imread('book.PNG');
    image('CData',img,'XData',[0 1 ],'YData',[3000 0 ])
hold on
plot(-res1(:,3),2000/100*res1(:,1),'DisplayName','2000 N/m');
plot(-res2(:,3),2000/100*res2(:,1),'DisplayName','1000 N/m');
% plot(-res3(:,3),2000/100*res3(:,1),'DisplayName','500');
plot(-res4(:,3),2000/100*res4(:,1),'DisplayName','500 N/m');
% plot(res(:,2),5*res(:,1),'DisplayName','x');
% plot(res0(:,1),res0(:,2),'*','DisplayName','AKChopra Book');
    axis([0 1 0 2000])
    legend
    title('Crisfield Von mises truss')
    xlabel('Disp (m)') 
    ylabel('Force (N)')  
    grid on
    hold off




