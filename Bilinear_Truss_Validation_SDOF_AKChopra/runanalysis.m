clc;clear
close all;
res = load('dyForce.txt');
res0 = load('dydisp_AKChopra_Result.txt');
figure(1)
%     img = imread('ThaiandKim0.PNG');
%     image('CData',img,'XData',[0 30],'YData',[0.08 -0.08])
hold on
% plot(res(:,1),res(:,3),'DisplayName','z');
plot(res(:,1),res(:,2),'DisplayName','MATLAB');
plot(res0(:,1),res0(:,2),'*','DisplayName','AKChopra Book');
%     axis([0 30 -0.080 0.080])
    legend
    title('T&K Q5 Example')
    xlabel('Time (s)') 
    ylabel('Force (kN)')  
    grid on
    hold off




