clc;clear;
close all;
system('"OpenSees.exe" Von_Misses_truss_Criesfield_Opensees.tcl','-echo')

%

force = load('dispopensees.txt');
% force1 = load('MALTAB.txt');
force2 = load('Disp.txt');
for i = 2
figure(i)
%     img = imread('book.PNG');
%     image('CData',img,'XData',[0 1],'YData',[1.5 0])
        hold on
        
plot(-force2(:,3),0.05*force2(:,1),'DisplayName',"disp")
plot(-force(:,3),force(:,1),'--','DisplayName',"Opensees")
plot(0,0,'Black','DisplayName',"Criesfield book")
%     axis([0 2 0 1])
    legend
    title('Von Misses truss Criesfield')
    xlabel('Disp (m)') 
    ylabel('Force (kN)')  
    grid on
    hold off
end
force = load('reactopensees.txt');
force2 = load('Disp.txt');
for i =2:4
figure(i)
%     img = imread('book.PNG');
%     image('CData',img,'XData',[0 1],'YData',[3 0])
        hold on        
plot(-force2(:,i-1+4),0.05*force2(:,1),'DisplayName',"disp")
plot(-force(:,i-1+4-2),force(:,1),'--','DisplayName',"Opensees")
plot(0,0,'Black','DisplayName',"Criesfield book")
%     axis([0 1 0 1])
    legend
    title('Von Misses truss Criesfield')
    xlabel('Disp (m)') 
    ylabel('Force (kN)')  
    grid on
    hold off
end