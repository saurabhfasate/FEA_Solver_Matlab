clc;clear;close all
import mlreportgen.dom.*
% =========================================================================
% Take TCL file as an Input
% =========================================================================
dataStr = 'SDOF_AKChopra_Opensees.tcl'; %file (better if it includes full path)
dataStr = strtrim(strsplit(fileread(dataStr),'\n')'); %separate entire file by line
dataStr(cellfun(@isempty,dataStr)) = []; % remove empty lines
dataStr = dataStr(:); %block of profile data (cell of chars)
dataStr = cellfun(@strsplit,dataStr,'UniformOutput',false); %char array broken up into sub-cells

% =========================================================================
% Define Node and Element
% =========================================================================
Ele = {};
fix = [];
for i=1:numel(dataStr)
    s = string(dataStr{i});
    if(s(1)=='#')
        continue
    elseif(s(1)=='node')
        nodetag = str2num(s(2));
        s = str2double(s);
        Node(nodetag) = Nodedef(s(3),s(4),s(5));
    elseif(s(1)=='element')
        if(s(2)=='CatenaryCable')
            eletag = str2num(s(3));
            I = str2num(s(4));
            J = str2num(s(5));
            s = str2double(s);
            w = s(6);
            E = s(7);
            A = s(8);
            S = s(9);
            rho = s(12);
            Tol = 1e-6;
            Nsubsteps = 20;
            L       = [Node(s(5)).x-Node(s(4)).x;Node(s(5)).y-Node(s(4)).y;Node(s(5)).z-Node(s(4)).z;];
            uL      = L;          % updated length
%             Ele(eletag) = CatenaryCableElement(I,J,w,E,A,S,rho,Tol,Nsubsteps,L,uL);
            Ele0 = CatenaryCableElement(I,J,w,E,A,S,rho,Tol,Nsubsteps,L,uL);
            Ele{eletag}=Ele0;
        elseif(s(2)=='Spring')
            eletag = str2num(s(3));
            I = str2num(s(4));
            J = str2num(s(5));
            s = str2double(s);
            K = s(6);
            Ele{eletag} = Spring(I,J,K);
        elseif(s(2).lower=='bilinearspring')
            eletag = str2num(s(3));
            I = str2num(s(4));
            J = str2num(s(5));
            s = str2double(s);
            K = s(6);
            epsy = s(7);
            b = s(8);
            Ele{eletag} = BilinearSpring(I,J,K,epsy,b);
        end
    elseif(s(1).lower=='mass')
        nodetag = str2num(s(2));
        s = str2double(s);
        Node(nodetag) = Nodalmass(Node(nodetag),s(3),s(4),s(5));
    elseif(s(1)=='fix')
        fix = [fix;str2num(s(2)) str2num(s(3)) str2num(s(4)) str2num(s(5))];      
    elseif(s(1)=='load')
        continue
%         Node(str2num(s(2))).px = str2double(s(3));
%         Node(str2num(s(2))).py = str2double(s(4));
%         Node(str2num(s(2))).pz = str2double(s(5));
    else
        printf("Unknown Data")
        disp(s)
    end
end
nn = length(Node);
ne= length(Ele);
Pext = zeros(3*nn,1);
nrn = length(fix(:,1));
restrain = b_c(fix,nn);             % Here restrain is a vector, we got 1 at restrained dof and zeros elsewhere

%--------------------------------------------------------------------------
% Structural DoF
%--------------------------------------------------------------------------
no_dof = 0;
no_rdof = 0;
for i=1:length(restrain)
    if restrain(i)==0
        no_dof = no_dof + 1;
        dof(no_dof) = i;
    else
        no_rdof = no_rdof + 1;
        rdof(no_rdof) = i;
    end
end
%%
% Define Dynamic Load (User Input)
% =======================
Pext(1) = 0;                % (Static Load)
r = 0*Pext;

N = 1;
Maxiter = 100;
Tol = 1e-9;

% [r,R_I,K0,Ele,Node] = NRMethod(Ele,Node,nn,ne,N,Tol,Maxiter,dof,Pext,r);
%_ _ _ _ _ _ _ _ __ _ _
%%
% Define Dynamic Load (User Input)
% =======================
P_dyna = [0            0;
       1.0000e-01   5.0000e+00;
       2.0000e-01   8.6603e+00;
       3.0000e-01   1.0000e+01;
       4.0000e-01   8.6603e+00;
       5.0000e-01   5.0000e+00;
       6.0000e-01   0.0;
       1.0          0.0];

% Natural Time Period
Tn = 1;
z = 0.05; % Damping (Rayleigh)

% Newmark Beta Parameters
gamma = 0.5;
beta = 0.25;

Maxiter = 100;
Tol = 1e-3;

Tmax = 1;
dt   = 0.1;
% _ _ _ _ _ _ _ _ 


N    = Tmax/dt;
for i = 1:N
    P_ramp(i)=interp1(P_dyna(:,1),P_dyna(:,2),dt*i);
end
P_ramp = [0,P_ramp];
PStatic = Pext;

[a0,a1] = RayleighDampingCoeff(Tn,z);
R_I = zeros(3*nn,1);
K0  = zeros(3*nn);
% _ _ _ _ _ _ _ _ 

Newmarkintegrator(Ele,Node,nn,ne,P_ramp,N,dt,a0,a1,K0,R_I,gamma,beta,dof,Maxiter,Tol,PStatic)
% runanalysis
% =====================
% End of Main()
% =====================

%--------------------------------------------------------------------------
% N-R Algorithm
%--------------------------------------------------------------------------
function [r,R_I,Kt,Ele,Node] = NRMethod(Ele,Node,nn,ne,N,Tol,maxiter,dof,Pext,r)

fileID = fopen('Disp.txt','w');
fprintf(fileID,'0 0 0\r\n');
fclose(fileID);
% fileID = fopen('ForceDeformation.txt','w');
% fclose(fileID);

    dR_E = Pext/N;
    R_E = 0*Pext;
    Kt = GSM(Ele,Node,nn,ne,zeros(3*nn,1));            
    R_u = 0*dR_E;
    for n=1:N
        R_E = R_E + dR_E;
        R_u = R_u + dR_E;
        dr_n = 0*r;
        j=1;
        while(j==1 || norm(dr)>Tol)
            [n,j]
            dr = 0*r;
            dr(dof) = Kt(dof,dof)^-1*R_u(dof);
            dr_n = dr_n + dr;
            [Ele,Node] = updatestate(Ele,Node,nn,ne,dr);
            R_I = statedetermination(Ele,Node,nn,ne,dr);
            Kt = GSM(Ele,Node,nn,ne,dr);            
            R_u(dof) = R_E(dof) - R_I(dof);
            j=j+1;
            if j>=maxiter %checking the amount of error at each iteration
                disp("NR Failed_ reached Maximum iteration")
                aaaaaa
                break
            end
        end
        r = r + dr_n;
%         R_I = statedetermination(Ele,nn,ne);
%         Kt = tangGSM(Ele,nn,ne);
                figure(1)
                clf(figure(1))
            for j = 1:nn
                hold on
                plot(Node(j).x_curr,Node(j).z_curr,'--or')
                hold off
            end
            for j = 1:ne
                hold on
                plot([Node(Ele{j}.I).x_curr,Node(Ele{j}.J).x_curr],[Node(Ele{j}.I).z_curr,Node(Ele{j}.J).z_curr])
                hold off
            end
        fileID = fopen('Disp.txt','a');
        fprintf(fileID,'%d %6.6f %6.6f\r\n',n,r(3*2-2),r(3*2));
        fclose(fileID);
    end
end


function  r = Newmarkintegrator(Ele,Node,nn,ne,P_ramp,N,dt,a0,a1,Kt,R_I,gamma,beta,dof,Maxiter,Tol,PStatic)
fileID = fopen('dyforce.txt','w');
fclose(fileID);

% _ _ _ _ _ _ _ _ _ _ _ _ _
% 0.Initialize
M  = GMM(Ele,Node,nn,ne);
C = a0*M+a1*Kt;

r_n   = zeros(3*nn,1);
vel_n = zeros(3*nn,1);
acc_n = zeros(3*nn,1);
% R_I = zeros(3*nn,1);
R_E = zeros(3*nn,1);
R_u = zeros(3*nn,1);
dR_En = zeros(3*nn,1);

    acc_n(dof) = M(dof,dof)^-1*(R_E(dof)+PStatic(dof)-C(dof,dof)*vel_n(dof)-R_I(dof));
    a = 1/beta/dt*M + gamma/beta*C;
    b = 1/beta/2*M + dt*(gamma/2/beta-1)*C;
    K_bar = 1/dt*a;
% _ _ _ _ _ _ _ _ _ _ _ _ _
% 1.For Each Step

    for n=1:N          % For Each Time Step
% _ _ _ _ _ _ _ _ _ _ _ _ _
% 1.1.Initialization
        dR_En(loadeddof) = P_ramp(n+1)-P_ramp(n);
        R_E = R_E + dR_En ;
        R_u = dR_En + a*vel_n + b*acc_n + PStatic ;

        dr = 0*r_n;
        dr_n = 0*r_n;
        dvel_n = 0*dr_n;
        dacc_n = 0*dr_n;
% _ _ _ _ _ _ _ _ _ _ _ _ _
% 1.2.Iterative technique
        j=1;
        while(j==1 || norm(R_u(dof))>Tol)
            K_cap = Kt + K_bar;
            dr(dof) = K_cap(dof,dof)^-1*R_u(dof);
            dr_n = dr_n + dr;
            dvel_n = gamma/beta/dt  *dr_n - gamma/beta   *vel_n + dt*(1-gamma/beta/2)*acc_n;
            dacc_n =     1/beta/dt^2*dr_n -     1/beta/dt*vel_n - 1/beta/2*acc_n;
            [Ele,Node] = updatestate(Ele,Node,nn,ne,dr);
            R_It = statedetermination(Ele,Node,nn,ne,dr);
            Kt = GSM(Ele,Node,nn,ne,dr);            
            R_u = R_E+PStatic - (M*(acc_n + dacc_n) + C*(vel_n+dvel_n) + R_It) ;
            j=j+1;
            if (j==Maxiter) %checking the amount of error at each iteration
                disp("Newmark algorithm Failed, Max iteration reached")
                Newmark_Failed_Max_iteration_reached
            end
        end
        r_n = r_n + dr_n;
        vel_n = vel_n + dvel_n;
        acc_n = acc_n + dacc_n;
        R_I = R_It;   
            figure(1)
                clf(figure(1))
            for j = 1:nn
                hold on
                plot(Node(j).x_curr,Node(j).z_curr,'--or')
                hold off
            end
            for j = 1:ne
                hold on
                plot([Node(Ele{j}.I).x_curr,Node(Ele{j}.J).x_curr],[Node(Ele{j}.I).z_curr,Node(Ele{j}.J).z_curr])
                hold off
            end
            
        fileID = fopen('dyforce.txt','a');
        fprintf(fileID,'%f %6.6f %6.6f \r\n',dt*n,r_n(loadeddof),r_n(3*2));
        fclose(fileID);
    end
end

%--------------------------------------------------------------------------
% get Global Stiffness Matrix
%--------------------------------------------------------------------------
function Kt = GSM(Ele,Node,nn,ne,dr)
Kt=zeros(3*nn);
for i=1:ne
    ele_dof = [3*Ele{i}.I-2,3*Ele{i}.I-1,3*Ele{i}.I-0,3*Ele{i}.J-2,3*Ele{i}.J-1,3*Ele{i}.J-0];     
K = getElementGSM(Ele{i},Node,dr(ele_dof));
    Kt(ele_dof,ele_dof)  = Kt(ele_dof,ele_dof) + K;
end
end
%--------------------------------------------------------------------------
% get Global Mass Matrix
%--------------------------------------------------------------------------
function M = GMM(Ele,Node,nn,ne)
M=zeros(3*nn);
for i=1:ne
    ele_dof = [3*Ele{i}.I-2,3*Ele{i}.I-1,3*Ele{i}.I-0,3*Ele{i}.J-2,3*Ele{i}.J-1,3*Ele{i}.J-0];     
    Mele = getElementMassMatrix(Ele{i});
M(ele_dof,ele_dof)  = M(ele_dof,ele_dof) + Mele;
end
for i=1:nn
    node_dof = [3*i-2,3*i-1,3*i-0];
    Mnode = getnodalMassMatrix(Node(i));
M(node_dof,node_dof)  = M(node_dof,node_dof) + Mnode;
end
end
%--------------------------------------------------------------------------
% get Rayleigh Coeff
%--------------------------------------------------------------------------
function [a0,a1] = RayleighDampingCoeff(Tn,z)
    if(length(Tn) == 1)
        w = 2*pi/Tn;
        a0 = 2*z*w;
        a1 = 2*z/w;
    elseif(length(Tn) == 2)
        w1 = 2*pi/Tn(1);
        w2 = 2*pi/Tn(2);
        a0 = 2*z*w1*w2/(w1+w2);
        a1 = 2*z/(w1+w2);
    end    
end
%--------------------------------------------------------------------------
% Function for Boundary condition
%--------------------------------------------------------------------------    
function [nrl]=b_c(fix,nn)
        nrl = zeros(1,3*nn);
    for j=1:length(fix(:,1))
        jn=fix(j,1);
        xr=fix(j,2);
        yr=fix(j,3);
        zr=fix(j,4);
        nrl(1,3*jn-2)=xr;
        nrl(1,3*jn-1)=yr;
        nrl(1,3*jn)=zr;
    end
return
end
%--------------------------------------------------------------------------
% Update State
%--------------------------------------------------------------------------
function [Ele,Node] = updatestate(Ele,Node,nn,ne,r)
    for j = 1:nn
        Node(j).x_curr   = Node(j).x_curr + r(3*j-2);
        Node(j).y_curr   = Node(j).y_curr + r(3*j-1);
        Node(j).z_curr   = Node(j).z_curr + r(3*j-0);
    end
    for j = 1:ne
        if(isprop(Ele{j},"L"))
        Ele{j}.L       = [Node(Ele{j}.J).x_curr-Node(Ele{j}.I).x_curr;Node(Ele{j}.J).y_curr-Node(Ele{j}.I).y_curr;Node(Ele{j}.J).z_curr-Node(Ele{j}.I).z_curr];
        Ele{j}.uL      = Ele{j}.L;
        Ele{j}.uP      = Ele{j}.P;
        end
    end

end
%--------------------------------------------------------------------------
% get Internal forces
%--------------------------------------------------------------------------
function R_I = statedetermination(Ele,Node,nn,ne,dr)
R_I=zeros(3*nn,1);
for i=1:ne
    ele_dof = [3*Ele{i}.I-2,3*Ele{i}.I-1,3*Ele{i}.I-0,3*Ele{i}.J-2,3*Ele{i}.J-1,3*Ele{i}.J-0];     
    RI  = getEndForce(Ele{i},Node,dr);
    R_I(ele_dof)  = R_I(ele_dof) + RI;
end
end

function val = loadeddof()
    val = 2*3-2;
end