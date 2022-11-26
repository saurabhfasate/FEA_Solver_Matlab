classdef BilinearCoroTruss < handle
   properties 
      I = 0
      J = 0
      K = 0
      ey = 0
      b = 0
      fy = 0
      f_trial = 0
      f_I = 0
      A = 0
      E = 0
      L0 = 0
   end
   methods
    function obj = BilinearCoroTruss(I,J,A,E,L0,epsy,b)
        obj.I = I;
        obj.J = J;
        obj.A = A;
        obj.E = E;
        obj.L0 = L0;
        obj.K = obj.A*obj.E/obj.L0;
        obj.ey = epsy;
        obj.b = b;
        obj.fy = obj.K*obj.ey;
        obj.f_trial = 0;
        obj.f_I = 0;
    end
    function K = getElementGSM(obj,Node,dr)
        [z,T,lf] = getTranformationMatrix(obj,Node);
        k = getLocalStiffnessMatrix(obj,Node,dr);
        K = T*k*T.' + obj.f_trial/lf*z*z.';
    end
    function k = getLocalStiffnessMatrix(obj,Node,dr)
        T = getTranformationMatrix(obj,Node);

        obj.f_trial = obj.f_trial+T.'*obj.K*dr;
        if(abs(obj.f_trial)> obj.fy)
            obj.f_trial = obj.f_trial/abs(obj.f_trial) * obj.fy;
            k = obj.b*obj.K;
        else
            k = obj.K;
        end
    end
    function [z,T,lf] = getTranformationMatrix(obj,Node)
        i = [Node(obj.I).x_curr,Node(obj.I).y_curr,Node(obj.I).z_curr];
        j = [Node(obj.J).x_curr,Node(obj.J).y_curr,Node(obj.J).z_curr];
        lf = norm(i-j);
        l = (Node(obj.J).x_curr-Node(obj.I).x_curr)/lf;
        m = (Node(obj.J).y_curr-Node(obj.I).y_curr)/lf;
        n = (Node(obj.J).z_curr-Node(obj.I).z_curr)/lf; 
        T = [-l -m -n l m n].';
        z = [n 0 -l -n 0 l].';
    end
    function M = getElementMassMatrix(obj)
      M = zeros(6);
    end
    function F = getEndForce(obj,Node,dr)
        T = getTranformationMatrix(obj,Node);

        obj.f_I= obj.f_I + T.'*obj.K*dr;
        if(abs(obj.f_I)> obj.fy)
            obj.f_I = obj.f_I/abs(obj.f_I) * obj.fy;
        end
        F = T*obj.f_I;
        
          
%         i0 = [Node(obj.I).x,Node(obj.I).y,Node(obj.I).z];
%         j0 = [Node(obj.J).x,Node(obj.J).y,Node(obj.J).z];
%         l0 = norm(i0-j0);
%         
%         i = [Node(obj.I).x_curr,Node(obj.I).y_curr,Node(obj.I).z_curr];
%         j = [Node(obj.J).x_curr,Node(obj.J).y_curr,Node(obj.J).z_curr];
%         lf = norm(i-j);
%         stretch = lf-l0
%         if(abs(stretch)<0.75)
%             AxialForce = obj.K*stretch;
%         else
%             AxialForce = sign(stretch)*obj.fy;            
%         end
%         AxialForce
%         T = getTranformationMatrix(obj,Node);
%         F = T*AxialForce;
    end
   end
end

% 
% %--------------------------------------------------------------------------
% % bi_linear model function
% %--------------------------------------------------------------------------    
% function stress=bi_linear(e,fy,E_elas,b)
%     E_lin =   b  * E_elas;
%     E_epp = (1-b)* E_elas;
%     fy_epp = (1-b)*fy;
%     
%     Et=E_epp+E_lin;
%     
%     stress = zeros(1,length(e));     % Initialise stress
%     stress_epp = zeros(1,length(e));     % Initialise stress
%     stress_lin = zeros(1,length(e));     % Initialise stress
%     
%     for i=2:length(e)
%         de = e(i)-e(i-1);
%         stress_lin(i) = stress_lin(i-1)+ E_lin*de;
%         stress_trial = stress_epp(i-1) + E_epp * de;
%         if abs(stress_trial)>=fy_epp
%             Et = E_lin;
%             stress_epp(i) = sign(stress_trial) * fy_epp;
%         elseif abs(stress_trial)<fy_epp
%             Et = E_epp + E_lin;
%             stress_epp(i) = stress_trial;
%         end
%        stress(i) = stress_epp(i) + stress_lin(i);
%     end
%     return
% end
% function K = getElementGSM(obj,Node,dr)
%         T = getTranformationMatrix(obj,Node);
%         k = getLocalStiffnessMatrix(obj,Node,dr);
%         K = T*k*T.';
%     end
%     function k = getLocalStiffnessMatrix(obj,Node,dr)
%         T = getTranformationMatrix(obj,Node);
% 
%         obj.f_trial = obj.f_trial+T.'*obj.K*dr;
%         if(abs(obj.f_trial)> obj.fy)
%             obj.f_trial = obj.f_trial/abs(obj.f_trial) * obj.fy;
%             k = obj.b*obj.K;
%         else
%             k = obj.K;
%         end
%     end
%     function T = getTranformationMatrix(obj,Node)
%         i = [Node(obj.I).x_curr,Node(obj.I).y_curr,Node(obj.I).z_curr];
%         j = [Node(obj.J).x_curr,Node(obj.J).y_curr,Node(obj.J).z_curr];
%         lf = norm(i-j);
%         l = (Node(obj.J).x_curr-Node(obj.I).x_curr)/lf;
%         m = (Node(obj.J).y_curr-Node(obj.I).y_curr)/lf;
%         n = (Node(obj.J).z_curr-Node(obj.I).z_curr)/lf; 
%         T = -[l m n -l -m -n].';
%     end
%     function M = getElementMassMatrix(obj)
%       M = zeros(6);
%     end
%     function F = getEndForce(obj,Node,dr)
%         T = getTranformationMatrix(obj,Node);
% 
%         obj.f_I= obj.f_I + T.'*obj.K*dr;
%         if(abs(obj.f_I)> obj.fy)
%             obj.f_I = obj.f_I/abs(obj.f_I) * obj.fy;
%         end
%         F = T*obj.f_I;
%     end
%    end
% end
% 
% % 
% % %--------------------------------------------------------------------------
% % % bi_linear model function
% % %--------------------------------------------------------------------------    
% % function stress=bi_linear(e,fy,E_elas,b)
% %     E_lin =   b  * E_elas;
% %     E_epp = (1-b)* E_elas;
% %     fy_epp = (1-b)*fy;
% %     
% %     Et=E_epp+E_lin;
% %     
% %     stress = zeros(1,length(e));     % Initialise stress
% %     stress_epp = zeros(1,length(e));     % Initialise stress
% %     stress_lin = zeros(1,length(e));     % Initialise stress
% %     
% %     for i=2:length(e)
% %         de = e(i)-e(i-1);
% %         stress_lin(i) = stress_lin(i-1)+ E_lin*de;
% %         stress_trial = stress_epp(i-1) + E_epp * de;
% %         if abs(stress_trial)>=fy_epp
% %             Et = E_lin;
% %             stress_epp(i) = sign(stress_trial) * fy_epp;
% %         elseif abs(stress_trial)<fy_epp
% %             Et = E_epp + E_lin;
% %             stress_epp(i) = stress_trial;
% %         end
% %        stress(i) = stress_epp(i) + stress_lin(i);
% %     end
% %     return
% end