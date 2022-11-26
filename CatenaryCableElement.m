classdef CatenaryCableElement < handle
   properties 
%       I
%       J
%       E
%       A
%       w
%       S
%       L
%       uL
%       rho
%       lambda
%       Tol
%       Nsubsteps
%       P
%       uP
%       Ti
%       Tj
%       f_trial
%       f_I
   end
   methods
    function obj = CatenaryCableElement(I,J,w,E,A,S,rho,Tol,Nsubsteps,L,uL)
        obj.I = I;
        obj.J = J;
        obj.w = w;
        obj.E = E;
        obj.A = A;
        obj.S = S;
        obj.rho = abs(w)/9.81;
        obj.Tol = Tol;
        obj.Nsubsteps = Nsubsteps;
        obj.f_trial = 0;
        obj.f_I = 0;
        obj.L = L;
        obj.uL = uL;
    end
    function K = getElementGSM(obj,Node,dr)
        K = getLocalStiffnessMatrix(obj,Node,dr);
    end
    function k = getLocalStiffnessMatrix(obj,Node,dr)
        RI = getEndForce(obj);
        F   = elemflexi(obj);
        k   = inv(F);
        k   = [-k  k ;
               k -k];
    end

    function M = getElementMassMatrix(obj)
      M = 0.5*obj.S*obj.rho*eye(6);
    end
    function RI = getEndForce(obj,Node,dr)
        obj = assumeP(obj);
    j=1;
    dlTol=obj.Tol;
    Maxiter=obj.Nsubsteps;
    while((j==1 || norm(dL) > dlTol))
        obj    = updatelength(obj);
        dL     = obj.L - obj.uL;
        F      = elemflexi(obj);
        dP     = F^-1 * dL;
        obj.P  = obj.P + dP;
        obj.Ti = norm(obj.P);
        obj.Tj = norm([-obj.P(1),-obj.P(2),obj.w*obj.S-obj.P(3)]);
        j=j+1;
        
        if (j==Maxiter)
            disp('WARNING -- Cable Element does not Convergence, Max');
            obj.P
            Cable_Element_does_not_Convergence
        end
    end
        obj    = updatelength(obj);
        RI = [obj.P(1),obj.P(2),obj.P(3),-obj.P(1),-obj.P(2),obj.w*obj.S-obj.P(3)].';
    end
    
    %--------------------------------------------------------------------------
    % Flexibility Matrix
    function F = elemflexi(obj)
    F(1,1) =  -1*obj.S/obj.E/obj.A - 1/obj.w*log((obj.Tj + obj.w*obj.S - obj.P(3))/(obj.Ti - obj.P(3))) + obj.P(1)^2/obj.w*(1/(obj.Ti*(obj.Ti-obj.P(3))) - 1/(obj.Tj*(obj.Tj+obj.S*obj.w-obj.P(3))));
    F(1,2) =  obj.P(1)*obj.P(2)/obj.w*(1/(obj.Ti*(obj.Ti-obj.P(3))) - 1/(obj.Tj*(obj.Tj+obj.S*obj.w-obj.P(3))));
    F(1,3) =  obj.P(1)/obj.w*(1/(obj.Tj) - 1/(obj.Ti));
    F(2,1) =  F(1,2);
    F(2,2) =  -1*obj.S/obj.E/obj.A - 1/obj.w*log((obj.Tj + obj.w*obj.S - obj.P(3))/(obj.Ti - obj.P(3))) + obj.P(2)^2/obj.w*(1/(obj.Ti*(obj.Ti-obj.P(3))) - 1/(obj.Tj*(obj.Tj+obj.S*obj.w-obj.P(3))));
    F(2,3) =  obj.P(2)/obj.w*(1/(obj.Tj) - 1/(obj.Ti));
    F(3,1) =  F(1,3);
    F(3,2) =  F(2,3);
    F(3,3) =  -1*obj.S/obj.E/obj.A - 1/obj.w*((obj.w*obj.S - obj.P(3))/(obj.Tj) + obj.P(3)/(obj.Ti));
    end
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % Update lengths
    function obj = updatelength(obj)
    obj.uL(1) = -obj.P(1)*obj.S/obj.E/obj.A - obj.P(1)/obj.w*log((obj.Tj + obj.w*obj.S - obj.P(3))/(obj.Ti - obj.P(3)));
    obj.uL(2) = -obj.P(2)*obj.S/obj.E/obj.A - obj.P(2)/obj.w*log((obj.Tj + obj.w*obj.S - obj.P(3))/(obj.Ti - obj.P(3)));
    obj.uL(3) = -obj.P(3)*obj.S/obj.E/obj.A + obj.w*obj.S^2/2/obj.E/obj.A + 1/obj.w*(obj.Tj - obj.Ti);
    end
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % Assume force
    function obj = assumeP(obj)
    obj = findlambda(obj);
    obj.P = [-obj.w*obj.L(1)/2/obj.lambda;
             -obj.w*obj.L(2)/2/obj.lambda;
              obj.w/2*(-obj.L(3)*(cosh(obj.lambda)/sinh(obj.lambda)) + obj.S)];

    obj.Ti = norm(obj.P);
    obj.Tj = norm([-obj.P(1),-obj.P(2),obj.w*obj.S-obj.P(3)]);
    end
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % Find Lambda
    function obj = findlambda(obj)
    if (obj.L(1)^2+obj.L(2)^2 == 0)
        obj.lambda=1e6;
        disp("Cable is vertical")
    elseif (obj.S^2 <= obj.L(1)^2+obj.L(2)^2+obj.L(3)^2)
        obj.lambda=0.2;
    elseif (obj.S^2 > obj.L(1)^2+obj.L(2)^2+obj.L(3)^2)
        obj.lambda=sqrt(3*((obj.S^2 - obj.L(3)^2)/(obj.L(1)^2+obj.L(2)^2) -1));
    end
    end
    %--------------------------------------------------------------------------
   end
end
