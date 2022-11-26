classdef Spring
   properties 
      I
      J
      K
   end
   methods
    function obj = Spring(Inode,Jnode,K)
        obj.I = Inode;
        obj.J = Jnode;
        obj.K = K;
    end
    function K = getElementGSM(obj,Node,dr)
        T = getTranformationMatrix(obj,Node);
        k = getLocalStiffnessMatrix(obj);
        K = T*k*T.';
    end
    function K = getLocalStiffnessMatrix(obj)
      K = obj.K;
    end
    function T = getTranformationMatrix(obj,Node)
        i = [Node(obj.I).x_curr,Node(obj.I).y_curr,Node(obj.I).z_curr];
        j = [Node(obj.J).x_curr,Node(obj.J).y_curr,Node(obj.J).z_curr];
        l0 = norm(i-j);
        l = (Node(obj.J).x_curr-Node(obj.I).x_curr)/l0;
        m = (Node(obj.J).y_curr-Node(obj.I).y_curr)/l0;
        n = (Node(obj.J).z_curr-Node(obj.I).z_curr)/l0; 
        T = [l m n -l -m -n].';
        
%         T = [l*l  l*m  l*n -l*l -l*m -l*n;
%              m*l  m*m  m*n -m*l -m*m -m*n;
%              n*l  n*m  n*n -n*l -n*m -n*n;
%             -l*l -l*m -l*n  l*l  l*m  l*n;
%             -m*l -m*m -m*n  m*l  m*m  m*n;
%             -n*l -n*m -n*n  n*l  n*m  n*n] 
%         Stop
    end
    function M = getElementMassMatrix(obj)
      M = zeros(6);
    end
    function F = getEndForce(obj,Node,dr)
        i0 = [Node(obj.I).x,Node(obj.I).y,Node(obj.I).z];
        j0 = [Node(obj.J).x,Node(obj.J).y,Node(obj.J).z];
        l0 = norm(i0-j0);
        
        i = [Node(obj.I).x_curr,Node(obj.I).y_curr,Node(obj.I).z_curr];
        j = [Node(obj.J).x_curr,Node(obj.J).y_curr,Node(obj.J).z_curr];
        lf = norm(i-j);
        stretch = lf-l0;
        AxialForce = obj.K*stretch;
%         l = (Node(obj.J).x_curr-Node(obj.I).x_curr)/lf;
%         m = (Node(obj.J).y_curr-Node(obj.I).y_curr)/lf;
%         n = (Node(obj.J).z_curr-Node(obj.I).z_curr)/lf; 
%         T = [l m n -l -m -n].';
        T = getTranformationMatrix(obj,Node);
        F = -T*AxialForce;
    end
   end
end

