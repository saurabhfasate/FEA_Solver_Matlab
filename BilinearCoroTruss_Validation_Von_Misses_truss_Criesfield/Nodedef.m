classdef Nodedef < handle
   properties 
      x
      y
      z
      x_curr
      y_curr
      z_curr
      mx =0
      my =0
      mz =0
      kx =0
      ky =0
      kz =0
   end
   methods
    function obj = Nodedef(x,y,z)
        obj.x = x;
        obj.y = y;
        obj.z = z;
        obj.x_curr = x;
        obj.y_curr = y;
        obj.z_curr = z;
    end
    function obj = Nodalmass(obj,mx,my,mz)
        obj.mx = mx;
        obj.my = my;
        obj.mz = mz;    
    end
    function obj = Nodalstiff(obj,kx,ky,kz)
        obj.kx = kx;
        obj.ky = ky;
        obj.kz = kz;   
    end
    function M = getnodalMassMatrix(obj)
        M = [obj.mx 0 0;
             0 obj.my 0;
             0 0 obj.mz];
    end
    function K = getNodalStiffMatrix(obj)
        K = [obj.kx 0 0;
             0 obj.ky 0;
             0 0 obj.kz];
    end
   end
end


