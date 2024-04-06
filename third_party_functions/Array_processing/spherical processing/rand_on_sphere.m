function varargout = rand_on_sphere(len, outputType)
arguments
    len (1,1) {mustBeNumeric, mustBeNonnegative, mustBeFinite} = 1;
    outputType (1,1) string {mustBeMember(outputType,["cart","sphere"])} = "sphere"
end
% Author: Tom Shlomo, ACLab BGU, 2020


ph  =  rand(len,1)*2*pi;
v = rand(len,1);
th = acos(2*v-1);

switch outputType
    case "sphere"
        if nargout<=1
            varargout = {[th ph]};
        elseif nargout==2
            varargout = {th,ph};
        else
            error("too many outputs")
        end
    case "cart"
        [x, y, z] = s2c(th, ph);
        if nargout<=1
            varargout = {[x y z]};
        elseif nargout==3
            varargout = {x,y,z};
        else
            error("illegal number of outputs");
        end
end

end

