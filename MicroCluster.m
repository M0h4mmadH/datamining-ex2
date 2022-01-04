classdef MicroCluster < matlab.mixin.Copyable
    properties
        n
        r
        c
    end
    methods
        function mc = MicroCluster(number,radius,centroid)
            mc.n = number;
            mc.r=radius;
            mc.c=centroid;
        end
    end
end