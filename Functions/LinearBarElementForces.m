function y = LinearBarElementForces(k,d)
%LinearBarElementForces   This function returns the element nodal 
%                         force vector given the element stiffness  
%                         matrix k and the element nodal displacement 
%                         vector d.
y = k * d;

