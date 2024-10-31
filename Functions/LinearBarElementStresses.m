function y = LinearBarElementStresses(k, d, A)
%LinearBarElementStresses   This function returns the element nodal 
%                           stress vector given the element stiffness  
%                           matrix k, the element nodal displacement 
%                           vector d, and the cross-sectional area A.
y = k * d/A;


