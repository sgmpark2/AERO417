function y = PlaneTrussElementStress(E,L,theta,d)
    %PlaneTrussElementStress   This function returns the element stress
    %                          given the modulus of elasticity E, the 
    %                          the length L, the angle theta (in 
    %                          degrees), and the element nodal 
    %                          displacement vector d.
    x = theta * pi/180;
    C = cos(x);
    S = sin(x);

    y = -E/L*[-C -S C S]* d;
end

