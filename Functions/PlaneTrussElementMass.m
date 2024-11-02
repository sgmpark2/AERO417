function y = PlaneTrussElementMass(rho,A,L)
    %PlaneTrussElementStiffness   This function returns the element 
    %                             lumped mass matrix for a plane truss   
    %                             element with a Density rho,  
    %                             cross-sectional area A, length L, and
    %                             The size of the element mass 
    %                             matrix is 4 x 4.


    y = (rho * A * L / 2) * [1, 0,0,0; 
                                   0, 1,0,0; 
                                   0,0,1, 0; 
                                   0, 0,0,1];


end