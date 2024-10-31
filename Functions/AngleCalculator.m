function output = AngleCalculator(x1,y1,x2,y2)
    
    %Calculate x and y length from coordinates
    xlength = (x1-x2);
    ylength = (y1-y2);

    %determine the angle using arctan
    output = atand(ylength/xlength);
    
    %if the value is below 0 sub 180 to make positive
    if output < 0
        output = 180 + output;
    end

end