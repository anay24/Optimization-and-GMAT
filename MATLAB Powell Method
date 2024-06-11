% DEFINE THE OBJECTIVE FUNCTION
func = @(x1,x2) 0.5*8*(sqrt(x1^2 + (10-x2)^2)-10)^2 + 0.5*1*(sqrt(x1^2 + ...
    (10+x2)^2)-10)^2 - 5*x1 - 5*x2;


% INITIAL POINT
loc = [0;
       0];

% INITIAL SEARCH DIRECTION MATRIX
sval = [1 0;   % DEFINE THE INITIAL SEARCH DIRECTION AS THE IDENTITY MATRIX
        0 1];

% COMPUTE THE FIRST MINIMUM USING GETALPHASTAR FUNCTION
loc(:,2) = getalphastar(loc(:,1), func, sval(:,1));

% ITERATE TO FIND THE MINIMUM
for q = 3:1000
    % COMPUTE THE MINIMUM FOR CURRENT SEARCH DIRECTION
    loc(:,q) = getalphastar(loc(:,q-1), func, sval(:,q-1));
    
    % UPDATE THE SEARCH DIRECTION MATRIX
    sval(:,q) = sval(:,q-2);
    
    % COMPUTE THE NEXT MINIMUM USING UPDATED SEARCH DIRECTION
    loc(:,q+1) = getalphastar(loc(:,q), func, sval(:,q));
    
    % UPDATE THE SEARCH DIRECTION MATRIX
    sval(:,q+1) = sval(:,q-1);
    
    % COMPUTE THE NEXT MINIMUM USING THE UPDATED SEARCH DIRECTION
    loc(:,q+2) = getalphastar(loc(:,q+1), func, sval(:,q+1)); 
    
    % UPDATE THE SEARCH DIRECTION TO BE THE NORMALIZED DIFFERENCE BETWEEN TWO CONSECUTIVE MINIMUM POINTS
    sval(:,q+2) = (loc(:,q+2) - loc(:,q)) / norm(loc(:,q+2) - loc(:,q)); 
    
    % CHECK FOR CONVERGENCE
    if abs(((loc(1,q) - loc(1,q-1)) / loc(1,q-1))) < 0.0035 && abs(((loc(2,q) - loc(2,q-1)) / loc(2,q-1))) < 0.0035
        % IF CONVERGENCE CRITERIA MET, BREAK THE LOOP
        totIterations = q;
        break
    end
end 

% PLOT THE CONTOUR OF THE OBJECTIVE FUNCTION
fcontour(func)
hold on 

% PLOT THE PATH OF ITERATIONS
plot(loc(1,:),loc(2,:)) 
hold off

finalpoint = [loc(end-1,end),loc(end,end)];
iterations = totIterations;
% fprintf(['Final Point: ' repmat(' %1.0f ',1,numel(xarrayoutput)) '\n'],xarrayoutput);
% fprintf(['Iterations: ' repmat(' %1.0f ',1,numel(iterations)) '\n'],iterations);

% FUNCTION TO FIND THE OPTIMAL STEP SIZE ALONG A SEARCH DIRECTION
function alpha = getalphastar(alphastar, func, sval)
    % DEFINE THE GOLDEN SECTION RATIO
    goldsec = (3 - sqrt(5)) / 2;
    
    % INITIAL BRACKET POINTS
    al = alphastar;
    alp = alphastar + 5 * sval;
    almin = alphastar - 5 * sval;
    alpmin = alphastar;
    
    % INITIAL POINTS
    al1 = (1 - goldsec) * al + goldsec * alp;
    al2 = goldsec * al + (1 - goldsec) * alp;
    al1min = (1 - goldsec) * almin + goldsec * alpmin;
    al2min = goldsec * almin + (1 - goldsec) * alpmin;
    
    % PERFORM ITERATIONS
    for q2 = 1:10
        % EVALUATE FUNCTION VALUES
        funcone = func(al1(1), al1(2));
        functwo = func(al2(1), al2(2));
        funconemin = func(al1min(1), al1min(2));
        functwomin = func(al2min(1), al2min(2));
        
        % UPDATE BRACKET POINTS
        if funcone > functwo
            al = al1;
            al1 = al2;
            al2 = goldsec * al + (1 - goldsec) * alp;
        else
            alp = al2;
            al2 = al1;
            al1 = (1 - goldsec) * al + goldsec * alp;
        end
        
        if funconemin > functwomin
            almin = al1min;
            al1min = al2min;
            al2min = goldsec * almin + (1 - goldsec) * alpmin;
        else
            alpmin = al2min;
            al2min = al1min;
            al1min = (1 - goldsec) * almin + goldsec * alpmin;
        end
    end
    
    % SELECT THE MINIMUM AMONG THE TWO DIRECTIONS
    functminselect = func(al1(1), al1(2));
    functminselect2 = func(al1min(1), al1min(2));
    
    if functminselect < functminselect2
        alpha = al1;
    else
        alpha = al1min;
    end
end
