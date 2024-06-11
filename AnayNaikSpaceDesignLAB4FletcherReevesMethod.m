% DEFINE THE OBJECTIVE FUNCTION
syms x1 x2
func = 0.5*8*(sqrt(x1^2 + (10-x2)^2)-10)^2 + 0.5*1*(sqrt(x1^2 + ...
    (10+x2)^2)-10)^2 - 5*x1 - 5*x2;

% INITIAL POINT
loc(:,1) = [0; % INITIALIZE THE LOCATION VECTOR WITH STARTING POINT
            0];

betaval = []; % INITIALIZE BETA VALUES VECTOR
iter = 2;
for q = 2:1000
syms x1 x2 

     if q == 2
    locone = [0; % COPY OF INITIAL POINT
          0];
    grad = double(subs(gradient(func,[x1,x2]),[x1;x2],locone)); % CALCULATE GRADIENT AT INITIAL POINT
    oggrad = grad; % STORE ORIGINAL GRADIENT
    sval(:,1) = -(grad/norm(grad)); % INITIAL SEARCH DIRECTION
    a = double(dot(grad,grad)); % COMPUTE INITIAL VALUE FOR 'a'
    [alph, ast] = getalphastar(locone,sval); % COMPUTE INITIAL ALPHA AND ASTAR
    second = alph(:,1); % STORE SECOND POINT

    loc(:,2) = second; % UPDATE LOCATION VECTOR WITH SECOND POINT
     continue
     end

    grad = double(subs(gradient(func,[x1 x2]),[x1;x2],loc(:,q-1))); % CALCULATE GRADIENT AT CURRENT POINT
    bval = dot(grad,grad); % COMPUTE BETA VALUE
    betaval(:,q-2) = bval/a; % STORE BETA VALUE
    a = bval; % UPDATE 'a' VALUE

    sval(:,q-1) = (-grad + betaval(:,q-2)*sval(:,q-2))/norm(-grad + betaval(:,q-2)*sval(:,q-2)); % COMPUTE SEARCH DIRECTION
    [alphagain,ast] = getalphastar(loc(:,q-1),sval(:,q-1)); % COMPUTE ALPHA AND ASTAR
    loc(:,q) = alphagain; % STORE NEW LOCATION

    slope = dot(sval(:,q-1),grad); % CALCULATE SLOPE

  if slope >= 0
     sval(:,q) = -(oggrad/norm(oggrad)); % RESET SEARCH DIRECTION
     disp("Slope has returned at one point")
  elseif ast == 0 
     break % TERMINATE IF ASTAR IS 0
   % CHECK FOR CONVERGENCE
  elseif abs(((loc(1,q) - loc(1,q-1)) / loc(1,q-1))) < 0.0035 && abs(((loc(2,q) - loc(2,q-1)) / loc(2,q-1))) < 0.0035
        % IF CONVERGENCE CRITERIA MET, BREAK THE LOOP
        totIterations = iter; % STORE TOTAL ITERATIONS
        break
  end
  iter = iter + 1;
end 

 %PLOT THE CONTOUR OF THE OBJECTIVE FUNCTION
 fcontour(func) % PLOT THE CONTOUR OF THE OBJECTIVE FUNCTION
 hold on 

% PLOT THE PATH OF ITERATIONS
plot(loc(1,:),loc(2,:)) % PLOT THE PATH OF ITERATIONS
hold off

finalpoint = [loc(end-1,end),loc(end,end)]; % FINAL POINT
iterations = totIterations; % TOTAL ITERATIONS
% fprintf(['Final Point: ' repmat(' %1.0f ',1,numel(xarrayoutput)) '\n'],xarrayoutput);
% fprintf(['Iterations: ' repmat(' %1.0f ',1,numel(iterations)) '\n'],iterations);


function [alph,ast] = getalphastar(alphastar, sval)

% DEFINE THE OBJECTIVE FUNCTION
func = @(x1,x2) 0.5*8*(sqrt(x1^2 + (10-x2)^2)-10)^2 + 0.5*1*(sqrt(x1^2 + ...
    (10+x2)^2)-10)^2 - 5*x1 - 5*x2;

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
        alph = al1;
        ast = norm(alph-alphastar);

    else
        alph = al1min;
        ast = norm(alph-alphastar);
    end
end
