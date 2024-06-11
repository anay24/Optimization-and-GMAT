
% INITIALIZE PARAMETERS
rpp = 1; % Initialize the parameter rpp
whilefunc = 0; % Initialize the condition for the while loop
amountFunc = 0; % Initialize the counter for function calls

% DEFINE FUNCTIONS
first = @(x1,x2) x1 + x2 - 1; % Define the first function
second = @(x1) x1; % Define the second function

initPoint = [0;0]; % Initialize the starting point

while whilefunc == 0 % Begin the while loop
    counterfunc = 0; % Initialize the counter for function evaluations
    initval = 0.011; % Initialize a value
    
    % Evaluate functions at the initial point
    firstc =  first(initPoint(1),initPoint(2)); 
    secc = second(initPoint(1));

    % Determine the appropriate function based on conditions
    if firstc > (-0.22) && secc > (-0.22)
        thefunc = @(x1,x2) (x1-1)^2 + (x2-1)^2 + rpp*((-1*(2*(-0.22)-(x1 + x2 -1))/(-0.22)^2)+(-1*(2*(-0.22)-x1)/(-0.22)^2));
    end
    if firstc <= (-0.22) && secc <= (-0.22)
        thefunc = @(x1,x2) (x1-1)^2 + (x2-1)^2 + rpp*((-1/(x1+x2-1))+(-1/x1));
    end
    if firstc <= (-0.22) && secc > (-0.22)
        thefunc = @(x1,x2) (x1-1)^2 + (x2-1)^2 + rpp*((-1/(x1+x2-1))+(-1*(2*(-0.22)-x1)/(-0.22)^2));
    end
    if firstc > (-0.22) && secc <= (-0.22)
        thefunc = @(x1,x2) (x1-1)^2 + (x2-1)^2 + rpp*((-1*(2*(-0.22)-(x1 + x2 -1))/(-0.22)^2)+(-1/x1));
    end

    % Perform function evaluation at the initial point
    thepoint = thefunc(initPoint(1),initPoint(2));
    endall = thepoint;

    % Compute the gradient of the function
    syms x1 x2
    functiongrad = gradient(thefunc(x1,x2));
    functiongrad = subs(functiongrad,[x1,x2],[initPoint(1),initPoint(2)]);
    ai = dot(functiongrad,functiongrad);

    % BEGIN LINE SEARCH
    whilefunc1 = false;
    whilefunc0 = false;

    while whilefunc1 == false
        Sval = -functiongrad;
        whilefunc2 = false;
        whilefunc3 = false;
        whilefunc4 = 0;
        compare = 0;

        while whilefunc2 == 0
            if whilefunc3 == 0
                al = 0;
                al2 = 1;
                lp = 0.38197*(al2-al)+al;
                lp2 = 1.61803*lp;
            end
            whilefunc3 = true;
            if whilefunc4 == 0
                funcl = thefunc(initPoint(1)+al*Sval(1),initPoint(2)+al*Sval(2));
                counterfunc = counterfunc+1;
                compare = compare+1;
                funcl2 = thefunc(initPoint(1)+al2*Sval(1),initPoint(2)+al2*Sval(2));
                counterfunc = counterfunc+1;
                compare = compare+1;
                funcl3 = thefunc(initPoint(1)+lp*Sval(1),initPoint(2)+lp*Sval(2));
                counterfunc = counterfunc+1;
                compare = compare+1;
                funcl4 = thefunc(initPoint(1)+lp2*Sval(1),initPoint(2)+lp2*Sval(2));
                counterfunc = counterfunc+1;
                compare = compare+1;
            end
            if whilefunc4 == 1
                funcl2 = funcl4;
                funcl4 = funcl3;
                funcl3 = thefunc(initPoint(1)+lp*Sval(1),initPoint(2)+lp*Sval(2));
                counterfunc = counterfunc+1;
                compare = compare+1;
            end
            if whilefunc4 == 2
                funcl = funcl3;
                funcl3 = funcl4;
                funcl4 = thefunc(initPoint(1)+lp2*Sval(1),initPoint(2)+lp2*Sval(2));
                counterfunc = counterfunc+1;
                compare = compare+1;
            end

            % Update the comparison variables
            if funcl3 < funcl4
                al2 = lp2;
                lp2 = lp;
                lp = al2+al-lp2;
                whilefunc4 = 1;
                ast = lp2;
            else
                al = lp;
                lp = lp2;
                lp2 = al2+al-lp;
                whilefunc4 = 2;
                ast = lp;
            end
            if compare >= -2*log(1)+3
                whilefunc2 = true;
            end
        end

        % Determine the optimal step size
        if funcl3 < funcl4
            alpha_3 = lp2;
            lp2 = lp;
            lp = al;
        else
            alpha_3 = al2;
        end

        % Compute values for convergence
        subtrac = thefunc(initPoint(1)+lp*Sval(1),initPoint(2)+lp*Sval(2));
        counterfunc = counterfunc+1;
        suctract = thefunc(initPoint(1)+lp2*Sval(1),initPoint(2)+lp2*Sval(2));
        counterfunc = counterfunc+1;
        subtractit = thefunc(initPoint(1)+alpha_3*Sval(1),initPoint(2)+alpha_3*Sval(2));
        counterfunc = counterfunc+1;
        t1 = (((subtractit-subtrac)/(alpha_3-lp))-((suctract-subtrac)/(lp2-lp)))/(alpha_3-lp2);
        t = ((suctract-subtrac)/(lp2-lp))-t1*(lp+lp2);
        ast = -t/(2*t1);
        if isnan(ast) == true
            ast = lp2;
        end

        % Check for convergence
        if abs(ast) < 0.001
            whilefunc1 = true;
            whilefunc0 = true;
        end

        % Update the initial point
        initPoint = initPoint+ast*Sval;

        % CHECK FOR FUNCTION
        while whilefunc0 == false
            firstc =  first(initPoint(1),initPoint(2));
            secc = second(initPoint(1));

            % Determine the appropriate function based on conditions
            if firstc > (-0.22) && secc > (-0.22)
                thefunc = @(x1,x2) (x1-1)^2 + (x2-1)^2 + rpp*((-1*(2*(-0.22)-(x1 + x2 -1))/(-0.22)^2)+(-1*(2*(-0.22)-x1)/(-0.22)^2));
            end
            if firstc <= (-0.22) && secc <= (-0.22)
                thefunc = @(x1,x2) (x1-1)^2 + (x2-1)^2 + rpp*((-1/(x1+x2-1))+(-1/x1));
            end
            if firstc <= (-0.22) && secc > (-0.22)
                thefunc = @(x1,x2) (x1-1)^2 + (x2-1)^2 + rpp*((-1/(x1+x2-1))+(-1*(2*(-0.22)-x1)/(-0.22)^2));
            end
            if firstc > (-0.22) && secc <= (-0.22)
                thefunc = @(x1,x2) (x1-1)^2 + (x2-1)^2 + rpp*((-1*(2*(-0.22)-(x1 + x2 -1))/(-0.22)^2)+(-1/x1));
            end

            % Compute the gradient
            functiongrad = gradient(thefunc(x1,x2));
            counterfunc = counterfunc+1;
            functiongrad = subs(functiongrad,[x1,x2],[initPoint(1),initPoint(2)]);
            bval = dot(functiongrad,functiongrad);
            Beta = bval/ai;
            Sval = -functiongrad+Beta*Sval;
            ai = bval;
            slope = dot(Sval,functiongrad);
            if slope >= 0
                whilefunc0 = true;
                whilefunc1 = true;
            end
            f_check = 0;
            endall = f_check;
        end
    end

    % UPDATE FUNCTION COUNTER
    amountFunc = amountFunc + counterfunc;

    % BREAK CONDITION FOR THE WHILE LOOP
    if first(initPoint(1),initPoint(2)) == 0 || first(initPoint(1),initPoint(2)) < 0 && second(initPoint(1)) == 0 || second(initPoint(1)) > 0
        break
    end
    rpp = rpp + 1; % Increment the parameter rpp
end

% DISPLAY FUNCTION CALLS
disp('function Calls:')
disp(counterfunc)
