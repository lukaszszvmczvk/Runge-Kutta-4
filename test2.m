function [] = test2()

fprintf("Test sprawdza wartość błędu globalnego dla 7 równań " + ...
    "różniczkowych,\nktórych rozwiązaniami są wielomiany kolejno stopnia" + ...
    "\nod 1 do 7. Wartości rozwiązań są przybliżane na przedziale długości " + ...
    "1" + ...
    "\ndla N = 10, czyli h = 0.1.\n\n")
pause;
range = zeros(7,2);
range(1,:) = [1 2];
range(2,:) = [-1 0];
range(3,:) = [2 3];
range(4,:) = [0 1];
range(5,:) = [2 3];
range(6,:) = [0 1];
range(7,:) = [3 4];

N = 10;

f{1} = {@(x) 1,@(x) 0, @(x) 1};
f{2} = {@(x) 1,@(x) 0, @(x) 0,@(x) 1};
f{3} = {@(x) 1,@(x) 0, @(x) 0,@(x) 0, @(x) 1};
f{4} = {@(x) x^2+4*x+3,@(x) 0, @(x) 0 @(x) 1};
f{5} = {@(x) x^2,@(x) 0, @(x) 0,@(x) 0, @(x) 1};
f{6} = {@(x) 1,@(x) 0, @(x) 0,@(x) 0, @(x) 0, @(x) 0, @(x) 0,@(x) 1};
f{7} = {@(x) 4*x^5,@(x) 0, @(x) 0,@(x) 1};

Fsolve{1} = @(x) x + 1;
Fsolve{2} = @(x) x.^2/2 + 2*x + 5/2;
Fsolve{3} = @(x) x.^3/6 - x.^2 +2*x - 1/3;
Fsolve{4} = @(x) x.^4/12 + 2*x.^3/3 + 3*x.^2/2;
Fsolve{5} = @(x) x.^5/60 + 2*x.^2/3-3*x+19/5;
Fsolve{6} = @(x) x.^5/720;
Fsolve{7} = @(x) 2*x.^7/21 - 483*x + 8699/7;

Y0{1} = 2;
Y0{2} = [1 1];
Y0{3} = [1 0 0];
Y0{4} = [0 0];
Y0{5} = [1 1 4];
Y0{6} = [0 0 0 0 0 0];
Y0{7} = [2 3];

eqn{1} = "y' = 1";
eqn{2} = "y'' = 1";
eqn{3} = "y^(3) = 1";
eqn{4} = "y'' = x^2 + 4x + 3";
eqn{5} = "y^(3) = x^2";
eqn{6} = "y^(6) = 1";
eqn{7} = "y'' = 4*x^5";

eqnsol{1} = "y = x + 1";
eqnsol{2} = "y = (1/2)*x^2 + 2*x + 5/2";
eqnsol{3} = "y = (1/6)*x^3 - x^2 + 2x - 1/3";
eqnsol{4} = "y = (1/12)*x^4 + (2/3)x^3 + (3/2)x^2";
eqnsol{5} = "y = (1/60)*x^5 + (2/3)*x^2 - 3x + 19/5";
eqnsol{6} = "y = (1/720)*x^6";
eqnsol{7} = "y = (2/21)*x^7 - 483x + 8699/7";

error = zeros(1,7);
for i = 1:7
    y0 = "[";
    for j = 1:length(Y0{i})
        if j == length(Y0{i})
            y0 = y0 + int2str(Y0{i}(j));
        else
            y0 = y0 + int2str(Y0{i}(j))+", ";
        end
    end
    y0 = y0 + "]";
    zakres = "[" + int2str(range(i,1)) + ", " + int2str(range(i,2)) + "]";

    [X,Y] = P1Z24_LSZ_runge_kutta(f{i},range(i,:),Y0{i},N);
    OgY = Fsolve{i}(X);
    error(i) = max(abs(OgY - Y'));
    fprintf("%d. Równanie, którego rozwiązaniem jest wielomian stopnia %d" + ...
        "\n",i,i);
    disp("Równanie: " + eqn{i});
    disp("Y0: " + y0);
    disp("Zakres: " + zakres);
    disp("Rozwiązanie analityczne: " + eqnsol{i});
    fprintf("Błąd globalny = %d\n\n",error(i));
    pause;
end
end

