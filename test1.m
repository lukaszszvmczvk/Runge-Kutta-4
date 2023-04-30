function [] = test1()

fprintf("Test dla 5 równań różniczkowych stopnia od 1 do 5 porównuje na \n" + ...
    "wykresie ich rozwiązanie analityczne z rozwiązaniem metodą \n" + ...
    "Rungego - Kutty rzędu 4. We wszystkich testach liczba kroków " + ...
    "całkowania \njest ustawiona na 100.\n\n");
pause;

N = 100;

f{1} = {@(x) 5.*exp(-x),@(x) 5, @(x) 1};
f{2} = {@(x) x^2,@(x) 0, @(x) 0,@(x) x};
f{3} = {@(x) sin(x),@(x) 0,@(x) 1, @(x) 0, @(x) 5};
f{4} = {@(x) sin(x),@(x) 0,@(x) 0,@(x) 1,@(x) 4,@(x) 1};
f{5} = {@(x) sin(x).*exp(-x),@(x) 0,@(x) -2,@(x) 0, @(x) 1,@(x) 0, @(x) 1};

Fsolve{1} = @(x) 5./(4*exp(x)) - 5*exp(-5.*x-12)./4;
Fsolve{2} = @(x) x.^3/6 - 11*x/2 + 16/3;
Fsolve{3} = @(x) sqrt(5)*sin(x/sqrt(5)) - 25*cos(x./sqrt(5))/4 ...
    + cos(x)/4 + 6;
Fsolve{4} = @(x) cos(x)./4 + ((26*sqrt(3)+45)/24)*exp(sqrt(3)*x-2*x) - ...
    ((26*sqrt(3)-45)/24)*exp(-sqrt(3)*x-2*x)+x-4;
Fsolve{5} = @(x) sin(sqrt(2)*x)./(12*sqrt(2)) - cos(sqrt(2)*x)./6 + ...
    sin(x)./(10*exp(x)) + cos(x)./(20*exp(x)) - (2/15)*exp(x) + 5/4;

range(1,:) = [-3 5];
range(2,:) = [1 5];
range(3,:) = [0 5];
range(4,:) = [0 5];
range(5,:) = [0 3];

Y0{1} = 0;
Y0{2} = [0 -5];
Y0{3} = [0 1 1];
Y0{4} = [0 0 0 0];
Y0{5} = [1 0 0 0 -1];

eqn{1} = "y' + 5*y = 5*exp(-x)";
eqn{2} = "x*y'' = x^2";
eqn{3} = "5*y^(3) + y' = sin(x)";
eqn{4} = "y^(4) + 4*y^(3) + y'' = sin(x)";
eqn{5} = "y^(5) + y^(3) - 2*y' = sin(x)*exp(-x)";

eqn_sol{1} = "y = 5/(4*exp(x)) - (5/4)*exp(-5x-12)";
eqn_sol{2} = "y = (1/6)*x^3 - (11/2)*x + 16/3";
eqn_sol{3} = "y = sqrt(5)*sin(x/sqrt(5)) - (25/4)*cos(x/sqrt(5))";
eqn_sol{3} = eqn_sol{3} + newline + " + (1/4)*cos(x) + 6";
eqn_sol{4} = "y = (1/4)*cos(x) + (1/24)*(26*sqrt(3)+45)*exp(sqrt(3)*x-2*x)";
eqn_sol{4} = eqn_sol{4} + newline + " - (1/24)*(26*sqrt(3)-45)*exp" + ...
    "(-sqrt(3)*x-2*x) + x - 4";
eqn_sol{5} = "y = (1/(12*sqrt(2))*sin(sqrt(2)*x) - (1/6)*cos(sqrt(2)*x)";
eqn_sol{5} = eqn_sol{5} + newline + "+ sin(x)/(10*exp(x)) + cos(x)/" + ...
    "(20*exp(x)) + (2/15)*exp(x) + 5/4";
for i = 1:5
    clf;
    [X,Y,h] = P1Z24_LSZ_runge_kutta(f{i},range(i,:),Y0{i},N);
    hold on;
    grid on;
    OgY = Fsolve{i}(X);
    plot(X,OgY,'-r',"LineWidth",2);
    plot(X,Y,"b--o");
    title("Równanie nr " + i);
    legend('rozwiązanie analityczne', 'rozwiązanie metodą Rungego - Kutty');
    xlabel('x');
    ylabel('y');

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
    disp(int2str(i) + ". Równanie różniczkowe rzędu: "+int2str(i));
    disp(eqn{i});
    fprintf("Rozw. analityczne: " + eqn_sol{i} + "\n");
    disp("Zakres: " + zakres);
    disp("Y0: " + y0);
    disp("h: " + h);
    fprintf("\n");
    pause;
end

end % function

