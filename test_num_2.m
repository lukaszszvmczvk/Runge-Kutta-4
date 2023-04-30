function [] = test_num_2()
fprintf("Test porównuje dokładność metody Rungego - Kutty rzędu 4-go\n" + ...
    "(wzór 3/8) z tą samą metodą określoną wzorem 'klasycznym' oraz\n" + ...
    "z metodą Adamsa - Bashfortha - Moultona rzędu 2-go. Zestawiane\nsą " + ...
    "ze sobą " + ...
    "błędy globalne tych metod dla 5 równań różniczkowych.\nRównania są " + ...
    "kolejno rzędu od 1 do 5. Ich rozwiązanie jest przybliżane\nna " + ...
    "przedziale [0,2] z h = 0.1 oraz h = 0.01.\n\n")
pause;

range = [0 2];
N = [20 200];

f{1} = {@(x) sin(x)+cos(x), @(x) -1, @(x) 1};
f{2} = {@(x) exp(-x),@(x) 0, @(x) 0, @(x) -3};
f{3} = {@(x) sin(x),@(x) 0,@(x) 1, @(x) 0, @(x) 5};
f{4} = {@(x) sin(x),@(x) 0,@(x) 0,@(x) 1,@(x) 4,@(x) 1};
f{5} = {@(x) sin(x).*exp(-x),@(x) 0,@(x) -2,@(x) 0, @(x) 1,@(x) 0, @(x) 1};

Fsolve{1} = @(x) exp(x) - cos(x);
Fsolve{2} = @(x) -1./(3*exp(x))+ 2.*x./3 + 1/3;
Fsolve{3} = @(x) sqrt(5)*sin(x/sqrt(5)) - 25*cos(x./sqrt(5))/4 ...
    + cos(x)/4 + 6;
Fsolve{4} = @(x) cos(x)./4 + ((26*sqrt(3)+45)/24)*exp(sqrt(3)*x-2*x) - ...
    ((26*sqrt(3)-45)/24)*exp(-sqrt(3)*x-2*x)+x-4;
Fsolve{5} = @(x) sin(sqrt(2)*x)./(12*sqrt(2)) - cos(sqrt(2)*x)./6 + ...
    sin(x)./(10*exp(x)) + cos(x)./(20*exp(x)) - (2/15)*exp(x) + 5/4;

Y0{1} = 0;
Y0{2} = [0 1];
Y0{3} = [0 1 1];
Y0{4} = [0 0 0 0];
Y0{5} = [1 0 0 0 -1];

eqn{1} = "y' - y = sin(x) + cos(x)";
eqn{2} = "-3y'' = exp(-x)";
eqn{3} = "5*y^(3) + y' = sin(x)";
eqn{4} = "y^(4) + 4*y^(3) + y'' = sin(x)";
eqn{5} = "y^(5) + y^(3) - 2*y' = sin(x)*exp(-x)";

H = zeros(1,2);

for i = 1:5
    fprintf("%d. Równanie postaci:\n",i);
    disp(eqn{i});
    [X1,Y1,h1] = P1Z24_LSZ_runge_kutta(f{i},range,Y0{i},N(1));
    [~,Y2] = runge_kutta_classic(modify_f(f{i}),range,Y0{i},N(1));
    Y3 = ABMPC(f{i},[0 Y0{i}]',range,N(1));
    Yog = Fsolve{i}(X1)';
    err1(1) = max(abs(Yog - Y1));
    err2(1) = max(abs(Yog - Y2));
    err3(1) = max(abs(Yog - Y3'));
    H(1) = h1;

    [X1,Y1,h2] = P1Z24_LSZ_runge_kutta(f{i},range,Y0{i},N(2));
    [~,Y2] = runge_kutta_classic(modify_f(f{i}),range,Y0{i},N(2));
    Y3 = ABMPC(f{i},[0 Y0{i}]',range,N(2));
    Yog = Fsolve{i}(X1)';
    err1(2) = max(abs(Yog - Y1));
    err2(2) = max(abs(Yog - Y2));
    err3(2) = max(abs(Yog - Y3'));
    H(2) = h2;
    tab = table;
    tab.h = H';
    tab.Runge_Kutta_3_8 = err1';
    tab.Runge_Kutta_klasyczny = err2';
    tab.Adams_Bashforth_Moulton = err3';
    disp(tab);
    fprintf("\n");
    pause;
end

end

