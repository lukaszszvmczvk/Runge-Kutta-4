function [] = test_num_1()

fprintf("Test sprawdza od jakiej wartości h błąd globalny przestaje maleć, " + ...
    "a zaczyna\n" + ...
    "wzrastać. Test został przeprowadzony na dwóch równaniach, jedynym\n" + ...
    "2 rzędu oraz drugim 3 rzędu. Rozwiązania rówanań są przybliżane na\n" + ...
    "przedziale [1,2]. Początkowa wartość h wynosi 0.1, a następnie\n" + ...
    "jest 10-krotnie zmniejszana w każdej iteracji testu.\n" + ...
    "Generowanie rozwiązań trwa kilka sekund, ponieważ dla bardzo małego h" + ...
    "\nilość wykonywanych operacji jest bardzo duża.\n\n");
pause;

range = [1 2];

f{1} = {@(x) x^2+x, @(x) 1, @(x) 0, @(x) 1};
f{2} = {@(x) -4*x^2, @(x) 0, @(x) 0,@(x) 0, @(x) 1};

F_sol{1} = @(x) -3*cos(1)*sin(x)+3*sin(1)*cos(x) + x.^2 + x-2;
F_sol{2} = @(x) -x.^5/15 + 7*x.^2/6 - 2*x + 19/10;

Y0{1} = [0 0];
Y0{2} = [1 0 1];

i = 1;
N = 10;
error = zeros(1,100);
H = zeros(1,100);
fprintf("1. Równanie rzędu 2 określone wzorem:\ny'' + y = x^2 + x\n");
pause;
while (true)
    [X,Y,h] = P1Z24_LSZ_runge_kutta(f{1},range,Y0{1},N);
    Yog = F_sol{1}(X);
    error(i) = max(abs(Y - Yog'));
    H(i) = h;
    N = N*10;
    if (i>=2)
        if(error(i-1)<error(i))
            break;
        end
    end
    i = i+1;
end
tab = table;
tab.h = H(1:i)';
tab.Blad_Globalny = error(1:i)';
disp(tab);
pause;

error = zeros(1,100);
H = zeros(1,100);
i = 1;
N = 10;
fprintf("2. Równanie rzędu 3 określone wzorem:\ny^(3) = -4*x^2\n");
pause;
while (true)
    [X,Y,h] = P1Z24_LSZ_runge_kutta(f{2},range,Y0{2},N);
    Yog = F_sol{2}(X);
    error(i) = max(abs(Y - Yog'));
    H(i) = h;
    N = N*10;
    if (i>=2)
        if(error(i-1)<error(i))
            break;
        end
    end
    i = i+1;
end
tab = table;
tab.h = H(1:i)';
tab.Blad_Globalny = error(1:i)';
disp(tab);
pause;
end

