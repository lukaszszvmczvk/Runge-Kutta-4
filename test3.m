function [] = test3()

fprintf("Test sprawdza wartość błędu globalnego w zależności od kroku\n" + ...
    "całkowania oraz przedstawia stosunek tych błędów.\n" + ...
    "Test został przeprowadzony na 3 równaniach różnych " + ...
    "rzędów\nrozwiązanych na przedziale [0,1]. Testowane są h = 0.1,\n" + ...
    "h = 0.01 oraz h = 0.001.\n\n");
pause;

N1 = [10 100 1000];

f{1} = {@(x) x*sin(x), @(x) 1, @(x) 0, @(x) 1};
f{2} = {@(x) 0, @(x) 0, @(x) 0, @(x) 1, @(x) 2, @(x) 1};
f{3} = {@(x) cos(x), @(x) cos(x), @(x) 1};

F_sol{1} = @(x) x.*sin(x)/4 -x.^2.*cos(x)/4;
F_sol{2} = @(x) (x+3)./exp(x)+4*x-1;
F_sol{3} = @(x) 1-2./exp(sin(x));

range = [0 1];

Y0{1} = [0 0];
Y0{2} = [2 2 1 0];
Y0{3} = (-1);

errors1 = zeros(3,3);

H1 = zeros(3,3);

ratio1 = ones(3,3);

for i = 1:3
    for j = 1:3
        [X,Y,h] = P1Z24_LSZ_runge_kutta(f{i},range,Y0{i},N1(j));
        Yog = F_sol{i}(X);
        errors1(i,j) = max(abs(Yog - Y'));
        H1(i,j) = h;
    end
end

for i = 2:3
    ratio1(:,i) = errors1(:,i)./errors1(:,1);
end

eqn{1} = "y'' + y = x*sin(x)";
eqn{2} = "y^(4) + 2*y^(3) + y'' = 0";
eqn{3} = "y' + cos(x)*y = cos(x)";


eqn_sol{1} = "y = (x*sin(x) - x^2*cos(x))/4";
eqn_sol{2} = "y = (x + 3)/exp(x) + 4x - 1";
eqn_sol{3} = "y = 1 - 2/exp(sin(x))";

for i = 1:3
    y0 = "[";
    for j = 1:length(Y0{i})
        if j == length(Y0{i})
            y0 = y0 + int2str(Y0{i}(j));
        else
            y0 = y0 + int2str(Y0{i}(j))+", ";
        end
    end
    y0 = y0 + "]";
    fprintf("%d. Równanie różniczkowe postaci:\n" + eqn{i} + "\n",i);
    fprintf("Y0: " + y0 + "\n");
    fprintf("Rozwiązanie analityczne: " + eqn_sol{i} + "\n");
    tab = table;
    tab.h = H1(i,:)';
    tab.BladGlobalny = errors1(i,:)';
    disp(tab);
    fprintf("Stosunek h1 : h2 : h3 =  1 : 0.1 : 0.01\n");
    fprintf("Stosunek err1 : err2 : err3 = %d : %d : %d\n", ratio1(i,1), ...
        ratio1(i,2), ratio1(i,3));
    pause;
    fprintf("\n");
end

end

