function [X,Y,h] = runge_kutta_classic(F,range,y0,N)
% Projekt 1, zadanie 24
% Łukasz Szymczyk, 320744
%
% Funkcja rozwiązuje równanie liniowe dowolnego rzędu metodą 
% Rungego - Kutty rzędu 4 określonej wzorem klasycznym.
%   F     - funkcja anonimowa postaci @(x,y) = [y(2);y(3);...;y(m);f(x,y)], 
%           gdzie y to wektor.
%   range - dwuelementowy wektor [a, b] określający przedział, na którym
%           przybliżamy szukaną funkcję
%   Y0    - m-elementowy wektor zawierający warunki początkowe równania
%           [y0, y'0, ... , ]
%   N     - liczba kroków
% Wyjście:
%   X     - wektor (rozmiar N+1x1) równoodległych węzłów, w których 
%           przybliżana jest wartość szukanej funkcji
%   Y     - wektor (rozmiar N+1x1) zawierający wyznacznone przybliżenia 
%           wartości szukanej funkcji w węzłach zawartych w wektorze X
%   h     - krok całkowania

% Zainicjowanie potrzebnych zmiennych
a = range(1); 
b = range(2);
h = (b-a)/N;
X = a:h:b;
h = X(2) - X(1);
d = length(y0);
N = size(X,2);
Y = zeros(N,d);
Y(1,:) = y0;

% Obliczanie wartości Y w kolejnych węzłach z X
for i = 1:N-1
   k1 = h*F(X(i),Y(i,:))'; 
   k2 = h*F(X(i) + h./2,Y(i,:) + k1./2)';
   k3 = h*F(X(i) + h/2,Y(i,:) + k2/2)';
   k4 = h*F(X(i) + h,Y(i,:) + k3)';
   
   Y(i+1,:) = Y(i,:)+(k1 + 2*k2 + 2*k3 + k4)./6;
end

% Wyłuskanie 1 kolumny, gdyż w niej znajduje sie rozwiązanie równania
Y = Y(:,1);

end % function
