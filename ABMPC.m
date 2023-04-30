function Y = ABMPC(Fx,Y0,xspan,n,usecor)
% Projekt 1, zadanie 24
% Łukasz Szymczyk, 320744
%
% Funkcja otrzymana od Mikołaja Wałachowskiego, wykorzystana w teście num.
%
% Rozwiązanie numeryczne liniowego równania różniczkowego dowolnego rzędu
% postaci a_m(x)*y^(m)+...+a_1(x)*y'+a_0(x)*y = b(x), przy użyciu metody
% predykator-korektor Adamsa-Bashfortha-Moultona 2-go rzędu
% Wejście:
%   Fx     - tablica komórek zawierających uchwyty do funkcji,
%          - w której F{1} = @(x) b(x) oraz F{i+2} = @(x) a_i(x)
%   Y0     - wektor zawierający kolejno przybliżenie początkowe x0 oraz
%          - oraz wartości y^(i) dla x0, przy 0 <= i <= m-1
%   xspan  - wektor postaci [a b], gdzie a i b to pierwszy i ostatni węzeł
%   n      - liczba iteracji, czyli n+1 ilość równoodległych węzłów,
%            w których przybliżamy funkcję
%   usecor - wartość logiczna determinująca, czy użyty zostanie korektor, 
%          - domyślnie true, czyli korektor włączony
% Wyjście:
%   Y      - wektor (rozmiar N+1x1) zawierający wyznacznone przybliżenia 
%           wartości szukanej funkcji


ord=length(Y0);
h=(xspan(2)-xspan(1))/n;
alf=[1.5,-0.5]; % wsp. A-B drugiego stopnia
mi=[0.5,0.5]; % wsp. A-M drugiego stopnia
Y=zeros(ord,n+1);
Y(:,1)=Y0;
% temp=zeros(ord,1);
% temp(1)=xspan(1)+h/2;
if nargin == 4
    usecor=true;
end
%stworzenie F takiego jak chcemy z Fx
     function out=Fvector(Y)
         out=Fx{1}(Y(1));
         for k=2:length(Y)
             out=out-Fx{k}(Y(1))*Y(k);
         end
         out=out/Fx{length(Y)+1}(Y(1));
     end
fcel=cell(ord,1);
fcel{1}=@(y) 1;
for i=2:ord-1
    fcel{i}=(@(y) y(i+1));
end
fcel{ord}=@Fvector;
F=@(y) cellfun(@(f)f(y),fcel);

%Zmodyfikowana metoda Eulera 
Y(:,2)=Y0+h*F(Y0+h/2*F(Y0));

for i=2:n
    Y0=Y(:,i-1);
    Y1=Y(:,i);
    % Predyktor 2 rzędu Adams-Bashforth
    Yp=Y1+h*(alf(1)*F(Y1)+alf(2)*F(Y0));
    % Korektor 2 rzędu Adams-Moulton
    if usecor
        Y(:,i+1)=Y1+h*(mi(1)*F(Yp)+mi(2)*F(Y1));
    else
        Y(:,i+1)=Yp;
    end
    
end
Y = Y(2,:);

end % function

