U1pp = 0;
U2pp = 0;
U3pp = 0;
U4pp = 0;
Y1pp = 0;
Y2pp = 0;
Y3pp = 0;
Tp = 0.5;

kk = 300;
t = (1:kk)';
U1 = U1pp*ones(kk, 1);
U2 = U2pp*ones(kk, 1);
U3 = U3pp*ones(kk, 1);
U4 = 1*ones(kk, 1);
Y1 = Y1pp*ones(kk, 1);
Y2 = Y2pp*ones(kk, 1);
Y3 = Y3pp*ones(kk, 1);

for k=7:kk
    [Y1(k),Y2(k),Y3(k)] = symulacja_obiektu2y_p4(U1(k-1),U1(k-2),U1(k-3),U1(k-4),...
    U2(k-1),U2(k-2),U2(k-3),U2(k-4),...
    U3(k-1),U3(k-2),U3(k-3),U3(k-4),...
    U4(k-1),U4(k-2),U4(k-3),U4(k-4),...
	Y1(k-1),Y1(k-2),Y1(k-3),Y1(k-4),...
	Y2(k-1),Y2(k-2),Y2(k-3),Y2(k-4),...
	Y3(k-1),Y3(k-2),Y3(k-3),Y3(k-4));
end

% plot_data = [t,Y1];
% save("zad1y.txt","plot_data","-ascii")
%wykresy
figure
plot(Y1, '--'); hold on; plot(Y2, '-'); hold on; plot(Y3, '-.', 'Color','magenta');
ylim([0 3]); % Ustawia zakres osi Y
legend("Y1", 'Y2', 'Y3')
title("Symulacja obiektu");
xlabel('Numer próbki');
ylabel("Wartość wyjścia");
matlab2tikz('zad2_u4.tex' , 'showInfo' , false)