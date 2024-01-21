n_in = 4;
n_out = 3;
kk = 300;
Y_all = [];
Y_s = {};
for n = 1:n_in
    Y1 = zeros(kk, 1);
    Y2 = zeros(kk, 1);
    Y3 = zeros(kk, 1);
    Uz = zeros(kk, 1);
    Uo = ones(kk, 1);
    U = [Uz Uz Uz Uz];
    U(:,n) = Uo;
    for k=7:kk
        [Y1(k),Y2(k),Y3(k)] = symulacja_obiektu2y_p4(U(k-1,1),U(k-2,1),U(k-3,1),U(k-4,1),...
        U(k-1,2),U(k-2,2),U(k-3,2),U(k-4,2),...
        U(k-1,3),U(k-2,3),U(k-3,3),U(k-4,3),...
        U(k-1,4),U(k-2,4),U(k-3,4),U(k-4,4),...
	    Y1(k-1),Y1(k-2),Y1(k-3),Y1(k-4),...
	    Y2(k-1),Y2(k-2),Y2(k-3),Y2(k-4),...
	    Y3(k-1),Y3(k-2),Y3(k-3),Y3(k-4));
    end
    Y_all = [Y_all Y1 Y2 Y3]

end

for m = 1:length(Y1)
    Y_s{m} =[Y_all(m,1) Y_all(m,4) Y_all(m,7) Y_all(m,10); Y_all(m,2) Y_all(m,5) Y_all(m,8) Y_all(m,11); Y_all(m,3) Y_all(m,6) Y_all(m,9) Y_all(m,12)];
end

Y_s = Y_s(7:end)

