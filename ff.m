function S = get_s(Upp, Ypp, sim_end)
    ny = 3;
    nu = 4;
    S = {};
    for u = 1:nu
        Y_step = [ones(1, sim_end + 1)*Ypp(1); ones(1, sim_end + 1) * Ypp(2); ones(1, sim_end + 1) * Ypp(3)];
        U = [ones(1, sim_end+1)*Upp(1); ones(1, sim_end+1)*Upp(2); ones(1, sim_end+1)*Upp(3); ones(1, sim_end+1)*Upp(4)];
        U(u, :) = ones(1, sim_end+1);
        U(u, 1) = Upp(u);
        for i = 1:sim_end+1
            [Y_step(1, i), Y_step(2, i),Y_step(3, i)] = symulacja_obiektu4y_p4( ...
                U(1, (max(1, i-1))), U(1, (max(1, i-2))), U(1, (max(1, i-3))), U(1, (max(1, i-4))), ...
                U(2, (max(1, i-1))), U(2, (max(1, i-2))), U(2, (max(1, i-3))), U(2, (max(1, i-4))), ...
                U(3, (max(1, i-1))), U(3, (max(1, i-2))), U(3, (max(1, i-3))), U(3, (max(1, i-4))), ...
                U(4, (max(1, i-1))), U(4, (max(1, i-2))), U(4, (max(1, i-3))), U(4, (max(1, i-4))), ...
                Y_step(1, (max(1,i-1))), Y_step(1, (max(1, i-2))), Y_step(1, (max(1, i-3))), Y_step(1, max(1,i-4)), ...
                Y_step(2, (max(1,i-1))), Y_step(2, (max(1, i-2))), Y_step(2, (max(1, i-3))), Y_step(2, max(1,i-4)), ...
                Y_step(3, (max(1,i-1))), Y_step(3, (max(1, i-2))), Y_step(3, (max(1, i-3))), Y_step(3, max(1,i-4)));
           
            s_i_1 = Y_step(1, i);
            s_i_2 = Y_step(2, i);
            s_i_3 = Y_step(3, i);
            S{i}(1,u) = s_i_1;
            S{i}(2,u) = s_i_2;
            S{i}(3,u) = s_i_3;
        end
    end
end