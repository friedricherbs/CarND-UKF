psi = atan2(vy_gt, vx_gt);
psi_d = diff(psi);
for i = 1:length(psi_d)
    curr_psi = psi_d(i);
    while curr_psi >= pi-0.001
        curr_psi = curr_psi - 2*pi;
    end
    while curr_psi <= -pi+0.001
        curr_psi = curr_psi + 2*pi;
    end
    psi_d(i) = curr_psi;
end