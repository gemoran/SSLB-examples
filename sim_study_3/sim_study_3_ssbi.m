
nrep = 50;
t = zeros(nrep, 1);

for r = 1:nrep
    Y = readmatrix(sprintf('data/rep%d/Y.txt',r));
    Y = transpose(Y);
    
    tic
    res = SSBiEM(Y, 15);
    t(r,1) = toc;

    writematrix(res.V, sprintf("data/rep%d/matlab/X.txt", r));
    writematrix(res.h, sprintf("data/rep%d/matlab/X_support.txt", r));

    writematrix(transpose(res.Z), sprintf("data/rep%d/matlab/B.txt", r));
    writematrix(transpose(res.g), sprintf("data/rep%d/matlab/B_support.txt", r));

end

writematrix(t, "results/SSBiEM_timing.txt")

 r = 1;
 Y = readmatrix(sprintf('data/rep%d/Y.txt',r));
 Y = transpose(Y);
    
 tic
 res = SSBiEM(Y, 30);
 t = toc;

 writematrix(res.V, sprintf("data/rep%d/matlab/X_K_large.txt", r));
 writematrix(res.h, sprintf("data/rep%d/matlab/X_support_K_large.txt", r));

 writematrix(transpose(res.Z), sprintf("data/rep%d/matlab/B_K_large.txt", r));
 writematrix(transpose(res.g), sprintf("data/rep%d/matlab/B_support_K_large.txt", r));
