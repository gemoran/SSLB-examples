Y = readmatrix('data/Y_quantile.txt');
    
tic
res = SSBiEM(Y, 50);
t = toc;

writematrix(res.V, "matlab/X_norm.txt");
writematrix(res.h, "matlab/X_norm_support.txt");

writematrix(transpose(res.Z), "matlab/B_norm.txt");
writematrix(transpose(res.g), "matlab/B_norm_support.txt");
