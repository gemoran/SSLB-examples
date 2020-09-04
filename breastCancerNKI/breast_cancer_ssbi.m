Y = readmatrix('data/Y_raw.txt');
    
tic
res = SSBiEM(Y, 30);
t = toc;

writematrix(res.V, "matlab/X.txt");
writematrix(res.h, "matlab/X_support.txt");

writematrix(transpose(res.Z), "matlab/B.txt");
writematrix(transpose(res.g), "matlab/B_support.txt");
