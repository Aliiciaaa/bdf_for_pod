% This code is for completing Table 4 and 5 in the article.
clear all;
close all;

nombre_base2 = 'results_POD_';
nr = 18;

num_metodos = 5;
Mss = [64,128,256,512,1024];
nM = length(Mss);

errL2 = zeros(num_metodos, nM);
errH1 = zeros(num_metodos, nM);

nm = zeros(num_metodos, nM);
for j = 1:nM
    M = Mss(j);
    nombre2 = [nombre_base2, num2str(M),'_r',num2str(nr),'.mat'];
    filePath = fullfile('..','01_compute_pod_approx', nombre2);
    load(filePath);
   
    errL2(2,j) = max(errL2CombinedBdf2(2));
    errL2(3,j) = max(errL2CombinedBdf3(2:3));
    errL2(4,j) = max(errL2CombinedBdf4(2:4));
    errL2(5,j) = max(errL2CombinedBdf5(2:5));
    
     
    errH1(2,j) = max(errH1CombinedBdf2(2));
    errH1(3,j) = max(errH1CombinedBdf3(2:3));
    errH1(4,j) = max(errH1CombinedBdf4(2:4));
    errH1(5,j) = max(errH1CombinedBdf5(2:5));

    %nm(:, j) = nm_iter;

end