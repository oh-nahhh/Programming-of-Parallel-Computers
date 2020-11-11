%Speed up plots

clear all; clc;
threads  = [1,2,4,8,11,14];


results_800  = [0.000328 0.000378 0.003892 0.050652 0.119161 0.161968];
results_8000  = [0.001730 0.001385 0.005117 0.078531 0.331402 0.885245];
speedup_800 = results_800(1) ./ results_800;
speedup_8000 = results_8000(1) ./ results_8000;
results_10000 = [0.003069 0.002083 0.004003 0.075439 0.376232  1.361676 ];
speedup_10000 = results_10000(1) ./ results_10000;
results_100000 = [0.020062, 0.017962 , 0.015491 , 0.092553, 0.543663, 2.464194];
speedup_100000 = results_100000(1) ./ results_100000;
results_1000000 = [0.22385, 0.160895 , 0.117940 , 0.138814, 0.388375, 0.735183];
speedup_1000000 = results_1000000(1) ./ results_1000000;
results_10000000 = [4.143608, 3.644617 , 2.054690, 0.873692, 1.153446, 2.704628 ];
speedup_10000000 = results_10000000(1) ./ results_10000000;

results_hej = [30.413834, 21.872313 , 16.373299, 10.330113, 9.250477, 15.913822];
speedup_100000000 = results_hej(1) ./ results_hej;


plot(threads, speedup_800, threads, speedup_8000,threads, speedup_10000,threads, speedup_100000,threads,speedup_1000000,threads,speedup_10000000, threads, speedup_100000000, 'LineWidth', 1.5);


set(gca, 'FontSize', 12)
title('Divide and Conquer parallelization')
legend('800','8,000', '10,000','100,000', '1,000,000', '10,000,000', '100,000,000')
ylabel('Speed-up')
xlabel('Nbr of threads')
