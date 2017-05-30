clear 
str1 = {'0.5', '1', '2', '4', '8'};
str2 = {'Duhamel', 'Polynomial', 'Eigenvalues'};
colors = {'b-', 'r--', 'k-.'};
close all
figure
hold on
for i = 1:length(str2)
    for j = 1:length(str1)
        filename = [str2{i}, '_', str1{j}, 'C.mat'];
        load(filename)
        cap = [0, cumsum(diff(ts)*cur)/3600];
        h(i) = plot(cap, v, colors{i});
    end
end

legend(h, str2)

xlabel('Specific capacity /A m^{-2}')
ylabel('Voltage /V')