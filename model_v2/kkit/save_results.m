function save_results(t, y, path, with_pump)
f_names={'ca_matlab';'raf_matlab';'araf_matlab';
    'mapk_matalb';'pmapk_matlab';'phsph_matlab';
    'ka_matlab';'pka_matlab';'pkc_matlab';
    'apkc_matlab'; 'aa_matlab'};
n_names = length(f_names);
y_plot = [y(:,1:4) y(:,6:7) y(:,9) y(:,11:14)];
if with_pump
    f_names{n_names+1}='pmca';
    y_plot=[y_plot y(:,17)];
end

for i=1 : length(f_names)
    dlmwrite(strcat(path,'/',f_names{i}),[t,y_plot(:,i)]);
end


end