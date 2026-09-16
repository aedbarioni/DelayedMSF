% load("datafile.mat")
figure();
imagesc(mtle(end:-1:1,1:1:end)) 
xlabel('\alpha')
ylabel('\beta')
title('Lyapunov Exponent')
colormap(bluewhitered),colorbar
xticks(1:(m/(2*re_fin)):m); xticklabels(-re_fin:10:re_fin);
yticks(1:(m/(4*im_fin)):m); yticklabels(-im_fin:5:re_fin);

hold on
contour_tres = 1e-2; % This must change depending on the number of rows and columns in the matrix mtle (in this case, the matrix is mxm)
contour(abs(mtle(end:-1:1,1:end)), [contour_tres contour_tres], 'k', 'LineWidth', 1.5);
hold off