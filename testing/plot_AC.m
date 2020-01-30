% plot functional connectivity matrix and coherence
% INPUT
% coh - coherence computed as the ratio between maximum weighted (with eigenvector) eigenvalue and sum of all weighted eigenvalues
% C   - functional connectivity matrix (channel X channel)

function plot_AC(coh,C)
    
    fig_C = figure;
    
    figure(fig_C)
    C(1:size(C,1)+1:end) = 0; 
    imagesc(C)
    title(sprintf('coh:%.4f',coh))
    