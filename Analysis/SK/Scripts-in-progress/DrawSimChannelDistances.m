% Generate log-normal distribution of IPIs based on MEC-4 spacing from
% Cueva et al., 2007
peak = 0.2;
width = 2.34;
rando = rng(40); %set same seed for reproducible realizations: 19, 26, 40 look nice
rand_IPI = lognrnd(peak,width,500,1);
% r_sum = cumsum(rand_IPI);
% histogram(rand_IPI,'Normalization','probability','BinWidth',0.2,'BinLimits',[0,15])

% Include only IPIs within the range plotted in Cueva et al. (less than
% 15um).
IPI_cutoff = 15;
cut_IPI = rand_IPI(rand_IPI < IPI_cutoff);

% Sum IPIs cumulatively to simulate puncta positions based on log-normal
% distribution. Include only positions up to the length of ALM (500um,
% about half the length of the body).
puncta_loc = cumsum(cut_IPI);
% sum(puncta_loc<500)
puncta_loc = puncta_loc(puncta_loc<500);

% Designate conditions for plotting a channel as open (filled circle) or
% closed (x).
% For simplicity, we'll just use the final distance at which channels are
% open for each displacement (i.e., open plus subconducting vs. closed).
% 5um stimulus: ~50um
% 10um stimulus: ~85um
dist_5 = 50;
dist_10 = 85;
open_5 = -puncta_loc( abs(puncta_loc) <= dist_5 );
closed_5 = -puncta_loc( abs(puncta_loc) > dist_5 );
open_5_y = ones(length(open_5),1);
closed_5_y = ones(length(closed_5),1);

open_10 = -puncta_loc( abs(puncta_loc) <= dist_10 );
closed_10 = -puncta_loc( abs(puncta_loc) > dist_10 );
open_10_y = ones(length(open_10),1)*0.8;
closed_10_y = ones(length(closed_10),1)*0.8;

figure();hold on;
scatter(open_5,open_5_y,'MarkerFaceColor','b','MarkerEdgeColor','k')
scatter(closed_5,closed_5_y,'kx')
scatter(open_10,open_10_y,'MarkerFaceColor','m','MarkerEdgeColor','k')
scatter(closed_10,closed_10_y,'kx')
ylim([-1 3])
xlim([-150 0])