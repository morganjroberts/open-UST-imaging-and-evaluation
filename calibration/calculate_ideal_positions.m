% script to calculate the ideal positions of the transducer elements

% close all
clearvars

r         = 0.11;    % array radius in m
p         = 2.54e-3; % element pitch in m
N_per_mod = 16; % number of elements per array module
Nmod      = 16; % number of modules
Ntd       = N_per_mod * Nmod; % total number of elements
delta_angle = 360 / Nmod;

element_positions = zeros(Ntd, 2); % initialise position vector

%% calculate positions of a single module

mod_len = ((N_per_mod) - 1) * p;
mod_x   = linspace(-mod_len/2, mod_len/2, N_per_mod)';
mod_y   = r * ones(N_per_mod, 1);

%% calculate positions of all modules

% loop over modules
for mdx = 1:Nmod
    angle = -(mdx - 1) * delta_angle;
    s = sind(angle);
    c = cosd(angle);
    
    % rotate module coordinates
    mod_x_new = mod_x * c - mod_y * s;
    mod_y_new = mod_x * s + mod_y * c;
   
    % save positions
    i_start = ((mdx - 1) * N_per_mod) + 1;
    i_end   = i_start + N_per_mod - 1;
    element_positions(i_start:i_end, :) = [mod_x_new, mod_y_new];
end

% calculate distances between pairs
distance_matrix = squareform(pdist(element_positions));

% calculate angle between tx/rx vector and normal axis of tx element
angle_TxRx_to_TxNorm = zeros(Ntd, Ntd);
for idx = 1:Ntd
    for jdx = 1:Ntd
        if idx == jdx
            angle_TxRx_to_TxNorm(idx, jdx) = NaN;
        else
            normal_angle  = floor((idx - 1)/N_per_mod) * delta_angle; % relative to positive y axis
            normal_vector = [sind(normal_angle), cosd(normal_angle), 0];
            pair_vector   = [element_positions(jdx,:) - element_positions(idx,:), 0];
            CosTheta      = dot(normal_vector,pair_vector) / (norm(normal_vector)*norm(pair_vector));
            angle_TxRx_to_TxNorm(idx, jdx) = real(acosd(-CosTheta));
        end    
    end
end

% calculate angle between rx/tx vector and normal axis of rx element
angle_RxTx_to_RxNorm = angle_TxRx_to_TxNorm';

% calculate angle between tx normal (beam axis) and rx normal (beam axis)
angle_TxNorm_to_RxNorm = zeros(Ntd, Ntd);
for idx = 1:Ntd
    for jdx = 1:Ntd
        if idx == jdx
            angle_TxNorm_to_RxNorm(idx, jdx) = NaN;
        else
            tx_normal_angle  = floor((idx - 1)/N_per_mod) * delta_angle; % relative to positive y axis
            tx_normal_vector = [sind(tx_normal_angle), cosd(tx_normal_angle), 0];

            rx_normal_angle  = floor((jdx - 1)/N_per_mod) * delta_angle; % relative to positive y axis
            rx_normal_vector = [sind(rx_normal_angle), cosd(rx_normal_angle), 0];            

            CosTheta      = dot(tx_normal_vector,rx_normal_vector) / (norm(tx_normal_vector)*norm(rx_normal_vector));
            angle_TxNorm_to_RxNorm(idx, jdx) = real(acosd(-CosTheta));
        end    
    end
end

%% Make mask for whether a receiver is directly opposite the transmitter


is_opposite = angle_TxRx_to_TxNorm;
for tdx = 1:Ntd
    is_opposite(tdx, is_opposite(tdx,:) ~= min(is_opposite(tdx,:)) ) = NaN;
end
is_opposite = ~isnan(is_opposite);

%% plot and save data

figure;
subplot(1, 3, 1);
hold on
plot(element_positions(:,1), element_positions(:,2), 'k.')
axis image
xlabel('x-position [m]');
ylabel('y-position [m]');

subplot(1, 3, 2);
imagesc(angle_TxRx_to_TxNorm);
C = colorbar;
xlabel('Receive Channel');
ylabel('Transmit Channel');
ylabel(C, 'Tx-Rx-Vector Angle to Tx-Normal [deg]');
axis image

subplot(1, 3, 3);
imagesc(angle_TxNorm_to_RxNorm);
C = colorbar;
xlabel('Receive Channel');
ylabel('Transmit Channel');
ylabel(C, 'Rx-Normal Angle to Tx-Normal [deg]');
axis image

set(gcf, 'position', [43 463 1725 378]);


save('ideal_element_positions.mat', ...
    'element_positions', 'distance_matrix', 'angle_TxRx_to_TxNorm', 'angle_RxTx_to_RxNorm', 'angle_TxNorm_to_RxNorm', 'is_opposite');