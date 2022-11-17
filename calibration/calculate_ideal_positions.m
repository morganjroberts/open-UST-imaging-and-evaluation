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
    i_end = i_start + N_per_mod - 1;
    element_positions(i_start:i_end, :) = [mod_x_new, mod_y_new];
end

% calculate distances between pairs
distance_matrix = squareform(pdist(element_positions));

% calculate angle between tx/rx vector and normal axis of tx element
angle_to_norm_matrix = zeros(Ntd, Ntd);
for idx = 1:Ntd
    for jdx = 1:Ntd
        if idx == jdx
            angle_to_norm_matrix(idx, jdx) = NaN;
        else
            normal_angle  = floor((idx - 1)/N_per_mod) * delta_angle; % relative to positive y axis
            normal_vector = [sind(normal_angle), cosd(normal_angle), 0];
            pair_vector   = [element_positions(jdx,:) - element_positions(idx,:), 0];
            CosTheta      = dot(normal_vector,pair_vector) / (norm(normal_vector)*norm(pair_vector));
            angle_to_norm_matrix(idx, jdx) = real(acosd(-CosTheta));
        end    
    end
end

%% plot and save data

figure;
subplot(1, 2, 1);
hold on
plot(element_positions(:,1), element_positions(:,2), 'k.')
axis image
xlabel('x-position [m]');
ylabel('y-position [m]');

subplot(1, 2, 2);
imagesc(angle_to_norm_matrix);
C = colorbar;
xlabel('Receive Channel');
ylabel('Transmit Channel');
ylabel(C, 'Receiver Angle to Transmitter Normal [deg]');
axis image

set(gcf, 'position', [581, 576, 1180, 423]);


save('256_element_UST_array_predicted_positions.mat', ...
    'element_positions', 'distance_matrix', 'angle_to_norm_matrix');