%Vikram Vijayakumar (02068559)
%MTH 565 Project 5
%https://github.com/msssm/lecture_files/blob/master/simulation_networks/small_world.m
%https://www.mathworks.com/matlabcentral/fileexchange/29786-random-regular-generator
%I had reviewed the above github link to obtain examples for clustering
%coefficients

n = 8000;               %Number of vetices, we will change this for n = 500 and 8000
k = 30;                
target_diameter = 5;   %Desired diameter of the small-world graph, changed it to 5 for n = 500 to enable faster results

%Create a regular lattice
G = zeros(n); %Initialize an empty adjacency matrix

%Connect each node to k/2 neighbors on either side
for i = 1:n
    for j = 1:k/2
        neighbor = mod(i + j - 1, n) + 1;
        G(i, neighbor) = 1;
        G(neighbor, i) = 1;
    end
end

%Initial measurements
initial_diameter = calculate_diameter(G);
initial_clustering = calculate_clustering_coefficient(G);
disp(['Initial Diameter = ', num2str(initial_diameter)]);
disp(['Initial Clustering Coefficient = ', num2str(initial_clustering)]);

%plot the initial graph
y = graph(G);
figure;
initial_plot = plot(y, 'Layout', 'circle');
title(['Lattice with n=', num2str(n), ', k=', num2str(k)]);

%Rewire edges and measure changes
q_values = [];    %Array to store percentage of rewired edges
diameters = [];   %Array to store diameter values
clusterings = []; %Array to store clustering coefficient values

total_edges = n * k / 2; %No of edges in the regular lattice
rewired_edges = 0;       %Counter for rewired edges

while calculate_diameter(G) > target_diameter
    %Randomly rewire an edge
    [G, rewired_edges] = rewire_edge(G, n, rewired_edges);

    q = (rewired_edges / total_edges) * 100;
    if floor(q) > length(q_values) %Calculate full percentage changes
        current_diameter = calculate_diameter(G);
        current_clustering = calculate_clustering_coefficient(G);

        %Plot d(q) and C(q) 
        q_values = [q_values, q];
        diameters = [diameters, current_diameter];
        clusterings = [clusterings, current_clustering];
        
        fprintf('q: %.1f%%, Diameter: %f, Clustering Coefficient: %f\n', q, current_diameter, current_clustering);
    end
end

figure;
plot(q_values, diameters / initial_diameter, '-o');
hold on;
plot(q_values, clusterings / initial_clustering, '-x');
xlabel('Percentage of Rewired Edges (q)');
ylabel('Normalized Diameter and Clustering Coefficient');
legend('d(q) / d(0)', 'C(q) / C(0)');
title('Small-World Network Simulation');
hold off;


%Function to Choose a random edge to rewire
function [G, rewired_edges] = rewire_edge(G, n, rewired_edges)
[row, col] = find(triu(G)); 
idx = randi(length(row));
u = row(idx);
v = col(idx);
%Remove the chosen edge
G(u, v) = 0;
G(v, u) = 0;
rewired_edges = rewired_edges + 1;

%Add a new random edge
while true
    new_u = randi(n);
    new_v = randi(n);
    if new_u ~= new_v && G(new_u, new_v) == 0
        G(new_u, new_v) = 1;
        G(new_v, new_u) = 1;
        break;
    end
end
end

%Fucntion to calculate diameter
function d = calculate_diameter(G)
graph_obj = graph(G);

%Calculate shortest path distances between all pairs
D = distances(graph_obj);

%Set infinite distances to zero to exclude disconnected pairs
D(D == Inf) = 0;

%Find the diameter as the longest shortest path
d = max(D(:));
end

%Function to calaculate lustering coefficient
function C = calculate_clustering_coefficient(G)
n = size(G, 1);
local_cluster = zeros(1, n);

for i = 1:n
    neighbors = find(G(i, :));
    k = length(neighbors);
    if k > 1
        subG = G(neighbors, neighbors);
        num_links = sum(sum(subG)) / 2;
        local_cluster(i) = (2 * num_links) / (k * (k - 1));
    else
        local_cluster(i) = 0;
    end
end
C = mean(local_cluster);
end