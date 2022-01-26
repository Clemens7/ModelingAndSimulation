
    tic;
    % Number of nodes in each direction
    numnodes_x = 240;
    numnodes_y = 80;

    % Number of timesteps over which to simulate.
    t_end = 20;

    % 3D array of nodes to store the vectors that represent the occupied
    % cells at each node.  
    %   0 - Cell unoccupied.
    %   1 - Cell occupied.
    %
    % 1st Index -- Node x-coordinate.
    % 2nd Index -- Node y-coordinate.
    % 3rd Index -- Cell number
    %
    % The elements of the occupancy vectors correspond to the cells in the 
    % following way:
    %
    %                2
    %                |  
    %            3 - O - 1
    %                |
    %                4
    nodes = zeros(numnodes_x, numnodes_y, 4);

    % Define the lattice velocities.
    c1 = [1; 0];
    c2 = [0; 1];
    c3 = [-1; 0];
    c4 = [0; -1];


    
    % Define a matrix that indicates where the flow obstacles are.
    %   0 - No obstacle present at that node.
    %   1 - Obstacle at the node.
    %
    % Don't forget to put 1's at the interior points, too!
    obstacle = zeros(numnodes_x, numnodes_y);

    % Insert a flat plate as the obstacle.
    for (j = 30:1:50)
        obstacle(80, j) = 1;
    end
        
    % Set up the simulation.
    for (i = 1:1:numnodes_x)
        for (j = 1:1:numnodes_y)  % Don't include the top and bottom walls.
            % Skip points on the obstacle boundary
            if (obstacle(i, j) ~= 1)               
                curr_cell = nodes(i, j, :);  % Get the cell for the current node.

                curr_cell(1) = 1;            % Put a particle in the cell flowing in the
                                             % rightward direction.
                                             
                nodes(i, j, :) = curr_cell;  % Reinsert the cell into the array.                             
            end
%               nodes(60,40,1) = 1;
%               nodes(70,40,3) = 1;
        end
    end
    

    % Carry out the simulation.
    for (t = 1:1:t_end)
        clf;
        % Carry out collisions at non-boundary nodes.
        for (i = 1:1:numnodes_x)
            for (j = 2:1:(numnodes_y - 1))  % Don't include the top and bottom walls.
                % Ensure that there's no obstacle in the way.
                if (obstacle(i, j) ~= 1)
                
                    % Extract the current cell.
                    cell_oc = nodes(i, j, :);

                    % Determine how many particles are in the cell.
                    numparts = sum(cell_oc);

                    % Determine and execute appropriate collision.
                    if (numparts == 2)
                        if ((cell_oc(1) == cell_oc(3)) && (cell_oc(2) == cell_oc(4)))
                            % Invert the cell contents.
                            nodes(i, j, :) = ~cell_oc;
                        else
                            nodes(i, j, :) = cell_oc;
                        end
                    else
                        nodes(i, j, :) = cell_oc;

                    end
                end
            end
        end
        
        % Carry out collisions along the top and bottom walls (no-slip).
        for (i = 1:1:numnodes_x)
            nodes(i, 1, :) = [nodes(i, 1, 3) nodes(i, 1, 4) nodes(i, 1, 1) nodes(i, 1, 2)];
            nodes(i, numnodes_y, :) = [nodes(i, numnodes_y, 3) nodes(i, numnodes_y, 4) nodes(i, numnodes_y, 1) nodes(i, numnodes_y, 2)];
        end
        
        % Carry out collisions at obstacle points (no-slip).
        for (i = 1:1:numnodes_x)
            for (j = 1:1:numnodes_y)
                if (obstacle(i, j) == 1)
                    nodes(i, j, :) = [nodes(i, j, 3) nodes(i, j, 4) nodes(i, j, 1) nodes(i, j, 2)];
                end
            end
        end
        
    
        % Create a new lattice which will hold the state of the current
        % lattice after propagation.
        n_nodes = zeros(numnodes_x, numnodes_y, 4);

        % Iterate over all the nodes, propagating the particles as we go.
        for (i = 1:1:numnodes_x)
            for(j = 1:1:numnodes_y)
                % Get the occupancy state of the current node.
                cell_oc = nodes(i, j, :);

                % Coordinates of the neighbor node.
                neighbor_x = 0;
                neighbor_y = 0;

                % Propagation in the 1-direction.
                neighbor_y = j;

                if (i == numnodes_x)
                    neighbor_x = numnodes_x;
                    n_cell_oc = n_nodes(neighbor_x, neighbor_y, :);
                    n_cell_oc(1) = 0;
                    n_nodes(neighbor_x, neighbor_y, :) = n_cell_oc; 
                else
                    neighbor_x = i + 1;
                    n_cell_oc = n_nodes(neighbor_x, neighbor_y, :);
                    n_cell_oc(1) = cell_oc(1);
                    n_nodes(neighbor_x, neighbor_y, :) = n_cell_oc;
                end

                %n_cell_oc = n_nodes(neighbor_x, neighbor_y, :);
                %n_cell_oc(1) = cell_oc(1);
                %n_nodes(neighbor_x, neighbor_y, :) = n_cell_oc;
                

                % Propagation in the 2-direction.
                neighbor_x = i;
                if (j ~= numnodes_y)
                    %neighbor_y = numnodes_y;
                    neighbor_y = j + 1;
                    n_cell_oc = n_nodes(neighbor_x, neighbor_y, :);
                    n_cell_oc(2) = cell_oc(2);
                    n_nodes(neighbor_x, neighbor_y, :) = n_cell_oc;
                    
                    
                %else
                    
                end



                

                % Propagation in the 3-direction.
                neighbor_y = j;

                if (i ~= 1)
                    %neighbor_x = numnodes_x; 
                    neighbor_x = i - 1;
                else
                    
                end

                n_cell_oc = n_nodes(neighbor_x, neighbor_y, :);
                n_cell_oc(3) = cell_oc(3);
                n_nodes(neighbor_x, neighbor_y, :) = n_cell_oc;
              

                % Propagation in the 4-direction.
                neighbor_x = i;

                if (j ~= 1)
                    %neighbor_y = numnodes_y;
                    neighbor_y = j - 1;
                    n_cell_oc = n_nodes(neighbor_x, neighbor_y, :);
                    n_cell_oc(4) = cell_oc(4);
                    n_nodes(neighbor_x, neighbor_y, :) = n_cell_oc;
                    
                    
                %else
                    
                end


                

            end
        end
        
        % Propagate the particles to their next nodes
        nodes = n_nodes;
        % create new particles
        for (j = 1:1:numnodes_y)
            nodes(1,j,1)= 1;

        end
        %nodes
        % Print the current time step every so often so we know that the 
        % program hasn't frozen or crashed.
        if (mod(t, t_end) == 0)
            disp(t);
        end
         
    
    % Subdivide the total domain into subdomains of size 32x32 for the
    % purposes of coarse-graining.  See pg.  51.
    grain_size = 1;
    grain_x = numnodes_x / grain_size;
    grain_y = numnodes_y / grain_size;
    
    % Pre-allocate vectors for the averaged velocities.
    av_vel_x_coords = zeros(1, grain_x * grain_y);
    av_vel_y_coords = zeros(1, grain_x * grain_y);
    av_vel_x_comps = zeros(1, grain_x * grain_y);
    av_vel_y_comps = zeros(1, grain_x * grain_y);

    %Pre-allocate vector for particle number

    par_num = zeros(grain_x, grain_y);
    average = zeros(grain_x,grain_y);
    streamfunc_u = zeros(grain_x,grain_y);
    streamfunc_v = zeros(grain_x,grain_y);
    % Iterate over the entire domain, averaging and storing the results as
    % we go.
    currval = 1;
    for (i = 1:1:grain_x)
            % Calculate the lower and upper x-boundaries.
            x_bd_l = (i - 1)*grain_size + 1;
            x_bd_u = i*grain_size;
        for (j = 1:1:grain_y)
            % Calculate the lower and upper y-boundaries.
            y_bd_l = (j - 1)*grain_size + 1;
            y_bd_u = j*grain_size;      
            
            % Get the number of particles moving in each direction in the
            % current subdomain.
            np = zeros(1, 4);
            np(1) = sum(sum(nodes(x_bd_l:1:x_bd_u, y_bd_l:1:y_bd_u, 1)));
            np(2) = sum(sum(nodes(x_bd_l:1:x_bd_u, y_bd_l:1:y_bd_u, 2)));
            np(3) = sum(sum(nodes(x_bd_l:1:x_bd_u, y_bd_l:1:y_bd_u, 3)));
            np(4) = sum(sum(nodes(x_bd_l:1:x_bd_u, y_bd_l:1:y_bd_u, 4)));

            par_num(i,j) = sum(np);

            
            % Compute the average velocity.
            av_vel = (1/(grain_size.^2))*(np(1)*c1 + np(2)*c2 + np(3)*c3 + np(4)*c4);
            
            % Store the velocity components.
            av_vel_x_comps(currval) = av_vel(1);
            av_vel_y_comps(currval) = av_vel(2);
            
            % Store the positional coordinates.
            av_vel_x_coords(currval) = i;
            av_vel_y_coords(currval) = j;
            average(i,j) = (av_vel(1)^2+av_vel(2)^2)^0.5;
            streamfunc_u(i,j) = av_vel(1);
            streamfunc_v(i,j) = av_vel(2);
            currval = currval + 1;
        end
    end
    [X,Y] = meshgrid(1:1:numnodes_x,1:1:numnodes_y);
    % Plot the average velocity field.
    contourf(transpose(par_num),[0  10 20  30],'LineStyle','none');
    colormap("cool")
    hold on;
    colorbar();
%     quiver(av_vel_x_coords, av_vel_y_coords, av_vel_x_comps, av_vel_y_comps);
%     hold on
    ys = 1:5:numnodes_x;
    xs = ones(1,length(ys));
    streamline(X,Y, transpose(streamfunc_u), transpose(streamfunc_v),xs,ys);
    
    
    % Plot the channel boundaries.
    hold on;
    plot([1; grain_x], [0.75; 0.75], 'k-');
    hold on;
    plot([1; grain_x], [grain_y + 0.25; grain_y + .25], 'k-');
    
    % Display the flow obstacle.
    obstacle_x = zeros(1, nnz(obstacle));
    obstacle_y = zeros(1, nnz(obstacle));
    k = 1;

    for (i = 1:1:numnodes_x)
        for (j = 1:1:numnodes_y)
            if (obstacle(i, j) == 1)
                obstacle_x(k) = 0.5 + (numnodes_x ./ (grain_size .* (numnodes_x - 1))) .* (i - 1);
                obstacle_y(k) = 0.5 + (numnodes_y ./ (grain_size .* (numnodes_y - 1))) .* (j - 1);
                k = k + 1;
            end
        end
    end
    
    hold on; 
    plot(obstacle_x, obstacle_y, 'r-');
    axis equal;    
    toc;
    shg;
end




