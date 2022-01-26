function fhp
% starting time of simulation
tic;

% time steps to simulate
t_steps = 3600; 
disp("Running the simulation with " +t_steps+ " steps")

% use FHP 1 collision
use_FHP1 = true;
%use_FHP1 = false;

% use FHP 2 collision rules
use_FHP2 = true;
%use_FHP2 = false;

% use FHP 3 collision rules
use_FHP3 = true;
%use_FHP3 = false;

% board size
board_x = 500;
board_y = 200;
% the board is build up top to bottom, left to right:
% (1,1) (2,1) (3,1)
% (1,2) (2,2) (3,2)
% (1,3) (2,3) (3,3)

% 7 states possible, 6 Directions and 1 standing still
% Directions possible
%    5   6
%     \ /
%  4 - 7 - 1
%     / \
%    3   2
%
% Collision rules [in] = [out]; or = 50% probability
% FHP  I:  (1,4) = (3,6) or (2,5)
%       => (2,5) = (1,4) or (3,6)
%       => (3,6) = (1,4) or (2,5)
%          (2,4,6) = (2,5,6)
%       => (1,3,5) = (1,3,5)
%
% FHP II:  (1,7) = (2,6)
%       => ...
%          (3,5) = (7,1)
%       => ...
%          (7,1,4) = (7,3,6) or (7,2,5)
%
% FHP III: (2,3,5,6) = (1,3,4,6) or (1,2,4,5)
%       => ...
%          (2,3,5,6,7) = (1,3,4,6,7) or (1,2,4,5,7)
%       => ...
%          (2,4,6,7) = (2,4,6,7)
%       => (1,3,5,7) = (1,3,5,7)
%          (3,4,6) = (2,6,7) or (1,2,5)
%       => ...
%          (3,5,7) = (1,3,6) or (1,2,5)
%       => ...
%          (2,3,4,5,6) = (1,2,4,6,7)
%       => ...
%          (1,3,4,5,7) = (1,2,3,5,6)
%       => ...
%          (3,4,6,7) = (1,2,4,6) or (1,2,5,7)
%       => ...
%          (1,3,4,5) = (1,3,6,7) or (1,2,5,7)
%       => ...



% initialize the board
% in each cell can (theoretically) be up to seven particles
% the 3rd dimensions can contain a 0 or a 1 for every direction (7 =
%  standing still)
board = zeros(board_x, board_y, 7);

% defining the lattice velocities [x, y]
c1 = [1; 0];
c2 = [cos(5*pi/3); sin(5*pi/3)];
c3 = [cos(4*pi/3); sin(4*pi/3)];
c4 = [-1; 0];
c5 = [cos(2*pi/3); sin(2*pi/3)];
c6 = [cos(pi/3); sin(pi/3)];
c7 = [0; 0];

% initialize obstacles array;
obstacles = zeros(board_x, board_y);

% defining the start and end positions of the obstacle [x y]
obstacle_start = [100 75];
obstacle_end   = [100 125];

% create the obstacle
% for x = obstacle_start(1):1:obstacle_end(1)
%     for y = obstacle_start(2):1:obstacle_end(2)
%         obstacles(x, y) = 1;
%     end
% end


% create a skew
y = 75;
for x = 100:1:150
   obstacles(x,y) = 1;
   y = y + 1;
end
    


% setting up the board
for x = 1:1:board_x
    % do not include the top and bottom borders
    for y = 2:1:(board_y-1)
        if obstacles(x, y) ~= 1
            % getting the cell of the coordinate
            cur_cell = board(x, y, :);
            
            % setting the cell to to the right
            cur_cell(1) = 1;
            
            % reinserting the cell made into the board
            board(x, y, :) = cur_cell;
        end
    end
end



% do the simulation
for t = 1:1:t_steps
    % calculate the collissions within the borders
    for x = 1:1:board_x
        % do not include top and bottom border
        for y = 2:1:(board_y-1)
            
            % ignore obstacle points
            if obstacles(x, y) ~= 1
                
                % get the current cell
                cur_cell = board(x, y, :);
                
                % initializing the new cell
                n_cell = cur_cell;
                
                % calculate the amount of particles in the cell
                cur_cell_sum = sum(cur_cell);
                
                % determine the type of collision
                if cur_cell_sum == 2
                    
                    % there is a static point and one other hitting that
                    % FHP 2:
                    % (4,7) = (2,6)
                    if  use_FHP2 && cur_cell(7) == 1
                        % finding the existing direction
                        dir = find(cur_cell, 1);
                        
                        % resetting the direction the particle came from
                        n_cell(dir) = 0;
                        
                        % setting the new directions of the particles
                        n_cell(handle_over_under(dir-2)) = 1;
                        n_cell(handle_over_under(dir+2)) = 1;
                        
                        % there is now no standing still particle anymore
                        n_cell(7)   = 0;
                    else
                        % check the different opposites
                        for dir = 1:1:6
                            
                            % FHP 1:
                            % (1,4) = (3,6) or (2,5)
                            if  use_FHP1 && (cur_cell(dir) == 1) && ...
                                    (cur_cell(dir) == cur_cell(handle_over_under(dir+3)))
                                % random chance of 50% between the options
                                r = rand;
                                % for that we rotate the cells
                                if r < 0.5
                                    n_cell(1) = cur_cell(2);
                                    n_cell(2) = cur_cell(3);
                                    n_cell(3) = cur_cell(4);
                                    n_cell(4) = cur_cell(5);
                                    n_cell(5) = cur_cell(6);
                                    n_cell(6) = cur_cell(1);
                                else
                                    n_cell(1) = cur_cell(6);
                                    n_cell(2) = cur_cell(1);
                                    n_cell(3) = cur_cell(2);
                                    n_cell(4) = cur_cell(3);
                                    n_cell(5) = cur_cell(4);
                                    n_cell(6) = cur_cell(5);
                                end
                                n_cell(7) = cur_cell(7);
                                
                                % exit loop as there was one pair found
                                break
                                
                            end
                                
                                % FHP 2:
                                % (3,5) = (1,7)
                            if use_FHP2 && (cur_cell(dir)==1) && ...
                                    (cur_cell(dir) == cur_cell(handle_over_under(dir+2)))
                                % as the loop starts at 1 it is always
                                %  going to find the lower indexed
                                %  direction first. Therefore it is
                                %  sufficient to only calculate 2 backwards
                                
                                % resetting the incoming
                                n_cell(dir) = 0;
                                n_cell(handle_over_under(dir+2)) = 0;
                                
                                % setting the outgoing and the stationary
                                n_cell(handle_over_under(dir-2)) = 1;
                                n_cell(7) = 1;
                                
                                % exit the loop as there was a pair found
                                break
                            end
                        end
                    end
                    
                elseif cur_cell_sum == 3
                    
                    % FHP 1:
                    % if either 1,3,5 are set or not set (then 2,4,6 is set)
                    if use_FHP1 && (cur_cell(1) == cur_cell(3)) && (cur_cell(3) == cur_cell(5))
                        % leave the state as it is
                        n_cell = cur_cell;
                        
                    else
                        for dir = 1:1:6
                            
                            %FHP 2:
                            % opposing particles hit with one additional standing
                            %  still in the middle
                            % ... indicates a new line
                            % (1,4,7) = (3,6,7) or (2,5,7)
                            if use_FHP2 && (cur_cell(7) == 1) && (cur_cell(dir) == 1) && ...
                                    (cur_cell(dir) == cur_cell(handle_over_under(dir+3)))
                                
                                % random chance of 50% between the options
                                r = rand;
                                % for that we rotate the cells
                                if r < 0.5
                                    n_cell(1) = cur_cell(2);
                                    n_cell(2) = cur_cell(3);
                                    n_cell(3) = cur_cell(4);
                                    n_cell(4) = cur_cell(5);
                                    n_cell(5) = cur_cell(6);
                                    n_cell(6) = cur_cell(1);
                                else
                                    n_cell(1) = cur_cell(6);
                                    n_cell(2) = cur_cell(1);
                                    n_cell(3) = cur_cell(2);
                                    n_cell(4) = cur_cell(3);
                                    n_cell(5) = cur_cell(4);
                                    n_cell(6) = cur_cell(5);
                                end
                                n_cell(7) = cur_cell(7);
                                
                                % exit loop as there was one pair found
                                break
                            end
                            
                            % FHP 3:
                            % (3,5,7) = (1,3,5) or (1,2,5)
                            if use_FHP3 && (cur_cell(7) == 1) && (cur_cell(dir) == 1) && ...
                                    (cur_cell(dir) == cur_cell(handle_over_under(dir+2)))
                                
                                % random chance between the two outcomes
                                r = rand;
                                
                                % resetting the cell
                                n_cell = [0 0 0 0 0 0 0];
                                if r < 0.5
                                    n_cell(dir) = 1;
                                    n_cell(handle_over_under(dir+3)) = 1;
                                    n_cell(handle_over_under(dir+4)) = 1;
                                else
                                    n_cell(handle_over_under(dir+2)) = 1;
                                    n_cell(handle_over_under(dir-1)) = 1;
                                    n_cell(handle_over_under(dir-2)) = 1;
                                    
                                end
                                
                                % exit loop as there was one pair found
                                break
                            end
                            
                            % FHP 3:
                            % (3,4,6) = (2,6,7) or (1,2,5)
                            if use_FHP3 && (cur_cell(dir) == 1) && ...
                                    (cur_cell(dir) == cur_cell(handle_over_under(dir+1))) && ...
                                    (cur_cell(dir) == cur_cell(handle_over_under(dir+3)))
                                
                                % random chance between the two outcomes
                                r = rand;
                                
                                % resetting the cell
                                n_cell = [0 0 0 0 0 0 0];
                                if r < 0.5
                                    n_cell(7) = 1;
                                    n_cell(handle_over_under(dir-1)) = 1;
                                    n_cell(handle_over_under(dir+3)) = 1;
                                else
                                    n_cell(handle_over_under(dir+2)) = 1;
                                    n_cell(handle_over_under(dir-1)) = 1;
                                    n_cell(handle_over_under(dir-2)) = 1;
                                    
                                end
                                
                                % exit loop as there was one pair found
                                break
                            end
                            
                        end
                    end
                    
                    % FHP 3: (only FHP 3 can have 4 particles in one cell)
                elseif cur_cell_sum == 4 && use_FHP3
                    
                    for dir = 1:1:6
                        
                        % going here for the directions not set as it
                        % is more understandable
                        % (2,3,5,6) = (1,3,4,6) or (1,2,4,5)
                        if (cur_cell(7) == 0) && (cur_cell(dir) == 0) && ...
                                (cur_cell(dir) == cur_cell(handle_over_under(dir+3)))
                            
                            % randomly choose one two options
                            r = rand;
                            
                            % rotate the directions by one
                            % clockwise
                            if r < 0.5
                                n_cell(1) = cur_cell(6);
                                n_cell(2) = cur_cell(1);
                                n_cell(3) = cur_cell(2);
                                n_cell(4) = cur_cell(3);
                                n_cell(5) = cur_cell(4);
                                n_cell(6) = cur_cell(5);
                                
                                % counterclockwise
                            else
                                n_cell(1) = cur_cell(2);
                                n_cell(2) = cur_cell(3);
                                n_cell(3) = cur_cell(4);
                                n_cell(4) = cur_cell(5);
                                n_cell(5) = cur_cell(6);
                                n_cell(6) = cur_cell(1);
                                
                            end
                            n_cell(7) = 0;
                               
                            % exit loop as there was one pair found
                            break
                        end
                        
                        % (2,4,6,7) = (2,4,6,7) 
                        if (cur_cell(7) == 1) && ...
                                (cur_cell(1) == cur_cell(3)) && ...
                                (cur_cell(3) == cur_cell(5))
                            
                            
                            % in this case let the cell as it is
                            n_cell = cur_cell;
                            
                            
                            % exit loop as there was one pair found
                            break
                        end
                        
                        % (3,4,6,7) = (1,2,3,6) or (1,2,5,7)
                        if (cur_cell(7)==1) && (cur_cell(dir)==1) && ...
                                (cur_cell(dir) == cur_cell(handle_over_under(dir+1))) && ...
                                (cur_cell(dir) == cur_cell(handle_over_under(dir+3)))
                            
                            
                            % random chance for both outcomes
                            r = rand;
                            
                            % resetting the new cell
                            n_cell = [0 0 0 0 0 0 0];
                            if r < 0.5
                                n_cell(handle_over_under(dir+1)) = 1;
                                n_cell(handle_over_under(dir-1)) = 1;
                                n_cell(handle_over_under(dir+3)) = 1;
                                n_cell(handle_over_under(dir-2)) = 1;
                            else
                                n_cell(7) = 1;
                                n_cell(handle_over_under(dir-1)) = 1;
                                n_cell(handle_over_under(dir+2)) = 1;
                                n_cell(handle_over_under(dir-2)) = 1;
                                
                            end
                            
                            % exit loop as there was one pair found
                            break
                        end
                        
                        % (1,3,4,5) = (1,3,6,7) or (1,2,5,7)
                        if (cur_cell(7) == 0) && (cur_cell(dir) == 1) && ...
                                (cur_cell(dir) == cur_cell(handle_over_under(dir+2))) && ...
                                (cur_cell(dir) == cur_cell(handle_over_under(dir+3))) && ...
                                (cur_cell(dir) == cur_cell(handle_over_under(dir+4)))
                            
                            % random chance for the output
                            r = rand;
                            
                            % resetting the new cell
                            n_cell = [0 0 0 0 0 0 0];
                            if r < 0.5
                                n_cell(handle_over_under(dir-1)) = 1;
                                n_cell(handle_over_under(dir+2)) = 1;
                                
                            else
                                n_cell(handle_over_under(dir+1)) = 1;
                                n_cell(handle_over_under(dir-2)) = 1;
                                
                            end
                            
                            n_cell(7) = 1;
                            n_cell(dir) = 1;
                            
                            % exit loop as there was one pair found
                            break
                        end
                        
                    end
                    
                    % FHP 3: (only FHP 3 can have 5 particles in one cell)
                elseif cur_cell_sum == 5 && use_FHP3
                    
                    for dir = 1:1:6
                        
                        % looking after the 0s here, as it is easier
                        % (2,3,5,6,7) = (1,3,4,6,7) or (1,2,4,5,7)
                        if (cur_cell(7)==1) && (cur_cell(dir)==0) && ...
                                (cur_cell(dir) == cur_cell(handle_over_under(dir+3)))
                            
                            
                            % randomly choose one two options
                            r = rand;
                            
                            % rotate the directions by one
                            % clockwise
                            if r < 0.5
                                n_cell(1) = cur_cell(6);
                                n_cell(2) = cur_cell(1);
                                n_cell(3) = cur_cell(2);
                                n_cell(4) = cur_cell(3);
                                n_cell(5) = cur_cell(4);
                                n_cell(6) = cur_cell(5);
                                
                                % counterclockwise
                            else
                                n_cell(1) = cur_cell(2);
                                n_cell(2) = cur_cell(3);
                                n_cell(3) = cur_cell(4);
                                n_cell(4) = cur_cell(5);
                                n_cell(5) = cur_cell(6);
                                n_cell(6) = cur_cell(1);
                                
                            end
                            n_cell(7) = 1;
                            
                            % exit loop as there was one pair found
                            break
                        end
                        
                        % going for the 0 here for obvious reasons
                        % (2,3,4,5,6) = (1,2,4,6,7)
                        if (cur_cell(7)==0) && (cur_cell(dir)==0)
                            
                            
                            % setting every direction + standing still
                            n_cell = [1 1 1 1 1 1 1];
                            
                            n_cell(handle_over_under(dir+2)) = 0;
                            n_cell(handle_over_under(dir-2)) = 0;
                            
                            % exit loop as there was one pair found
                            break
                        end
                        
                        % the reversed rule to the previous one
                        % (1,3,4,5,7) = (1,2,3,5,6)
                        if (cur_cell(7)==1) && (cur_cell(dir)==0) && ...
                                (cur_cell(dir) == cur_cell(handle_over_under(dir-2)))
                            
                            % setting every direction without standing still
                            n_cell = [1 1 1 1 1 1 0];
                            
                            n_cell(handle_over_under(dir+2)) = 0;
                            
                            % exit loop as there was one pair found
                            break
                            
                        end
                        
                        
                    end
                    
                    
                    % no collision
                else
                    
                end
                % an obstacle is there
            else
                % inverting the velocities
                n_cell(1) = cur_cell(4);
                n_cell(2) = cur_cell(5);
                n_cell(3) = cur_cell(6);
                n_cell(4) = cur_cell(1);
                n_cell(5) = cur_cell(2);
                n_cell(6) = cur_cell(3);
                n_cell(7) = cur_cell(7);
            end
            % set the new cell (if nothing found, this is the same)
            board(x, y, :) = n_cell;
            
        end
    end
    
    % Calculate the collissions with the borders
    for x = 1:1:board_x
        % inverting the velocities
        board(x, 1, :) = [board(x, 1, 4) board(x, 1, 5) board(x, 1, 6) ...
            board(x, 1, 1) board(x, 1, 2) board(x, 1, 3) ...
            board(x, 1, 7)];
        
        board(x, board_y, :) = [board(x, board_y, 4) board(x, board_y, 5) ...
            board(x, board_y, 6) board(x, board_y, 1) ...
            board(x, board_y, 2) board(x, board_y, 3) ...
            board(x, board_y, 7)];
    end
    
    
    % create the new board for the next step
    n_board = zeros(board_x, board_y, 7);
    
    % moving every particle
    for x = 1:1:board_x
        % also include the borders
        for y = 1:1:board_y
            
            % get the current cell
            cur_cell = board(x, y, :);
            
            % Step into direction 1:
            if x == board_x
                n_x = 1;
            else
                n_x = x + 1;
            end
            n_y = y;
            
            % get the new cell
            n_cell = n_board(n_x, n_y, :);
            % set the 1. component
            n_cell(1) = cur_cell(1);
            
            % if the comes back into the left side it should reset to only
            % going into direction 1
            if n_x == 1 && y ~= 1 && y ~= board_y
                n_cell = [1 0 0 0 0 0 0];
            end 
            % set the new board
            n_board(n_x, n_y, :) = n_cell;
            
            
            
            % Step into direction 2:
            % if in the border, do not move downwards
            if( y ~= board_y)
                n_y = y + 1;
                
                if mod(y,2) == 0
                    if x == board_x
                        n_x = 1;
                    else
                        n_x = x + 1;
                    end
                else
                    n_x = x;
                end
                
                % get the new cell
                n_cell = n_board(n_x, n_y, :);
                % set the 1. component
                n_cell(2) = cur_cell(2);
                
                % if the comes back into the left side it should reset to only
                % going into direction 1
                if n_x == 1 && y ~= 1 && y ~= board_y
                    n_cell = [1 0 0 0 0 0 0];
                end
                % set the new board
                n_board(n_x, n_y, :) = n_cell;
                
            end
            
            
            % Step into direction 3:
            % if in the border, do not move downwards
            if( y ~= board_y)
                n_y = y + 1;
                
                if mod(y,2) == 1
                    if x == 1
                        n_x = board_x;
                    else
                        n_x = x - 1;
                    end
                else
                    n_x = x;
                end
                
                % get the new cell
                n_cell = n_board(n_x, n_y, :);
                % set the 1. component
                n_cell(3) = cur_cell(3);
                % set the new board
                n_board(n_x, n_y, :) = n_cell;
                
            end
            
            % Step into direction 4:
            if x == 1
                n_x = board_x;
            else
                n_x = x - 1;
            end
            n_y = y;
            
            % get the new cell
            n_cell = n_board(n_x, n_y, :);
            
            % set the 1. component
            n_cell(4) = cur_cell(4);
            % set the new board
            n_board(n_x, n_y, :) = n_cell;
            
            % Step into direction 5:
            % if in the border, do not move upwards
            if y ~= 1
                n_y = y - 1;
                
                if mod(y,2) == 1
                    if x == 1
                        n_x = board_x;
                    else
                        n_x = x - 1;
                    end
                else
                    n_x = x;
                end
                
                
                % get the new cell
                n_cell = n_board(n_x, n_y, :);
                % set the 1. component
                n_cell(5) = cur_cell(5);
                % set the new board
                n_board(n_x, n_y, :) = n_cell;
                
                
            end
            
            % Step into direction 6:
            % if in the border, do not move upwards
            if y ~= 1
                n_y = y - 1;
                
                if mod(y,2) == 0
                    if x == board_x
                        n_x = 1;
                    else
                        n_x = x + 1;
                    end
                else
                    n_x = x;
                end
                
                % get the new cell
                n_cell = n_board(n_x, n_y, :);
                % set the 1. component
                n_cell(6) = cur_cell(6);
                
                % if the comes back into the left side it should reset to only
                % going into direction 1
                if n_x == 1 && y ~= 1 && y ~= board_y
                    n_cell = [1 0 0 0 0 0 0];
                end
                % set the new board
                n_board(n_x, n_y, :) = n_cell;
                
            end
            
            
            % Step into direction 7:
            n_x = x;
            n_y = y;
            
            % get the new cell
            n_cell = n_board(n_x, n_y, :);
            % set the 1. component
            n_cell(7) = cur_cell(7);
            % set the new board
            n_board(n_x, n_y, :) = n_cell;
            
            
            
            
        end
    end
    
    % set the new board
    board = n_board;
    
    % provide feedback to the user that system is not frozen
    if mod(t, 5) == 0
        disp(t)
    end
    
%include the output in the simulation and therefore plot it everytime
end

    % divide for course-graining
    divider = 5;
    sub_board_x = board_x/divider;
    sub_board_y = board_y/divider;

    % initialize the average velocity vectors
    avg_velocity_x_coordinates = zeros(1, sub_board_x * sub_board_y);
    avg_velocity_y_coordinates = zeros(1, sub_board_x * sub_board_y);
    avg_velocity_horizontal  = zeros(1, sub_board_x * sub_board_y);
    avg_velocity_vertical    = zeros(1, sub_board_x * sub_board_y);

    % initialize the particle density matrix
    density = zeros(sub_board_x, sub_board_y);

    % streamline parameter
    avg_streamline_horizontal = zeros(sub_board_x, sub_board_y);
    avg_streamline_vertical   = zeros(sub_board_x, sub_board_y);
    streamline_ys = 1:5:sub_board_x;
    streamline_xs = ones(1,length(streamline_ys));

    % introducing a counter
    cur_value = 1;
    % running over all the data
    for x = 1:1:sub_board_x
        % calculate the boundaries of x
        x_b_do = (x-1)*divider + 1;
        x_b_up = x*divider;
        for y = 1:1:sub_board_y
            %calculate the boundaries of y
            y_b_do = (y-1)*divider + 1;
            y_b_up = y*divider;

            % get the numbers of each direction
            nums = zeros(1, 7);

            nums(1) = sum(sum(board(x_b_do:1:x_b_up, y_b_do:1:y_b_up, 1)));
            nums(2) = sum(sum(board(x_b_do:1:x_b_up, y_b_do:1:y_b_up, 2)));
            nums(3) = sum(sum(board(x_b_do:1:x_b_up, y_b_do:1:y_b_up, 3)));
            nums(4) = sum(sum(board(x_b_do:1:x_b_up, y_b_do:1:y_b_up, 4)));
            nums(5) = sum(sum(board(x_b_do:1:x_b_up, y_b_do:1:y_b_up, 5)));
            nums(6) = sum(sum(board(x_b_do:1:x_b_up, y_b_do:1:y_b_up, 6)));
            nums(7) = sum(sum(board(x_b_do:1:x_b_up, y_b_do:1:y_b_up, 7)));

            % calculate average numsity
            avg_velocity = (1/divider.^2) * (nums(1)*c1 + nums(2)*c2 + ...
                nums(3)*c3 + nums(4)*c4 + nums(5)*c5 + nums(6)*c6 + nums(7)*c7);

            % calculate the density
            density(x, y) = sum(nums) / ( 7 * (x_b_up - x_b_do) * ...
                (y_b_up - y_b_do) );

            % storing the average velocity
            avg_velocity_horizontal(cur_value) = avg_velocity(1);
            avg_velocity_vertical(cur_value)   = avg_velocity(2);

            % storing the coordinates
            avg_velocity_x_coordinates(cur_value) = x;
            avg_velocity_y_coordinates(cur_value) = y;

            % storing the streamline values
            avg_streamline_horizontal(x,y) = avg_velocity(1);
            avg_streamline_vertical(x,y) = avg_velocity(2);

            % incrementing the counter
            cur_value = cur_value + 1;
        end
    end

    % set te size of the figure
    set(gcf, 'position', [10, 10, 1500, 900])
    
    % create a tiled layout for showing 2 graphs at the same time
    tiledlayout(2,1)

    % plot the avegare velocity vectors in the first tile
    nexttile
    title('Flow')
    quiver(avg_velocity_x_coordinates, avg_velocity_y_coordinates, ...
        avg_velocity_horizontal, avg_velocity_vertical);

    % plotting the streamline into the avg velocity vectors graph
    streamline(avg_streamline_horizontal.', avg_streamline_vertical.', streamline_xs, streamline_ys);

    % Plot the channel boundaries.
    hold on;
    plot([1; sub_board_x], [0.75; 0.75], 'k-');
    hold on;
    plot([1; sub_board_x], [sub_board_y + 0.25; sub_board_y + .25], 'k-');

    %plot the obstacle
    obstacle_x = zeros(1, nnz(obstacles));
    obstacle_y = zeros(1, nnz(obstacles));
    k = 1;
    for x = 1:1:board_x
        for y = 1:1:board_y
            if obstacles(x, y) == 1
                obstacle_x(k) = 0.5 + (board_x ./ (divider .* (board_x - 1))) .* (x - 1);
                obstacle_y(k) = 0.5 + (board_y ./ (divider .* (board_y - 1))) .* (y - 1);
                k = k + 1;
            end
        end
    end

    hold on;
    plot(obstacle_x, obstacle_y, 'r-');
    axis equal;


    % plotting the particle density in the next tile
    nexttile
    title('Density')
    contourf(density.');
    colormap("cool");
    colorbar();

    % plotting the streamline into the density graph
    streamline(avg_streamline_horizontal.', avg_streamline_vertical.', streamline_xs, streamline_ys);

    % ending time of simulation
    toc;

    % show the current figure
    shg;
    % for the continious graph showing
    %end
end

% handle over- and underflow
function n_dir = handle_over_under(dir)
n_dir = dir;
if dir > 6
    n_dir = dir - 6;
elseif dir < 1
    n_dir = dir + 6;
end

end
