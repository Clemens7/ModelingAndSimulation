clear
clc
l = 50;
h = 20;
A = strings(h,l);
A(:,1) = 'r';
A(10,10) = 'wall';
A(1,:) = 'wall';
A(end,:) = 'wall';
%A(70:140,100:105) = 'wall';
%A(4:5,4:5) = 'wall';
B = A; %B is update matrix
C = B; %C is used to reset update matrix at end of loop

for r=1:100 %loop for solving until we have a steady state
    for i=2:h-1 %loop for propagation in h
        for j=2:l-1 %loop for propagation in l
            if A(i,j) ~= 'wall'
                if contains(A(i-1,j), 'd') & contains(A(i+1,j), 'u'); 
                    B(i,j) = strcat(B(i,j),'lr'); 
                elseif contains(A(i+1,j), 'u') & not(contains(A(i-1,j), 'd')) & A(i-1,j)~= 'wall';;
                    B(i,j) = strcat(B(i,j),'u');
                elseif contains(A(i-1,j), 'd') & not(contains(A(i+1,j), 'u')) & A(i+1,j)~= 'wall';;
                    B(i,j) = strcat(B(i,j),'d');  
                elseif A(i-1,j) == 'wall' & contains(A(i,j), 'u');
                    B(i,j) = strcat(B(i,j),'d');
                elseif A(i+1,j) == 'wall' & contains(A(i,j), 'd');
                    B(i,j) = strcat(B(i,j),'u');                    
                end
                if contains(A(i,j+1), 'l')&contains(A(i,j-1), 'r'); 
                    B(i,j) = strcat(B(i,j),'du');
                elseif contains(A(i,j+1), 'l') & not(contains(A(i,j-1), 'r')) & A(i,j-1)~= 'wall';;
                    B(i,j) = strcat(B(i,j),'l');        
                elseif contains(A(i,j-1), 'r') & not(contains(A(i,j+1), 'l')) & A(i,j+1)~= 'wall';
                    B(i,j) = strcat(B(i,j),'r');  
                elseif contains(A(i,j+1), 'wall') & contains(A(i,j), 'r');
                    B(i,j) = strcat(B(i,j),'l'); 
                elseif contains(A(i,j-1), 'wall') & contains(A(i,j), 'l');
                    B(i,j) = strcat(B(i,j),'r');
                end
            end
        
        end
    end
    A = B;
    
    B = C;
end

