function epsilon = cenProlong(r)

M = size(r,1)-2; N = size(r,2)-2;
epsilon = zeros(2*M+2,2*N+2);


for j = 2:N+1
    for i = 2:M+1
        epsilon(2*i-2,2*j-2) = r(i,j);
        epsilon(2*i-1,2*j-2) = r(i,j);
        epsilon(2*i-2,2*j-1) = r(i,j);
        epsilon(2*i-1,2*j-1) = r(i,j);
    end
    %CHANGE BOUNDARY CONDITIONS FOR EACH PROBLEM
    epsilon = bcGS(epsilon);
end