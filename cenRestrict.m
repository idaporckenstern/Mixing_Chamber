function epsilon = cenRestrict(r)

M = size(r,1)/2-1; N = size(r,2)/2-1;
epsilon = zeros(M+2,N+2);

for j = 2:N+1
    for i = 2:M+1
        epsilon(i,j) = 1/4*(r(2*i-2,2*j-2)+r(2*i-1,2*j-2)+r(2*i-2,2*j-1)+r(2*i-1,2*j-1));
    end
end
%CHANGE BOUNDARY CONDITIONS FOR EVERY PROBLEM
epsilon = bcGS(epsilon);
end