% Dimitrios Piretzidis (2021). Surface Spherical Harmonic Functions Visualization (https://www.mathworks.com/matlabcentral/fileexchange/44869-surface-spherical-harmonic-functions-visualization), MATLAB Central File Exchange. Retrieved May 29, 2021.

function [ output ] = spherical_harmonics(n,m,theta,lamda)
[lamda, theta] = meshgrid(lamda,theta);
output = zeros(size(lamda));
for i = 1: size(lamda,1)
    for j = 1:size(lamda,1)
        output(i,j) = NLegendre_Fun(n, m, cos(theta(i,j))) * cos(m * lamda(i,j));
    end
end
end

function [output] = NLegendre_Fun(n,m,t)
output=zeros(n+1,m+1);
output(1,1)=1;
output(2,1)=sqrt(3)*t;
output(2,2)=sqrt(3)*sqrt(1-t^2);
%Recursively compute the elements of the diagonal until degree m
for i=3:m+1
    nn=i-1;
    output(i,i)=sqrt((2*nn+1)/(2*nn))*sqrt(1-t^2)*output(i-1,i-1);
end
%Recursively compute P(m+1,m)
if n>m
    nn=m+1;
    output(m+2,m+1)=sqrt(2*nn+1)*t*output(m+1,m+1);
end
%Recursively compute P(n,m)
if n>m+1
    for i=m+3:n+1
        nn=i-1;
        output(i,m+1)=sqrt(((2*nn-1)*(2*nn+1))/((nn-m)*(nn+m)))*t*output(i-1,m+1) - sqrt(((2*nn+1)*(nn+m-1)*(nn-m-1))/((2*nn-3)*(nn+m)*(nn-m)))*output(i-2,m+1);
    end
end
output=output(n+1,m+1);
end
