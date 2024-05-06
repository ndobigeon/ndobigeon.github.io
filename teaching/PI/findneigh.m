%f indneigh.m
function A = findneigh(i,j,m)
%find neighbors of the (i,j)-th site of an m by m grid
l = 1;
if i==l
    if j==l
        A = [l,l;l,2;2,l];
    elseif j==m
        A = [l,m;l,m-l;2,m];
    else
    A = [l,j;l,j-l;l,j+l;2,j];
    end
elseif i==m
    if j==l
        A = [m,l;m,2;m-l,l];
    elseif j==m
        A = [m,m;m,m-l;m-l,m];
    else
        A = [m,j;m,j-l;m,j+l;m-l,j];
    end
else
    if j==l
        A = [i,l;i,2;i-l,l;i+i,l];
    elseif j==m
        A = [ i , m ; i , m - l ; i + l , m ; i - l , m ] ;
    else
        A = [ i , j ; i , j - l ; i , j + l ; i + l , j ; i - l , j ] ;
    end
end