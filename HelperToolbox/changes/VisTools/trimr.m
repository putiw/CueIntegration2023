function Y=trimr(X, N)%% TRIMR	trim N columns from the right of an array%% usage:%   Y = TRIMR(X, N) trims N elements from the right of the array,%% see also: trimb, trimt, triml, trim%% Lawrence K. Cormack% history:% ??/??/19??    lkc     wrote it% 07/04/2002    lkc     added and formated comments[R C]=size(X);s=C-N+1;X(:,s:C)=[];Y=X;