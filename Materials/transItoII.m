function AII=transItoII(AI)
% To transform a property tensor A in type I into type II
%
M=sparse(8,8);
M(1,1)=1;
M(2,5)=1;
M(3,2)=1;
M(4,6)=1;
M(5,3)=1;
M(6,7)=1;
M(7,4)=1;
M(8,8)=1;
%
N=sparse(8,8);
N(1,1)=1;
N(2,3)=1;
N(3,5)=1;
N(4,7)=1;
N(5,2)=1;
N(6,4)=1;
N(7,6)=1;
N(8,8)=1;
%
I=eye(8,8);
%
P=I+M-N;
AII=P'*AI*P;
end
