%{ 
check if relevant set functions are submodular. 
turns out area and perimeter are. but isoperimetric ratio is not.
neither is the squared version.
%}
N = 1000;
eps = .000001;
isSubmod = @(S,T,f) (f(S) + f(T)) - (f(union(S,T)) + f(intersect(S,T)));
pp = @(S) area(S)/perimeter(S); % def not submod
%pp2 = @(S) area(S)^2/perimeter(S);
pp3 = @(S) area(S)/perimeter(S)^2; % def not submod
%ppf = @(S) perimeter(S)/area(S);
%pp2f = @(S) perimeter(S)/area(S)^2;
%pp3f = @(S) perimeter(S)^2/area(S);
Ss = cell(N,1);
Ts = cell(N,1);
cs1s = cell(N,1);
rs1s = cell(N,1);
cs2s = cell(N,1);
rs2s = cell(N,1);
checks = zeros(N,4);
for i=1:N
    display(i);
    n1 = randi(20);
    cs1 = rand(n1,2);
    rs1 = randn(n1,1)/2;
    cs1s{i} = cs1;
    rs1s{i} = rs1;

    n2 = randi(20);
    cs2 = rand(n2,2);
    rs2 = randn(n2,1)/2;
    cs2s{i} = cs2;
    rs2s{i} = rs2;
    
    S = union(circlepoly_v3(rs1,cs1,.0001));
    T = union(circlepoly_v3(rs2,cs2,.0001));
    Ss{i} = S;
    Ts{i} = T;
    
    checks(i,1) = isSubmod(S,T,@perimeter);
    checks(i,2) = isSubmod(S,T,@area);
    checks(i,3) = isSubmod(S,T,pp);
    checks(i,4) = isSubmod(S,T,pp3);
%     checks(i,5) = isSubmod(S,T,pp2);
%     checks(i,6) = isSubmod(S,T,ppf);
%     checks(i,7) = isSubmod(S,T,pp2f);
%     checks(i,8) = isSubmod(S,T,pp3f);
end
triminds = any(isnan(checks),2);
checksTrimmed = checks; 
SsTrimmed = Ss;
TsTrimmed = Ts;
cs1sTrimmed = cs1s;
rs1sTrimmed = rs1s;
cs2sTrimmed = cs2s;
rs2sTrimmed = rs2s;
checksTrimmed(triminds,:)=[];
SsTrimmed(triminds,:)=[];
TsTrimmed(triminds,:)=[];
cs1sTrimmed(triminds)=[];
rs1sTrimmed(triminds)=[];
cs2sTrimmed(triminds)=[];
rs2sTrimmed(triminds)=[];
sum(checksTrimmed > -1e-6)./size(checksTrimmed,1)

failed1 = find(~(checksTrimmed(:,3) > -1e-6));
failed2 = find(~(checksTrimmed(:,4) > -1e-6));

ii = failed2(4);
S = SsTrimmed{ii}; T = TsTrimmed{ii};
cs1 = cs1sTrimmed{ii};
rs1 = rs1sTrimmed{ii};
cs2 = cs2sTrimmed{ii};
rs2 = rs2sTrimmed{ii};
fh=figure; hold all; axis equal; 
fh.WindowStyle='docked'
plot(S,'facecolor','g');
plot(T,'facecolor','r');

(pp(S) + pp(T)) - (pp(union(S,T)) + pp(intersect(S,T)))
(pp3(S) + pp3(T)) - (pp3(union(S,T)) + pp3(intersect(S,T)))


