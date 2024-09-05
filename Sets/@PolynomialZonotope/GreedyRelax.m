
function PZ = GreedyRelax(PZ,eta)

% % Construct Graph:
% Graph.V = 1:PZ.dims.h;
% Graph.E = [];
% G.Weights = [];
% for i = 1:PZ.dims.h
%     for j = 1:PZ.dims.h
%         k = find(PZ.E(:,i)-PZ.E(:,j));
%         if isscalar(k) && (mod(PZ.E(k,i),2) == 0 && mod(PZ.E(k,j),2) == 0 || PZ.E(k,j) == 0) && (PZ.E(k,i) > PZ.E(k,j) && PZ.E(k,j) >= eta)
%             Graph.E = [Graph.E, [i;j]];
%             G.Weights = [G.Weights exp(PZ.E(k,i)*log(PZ.E(k,j)/PZ.E(k,i))/(PZ.E(k,i)-PZ.E(k,j)))];
%         end
%     end
% end

% Decompose PZ:
c = sum(PZ.G(:,all(PZ.E==0,1)),2);
GI = PZ.G(:,all(sum(PZ.E,1)==1) & all(sum(~PZ.E,1)==PZ.dims.p-1));
G = PZ.G(:,~(all(sum(PZ.E,1)==1) & all(sum(~PZ.E,1)==PZ.dims.p-1)) & ~all(PZ.E==0,1));
E = PZ.E(:,~(all(sum(PZ.E,1)==1) & all(sum(~PZ.E,1)==PZ.dims.p-1)) & ~all(PZ.E==0,1));
id = PZ.id;
h = size(E,2);

% Construct Greedy Graph:
Graph.V = 1:h;
Graph.E = [];
for i = 1:h
    MinWeight = 1;
    MinEdge = [i;i];
    for j = 1:h
        k = find(E(:,i)-E(:,j));
        if isscalar(k) && (mod(E(k,i),2) == mod(E(k,j),2) || E(k,j) == 0) && (E(k,i) > E(k,j) && E(k,j) >= eta) && (exp(E(k,i)*log(E(k,j)/E(k,i))/(E(k,i)-E(k,j))) < MinWeight)
            MinWeight = exp(E(k,i)*log(E(k,j)/E(k,i))/(E(k,i)-E(k,j)));
            MinEdge = [i;j];
        end
    end
    if MinWeight ~= 1
        Graph.E = [Graph.E MinEdge];
    end
end

while ~isempty(Graph.E)

    LeafEdges = Graph.E(:,~ismember(Graph.E(1,:),Graph.E(2,:)));
    if isempty(LeafEdges)
        warning('Leaf Edges not found');
        break;
    end
    while ~isempty(LeafEdges)
        [c,G,GI,E] = Relax(c,G,GI,E,LeafEdges(1,1),LeafEdges(2,1));
        % Update the Graph:
        Graph.V = setdiff(Graph.V,LeafEdges(1,1));
        Graph.E(:,LeafEdges(1,1)==Graph.E(1,:)) = [];
        Graph.E(Graph.E>LeafEdges(1,1)) = Graph.E(Graph.E>LeafEdges(1,1))-1;
        LeafEdges(LeafEdges>LeafEdges(1,1)) = LeafEdges(LeafEdges>LeafEdges(1,1))-1;
        LeafEdges(:,1) = [];
    end
    % PZ = Relax(PZ,LeafEdges(1,1),LeafEdges(2,2));
    % Graph.E(:,all(LeafEdges(:,1)==Graph.E,1)) = []; % Remove Edges

end

% Construct PZ:
PZ = PolynomialZonotope([c GI G],[blkdiag(zeros(size(c,2),size(c,2)),eye(size(GI,2)),E)],[UniqueID(size(c,2)+size(GI,2)) id]);

end

function [c,G,GI,E] = Relax(c,G,GI,E,i,j)

k0 = find(E(:,i)-E(:,j));
if length(k0) ~= 1
    error('Invalid Indices for Relaxation');
end

etaE = exp(E(k0,i)*log(E(k0,j)/E(k0,i))/(E(k0,i)-E(k0,j)));
if all(mod(E(:,j),2)==0)
    etaG = etaE/2;
else
    etaG = etaE;
end
etac = etaE-etaG;

c = c + etac*G(:,i);
GI = [GI etaG*G(:,i)];
% E = [E(:,setdiff(1:size(G,2),[i,j])) E(:,j)];
% G = [G(:,setdiff(1:size(G,2),[i,j])) G(:,j)+G(:,i)];
E = [E(:,setdiff(1:j-1,i)) E(:,j) E(:,setdiff(j+1:end,i))];
G = [G(:,setdiff(1:j-1,i)) G(:,j)+G(:,i) G(:,setdiff(j+1:end,i))];

end
