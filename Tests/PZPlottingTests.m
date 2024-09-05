
G = [2 1 0;0 1 2];
E = [2 1 0;0 1 2];

PZ = PolynomialZonotope(G,E);

N = 10000;
r = sqrt(rand(1,N));
t1 = rand(1,N)*2*pi;
alpha = [r.*cos(t1);r.*sin(t1)];
P = zeros(2,N);
for i = 1:N
    P(:,i) = G*prod(alpha(:,i).^E,1)';
end

figure; hold on;
plot(PZ);
scatter(P(1,:),P(2,:),'b.');

% Compute the Polynomial Map:
E = PolynomialZonotope(Ellipsoid([1 0;0 1],[0;0]));
% doesnt work very well

%% Try Taking the Polymap of a polygon:


%% Cut Splits:
Splits{1} = [-1 1;-1 1];
r = MaxSplitDirection(PZ);
PZSplit = SplitMax(PZ);
NewSplits{1} = Splits{1}; 
NewSplits{2} = Splits{1};
NewSplits{1}(r,:) = [mean(Splits{1}(r,:)) Splits{1}(r,2)];
NewSplits{2}(r,:) = [Splits{1}(r,1) mean(Splits{1}(r,:))];
Splits = NewSplits;
% PolySplit = Polygon(PZSplit);
PolySplit = arrayfun(@(obj) obj.Polygon,PZSplit);
PolyFull = SetOperation(PolySplit,'union');

for i = 1:15
    
    PZSplitNew = [];
    PolySplitNew = [];
    NewSplits = {};
    for j = 1:length(PZSplit)
    
        if any(sum(ismembertol(PolySplit(j).Vertices,PolyFull.Vertices,1e-10),2) == 2)
            % [SplitMaxers,r] = SplitMax(PZSplit(j));
            % PZSplitNew = [PZSplitNew SplitMaxers];
            % PolySplitNew = [PolySplitNew arrayfun(@(obj) obj.Polygon,PZSplitNew(end-1:end))];
            % NewSplits{end+1} = Splits{j}; NewSplits{end}(r,:) = [mean(Splits{j}(r,:)) Splits{j}(r,2)];
            % NewSplits{end+1} = Splits{j}; NewSplits{end}(r,:) = [Splits{j}(r,1) mean(Splits{j}(r,:))];
            % if all(vecnorm(NewSplits{end-1}))
            r = MaxSplitDirection(PZSplit(j));
            Split = Splits{j};
            Split(r,:) = [mean(Splits{j}(r,:)) Splits{j}(r,2)];
            if any(vecnorm([Split(1,1) Split(2,1);Split(1,2) Split(2,1);Split(1,2) Split(2,2);Split(1,1) Split(2,2)],2,2) <=1)
                PZSplitNew = [PZSplitNew Subset(PZSplit(j),r,0,1)];
                PolySplitNew = [PolySplitNew PZSplitNew(end).Polygon];
                NewSplits = [NewSplits, Split];
            end
            Split(r,:) = [Splits{j}(r,1) mean(Splits{j}(r,:))];
            if any(vecnorm([Split(1,1) Split(2,1);Split(1,2) Split(2,1);Split(1,2) Split(2,2);Split(1,1) Split(2,2)],2,2) <=1)
                PZSplitNew = [PZSplitNew Subset(PZSplit(j),r,-1,0)];
                PolySplitNew = [PolySplitNew PZSplitNew(end).Polygon];
                NewSplits = [NewSplits, Split];
            end

        else
            PZSplitNew = [PZSplitNew PZSplit(j)];
            PolySplitNew = [PolySplitNew PolySplit(j)];
            NewSplits = [NewSplits Splits{j}];
        end

    end
    PolyFull = SetOperation(PolySplitNew,'union');
    PZSplit = PZSplitNew;
    PolySplit = PolySplitNew;
    Splits = NewSplits;

    hold on;
    for k = 1:length(PolySplit)
        plot(PolySplit(k));
    end

    % hold on;
    % for k = 1:length(Splits)
    %     plot(polyshape([Splits{k}(1,1) Splits{k}(1,2) Splits{k}(1,2) Splits{k}(1,1)],[Splits{k}(2,1) Splits{k}(2,1) Splits{k}(2,2) Splits{k}(2,2)]));
    % end

end

figure;
p = plot(PolyFull);

%% Alternative Circle Version:
G = [1 0 -1/4 0 -1/32 0 -1/128 0 -5/2048 0 -7/8192 0;0 1 0 -1/4 0 -1/32 0 -1/128 0 -5/2048 0 -7/8192];
E = [1 0 1 2 1 4 1 6 1 8 1 10;0 1 2 1 4 1 6 1 8 1 10 1];
PZ = PolynomialZonotope(G(:,1:6),E(:,1:6));
plot(PZ);

pz = polyZonotope([0;0],G(:,1:6),[],E(:,1:6));
plot(polyMap(pz,[2 1 0;0 1 2],[2 1 0;0 1 2]),[1,2],'splits',15);