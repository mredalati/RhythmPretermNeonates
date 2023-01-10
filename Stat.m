%% t-test with Bonferroni correction

% Duple/Triple Rhythm
load('DupleTripleRhythmMat.mat');
[h1,p1,sigPairs1] = ttest_bonf(DupleRhyt,[1 6;2 6;3 6;4 6;5 6]);

% Quadruple Rhythm
load('QuadrupleRhythmMat.mat');
[h2,p2,sigPairs2] = ttest_bonf(QuadRhyt,[1 11;2 11;3 11;4 11;5 11;6 11;7 11;8 11;9 11;...
    10 11]);

%% Rayleigh test
load('Phase.mat');
H1 = zeros(128,size(Phase{1},3));
H2 = zeros(128,size(Phase{2},3));
for i=1:128
    for j=1:size(Phase{1},3)
        P1(i,j) = circ_rtest(Phase{1}(i,:,j));   % Duple\Triple Rhythm
    end
    for j=1:size(Phase{2},3)
        P2(i,j) = circ_rtest(Phase{2}(i,:,j));   % Quadruple Rhythm
    end
end
H1(P1<.05) = 1;
H2(P2<.05) = 1;