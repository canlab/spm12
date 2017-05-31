function DK = define_Kurtosistensor(Asym,MSK2)
% S.Mohammadi 22/10/2012

DK              = zeros(3,3,3,3,numel(MSK2));

% define Kurtosis
DK(1,1,1,1,:)   = Asym(MSK2,7); % 1
DK(2,2,2,2,:)   = Asym(MSK2,8); % 2
DK(3,3,3,3,:)   = Asym(MSK2,9); % 3

DK(1,1,1,2,:)   = Asym(MSK2,10); % 4
DK(1,1,2,1,:)   = Asym(MSK2,10); % 4
DK(1,2,1,1,:)   = Asym(MSK2,10); % 4
DK(2,1,1,1,:)   = Asym(MSK2,10); % 4

DK(1,1,1,3,:)   = Asym(MSK2,11); % 5
DK(1,1,3,1,:)   = Asym(MSK2,11); % 5
DK(1,3,1,1,:)   = Asym(MSK2,11); % 5
DK(3,1,1,1,:)   = Asym(MSK2,11); % 5

DK(2,2,2,1,:)   = Asym(MSK2,12); % 6
DK(2,2,1,2,:)   = Asym(MSK2,12); % 6
DK(2,1,2,2,:)   = Asym(MSK2,12); % 6
DK(1,2,2,2,:)   = Asym(MSK2,12); % 6

DK(2,2,2,3,:)   = Asym(MSK2,13); % 7
DK(2,2,3,2,:)   = Asym(MSK2,13); % 7
DK(2,3,2,2,:)   = Asym(MSK2,13); % 7
DK(3,2,2,2,:)   = Asym(MSK2,13); % 7

DK(3,3,3,1,:)   = Asym(MSK2,14); % 8
DK(3,3,1,3,:)   = Asym(MSK2,14); % 8
DK(3,1,3,3,:)   = Asym(MSK2,14); % 8
DK(1,3,3,3,:)   = Asym(MSK2,14); % 8

DK(3,3,3,2,:)   = Asym(MSK2,15); % 9
DK(3,3,2,3,:)   = Asym(MSK2,15); % 9
DK(3,2,3,3,:)   = Asym(MSK2,15); % 9
DK(2,3,3,3,:)   = Asym(MSK2,15); % 9

DK(1,1,2,2,:)   = Asym(MSK2,16); % 10
DK(1,2,1,2,:)   = Asym(MSK2,16); % 10
DK(1,2,2,1,:)   = Asym(MSK2,16); % 10
DK(2,1,2,1,:)   = Asym(MSK2,16); % 10
DK(2,1,1,2,:)   = Asym(MSK2,16); % 10
DK(2,2,1,1,:)   = Asym(MSK2,16); % 10

DK(1,1,3,3,:)   = Asym(MSK2,17); % 11
DK(1,3,1,3,:)   = Asym(MSK2,17); % 11
DK(1,3,3,1,:)   = Asym(MSK2,17); % 11
DK(3,1,3,1,:)   = Asym(MSK2,17); % 11
DK(3,1,1,3,:)   = Asym(MSK2,17); % 11
DK(3,3,1,1,:)   = Asym(MSK2,17); % 11

DK(2,2,3,3,:)   = Asym(MSK2,18); % 12
DK(2,3,2,3,:)   = Asym(MSK2,18); % 12
DK(2,3,3,2,:)   = Asym(MSK2,18); % 12
DK(3,2,3,2,:)   = Asym(MSK2,18); % 12
DK(3,2,2,3,:)   = Asym(MSK2,18); % 12
DK(3,3,2,2,:)   = Asym(MSK2,18); % 12

DK(1,1,2,3,:)   = Asym(MSK2,19); % 13
DK(1,1,3,2,:)   = Asym(MSK2,19); % 13
DK(3,1,1,2,:)   = Asym(MSK2,19); % 13
DK(2,1,1,3,:)   = Asym(MSK2,19); % 13
DK(2,3,1,1,:)   = Asym(MSK2,19); % 13
DK(3,2,1,1,:)   = Asym(MSK2,19); % 13
DK(1,3,1,2,:)   = Asym(MSK2,19); % 13
DK(1,2,1,3,:)   = Asym(MSK2,19); % 13
DK(3,1,2,1,:)   = Asym(MSK2,19); % 13
DK(2,1,3,1,:)   = Asym(MSK2,19); % 13
DK(1,2,3,1,:)   = Asym(MSK2,19); % 13
DK(1,3,2,1,:)   = Asym(MSK2,19); % 13

DK(2,2,1,3,:)   = Asym(MSK2,20); % 14
DK(2,2,3,1,:)   = Asym(MSK2,20); % 14
DK(3,2,2,1,:)   = Asym(MSK2,20); % 14
DK(1,2,2,3,:)   = Asym(MSK2,20); % 14
DK(1,3,2,2,:)   = Asym(MSK2,20); % 14
DK(3,1,2,2,:)   = Asym(MSK2,20); % 14
DK(3,2,1,2,:)   = Asym(MSK2,20); % 14
DK(2,3,2,1,:)   = Asym(MSK2,20); % 14
DK(1,2,3,2,:)   = Asym(MSK2,20); % 14
DK(2,1,2,3,:)   = Asym(MSK2,20); % 14
DK(2,1,3,2,:)   = Asym(MSK2,20); % 14
DK(2,3,1,2,:)   = Asym(MSK2,20); % 14

DK(3,3,1,2,:)   = Asym(MSK2,21); % 15
DK(3,3,2,1,:)   = Asym(MSK2,21); % 15
DK(2,3,3,1,:)   = Asym(MSK2,21); % 15
DK(1,3,3,2,:)   = Asym(MSK2,21); % 15
DK(1,2,3,3,:)   = Asym(MSK2,21); % 15
DK(2,1,3,3,:)   = Asym(MSK2,21); % 15
DK(1,3,2,3,:)   = Asym(MSK2,21); % 15
DK(2,3,1,3,:)   = Asym(MSK2,21); % 15
DK(3,2,3,1,:)   = Asym(MSK2,21); % 15
DK(3,1,3,2,:)   = Asym(MSK2,21); % 15
DK(3,1,2,3,:)   = Asym(MSK2,21); % 15
DK(3,2,1,3,:)   = Asym(MSK2,21); % 15