%%
A = importdata('/Users/nicolas/OneDrive/Documents/Honours/Scripts/ratomatic/Rat9Scan2/20180830_PTSD_Rat9Scan2_v1_gatetags3.log');
%A = importdata('/Users/nicolas/OneDrive/Documents/Honours/Scripts/ratomatic/Rat7Scan2/20180803_PTSD_Rat7Scan1_v1_gatetags3.log');
M = Load_mPET_motion_file('/Users/nicolas/OneDrive/Documents/Honours/Scripts/ratomatic/Rat9Scan2/LORebin_Rat9Scan2.dat',1);
 
figure('name','rat9Scan2'); 
 subplot(2,1,1); plot(diff(A(:,1))), title('gate interval')
 subplot(2,1,2); plot(diff(A(:,2))), title('time interval')

 %%
gate8glitchloc = find(A(:,3) > 80); %find gate8 missing data location
gate8glitchval = A(gate8glitchloc,3); %find ms interval length
 
for i = size(gate8glitchloc,1):-1:1
    filler_size = floor(gate8glitchval(i)/33);
    filler_insert = zeros(filler_size,3); filler_insert(:,3) = 33;
    A = [A(1:gate8glitchloc(i)-filler_size,:); filler_insert; A(gate8glitchloc(i)+1:end,:)];
end
 
numtoadd = size(M,1) - size(A,1); %syncs streams from the back
A = [zeros(numtoadd,3); A];

figure('name','Rat9Scan2 Synced');
 plot(A(:,3))
 hold on
 plot(diff(M(:,2))-5)
 ylim([0 60])
 
%%
figure; plot(A(:,2:3)-A(1,2:3))
G = A(:,2:3)-A(1,2:3);
figure; plot([G(:,1)-G(:,2)])

figure; plot(diff(A(:,3)))
figure; plot(diff(A(:,2)))