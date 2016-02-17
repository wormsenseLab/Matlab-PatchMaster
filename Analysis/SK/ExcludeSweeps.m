% ExcludeSweeps.m





figure()
for i = Files(:,1):Files(:,end)
%plot(Time,ASubtract(:,i));
%title('Current - with leak substraction')
%hold on
%plot(Time,A(:,i));
%subplot(4,2,i)
%figure()
subplot(7,5,i-Files(:,1)+1)
plot(Time,A(:,i))
hold on
plot(Time,ASubtract(:,i))
title(MeanIndentation(i))
end
