function out=RICSfilter(frame, filterwindow, startframe, endframe)

frame2=zeros(size(frame,1), size(frame,2), size(frame,3)+1-filterwindow);


for i=filterwindow/2:endframe+1-startframe-filterwindow/2
blockmean =mean(frame(:,:,i-filterwindow/2+1:i+filterwindow/2),3);
frame2(:,:,i-filterwindow/2+1)=frame(:,:,i)-blockmean;
end

out=frame2;