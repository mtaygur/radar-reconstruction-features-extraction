function [feature_vector] =feat_extract(input)

%close all;
a=input;
threshold=mean(sqrt(a.^2));
a_clipped=zeros(1,length(a));
for i=1:length(a)
    if a(i)>threshold
        a_clipped(i)=a(i)-threshold;
    elseif a(i)<-threshold
        a_clipped(i)=a(i)+threshold;
    else
    end
end

i=2;
feature_matrix_index=0;
feature_matrix=zeros(1,6);
figure;plot(a);hold on;plot([1,length(a)],[threshold threshold],'g');plot([1,length(a)],[-threshold -threshold],'g');hold on
%figure(2);plot(a_clipped);hold on;
while (i<length(a))
    if(abs(a_clipped(i))>eps)
        j=0;
        no_of_zeros=0;
        while (i+j+10)<length(a)  &&  any(a_clipped((i+j):(i+j+2))*a_clipped(i)>0)
            j=j+1;
            if abs(a_clipped(i+j))==0
                no_of_zeros=no_of_zeros+1;
            end
        end
        if j>9 && no_of_zeros<(j/4) %minimum width of peaks
            %figure(2);plot(i,a_clipped(i),'ko');hold on;
            %figure(2);plot(i+j+2,a_clipped(i+j+2),'g*');hold on;
            [val,indx]=max(abs(a_clipped(i:(i+j))));
            if val>(max(a_clipped)-min(a_clipped))*0.1
                feature_matrix_index=feature_matrix_index+1;
                %figure(2);plot(indx+i-1,a_clipped(indx+i-1),'r*');hold on;
                feature_matrix(feature_matrix_index,1)=indx+i-1;
                feature_matrix(feature_matrix_index,6)=a(indx+i-1);
            end
            i=i+j+1;
        end
    end
    i=i+1;
end


%if feature_matrix_index>7     % reduce peak number to 7
%   [~,ix]=sort(abs(feature_matrix(:,6)),'descend');
% feature_matrix=feature_matrix(ix(1:7),:);
%end
feature_matrix_index=length(feature_matrix(:,1));
if (feature_matrix(1,1)>eps)
    for i=1:feature_matrix_index
        j=0;
        while (feature_matrix(i,1)-j)<length(a) && (feature_matrix(i,1)-j-1)>0 && a(feature_matrix(i,1)-j)*a(feature_matrix(i,1))>eps  && a(feature_matrix(i,1)-j-1)*a(feature_matrix(i,1))>eps && (feature_matrix(i,1)-j-2)>0 
            j=j+1;   %find left-hand width
        end         
        m= a(feature_matrix(i,1)-j)- a(feature_matrix(i,1)-j-1);
        b=a(feature_matrix(i,1)-j-1)-m*(feature_matrix(i,1)-j-1);       
        j=(feature_matrix(i,1))+b/m;
        feature_matrix(i,2)=j;
        feature_matrix(i,4)=feature_matrix(i,6)/(j);  %find tangent
        j=0;
        while (feature_matrix(i,1)+j+1)<length(a) && (feature_matrix(i,1)+j)>0  && (a(feature_matrix(i,1)+j)*a(feature_matrix(i,1))>eps  && a(feature_matrix(i,1)+j+1)*a(feature_matrix(i,1))>eps) && (feature_matrix(i,1)+j+1)<length(a)
            j=j+1;  %find right-hand width
        end       
        m= a(feature_matrix(i,1)+j+1)- a(feature_matrix(i,1)+j);
        b=a(feature_matrix(i,1)+j)-m*(feature_matrix(i,1)+j);       
        j=-(feature_matrix(i,1))-b/m;
        feature_matrix(i,3)=j;
        feature_matrix(i,5)=feature_matrix(i,6)/(j); %find tangent
    end    
    
    [~,indx]=max(abs(feature_matrix(:,6)));
    
    
    
    if indx >= 3
        if indx>3
            feature_matrix(1:indx-3,:)=zeros(indx-3,6);
        end
        for k=indx:-1:indx-1
            if (feature_matrix(k,1)-feature_matrix(k-1,1))>1.5*(feature_matrix(k,2)+feature_matrix(k-1,3))
                feature_matrix(1:k-1,:)=zeros(k-1,6);
                break;
            end
        end
    elseif indx==2
        if (feature_matrix(indx,1)-feature_matrix(indx-1,1))>1.5*(feature_matrix(indx,2)+feature_matrix(indx-1,3))
            feature_matrix(indx-1,:)=zeros(1,6);
        end
    end
    
    if indx+2 <= length(feature_matrix(:,1))
        if indx+2 < length(feature_matrix(:,1))
           feature_matrix(indx+3:end,:)=zeros(length(feature_matrix(indx+3:end,1)),6);
        end
        for k=indx:indx+1
            if (feature_matrix(k+1,1)-feature_matrix(k,1))>1.5*(feature_matrix(k+1,2)+feature_matrix(k,3))
                feature_matrix(k+1:end,:)=zeros(length(feature_matrix(k+1:end,1)),6);
                break;
            end
        end
    elseif indx==length(feature_matrix(:,1))-1
        if (feature_matrix(indx+1,1)-feature_matrix(indx,1))>1.5*(feature_matrix(indx+1,2)+feature_matrix(indx,3))
                feature_matrix(indx+1,:)=zeros(1,6);
        end        
    end
    feature_matrix
    i=1;
    while (feature_matrix(i,1)==0)
        i=i+1;
    end
     j=length(feature_matrix(:,1));
    while (feature_matrix(j,1)==0) 
        j=j-1;
    end    
    feature_matrix=feature_matrix(i:j,:);
    
    
    %figure(1);
    plot(feature_matrix(:,1),feature_matrix(:,6),'r*');hold on;axis([0,1023,1.1*min(a),1.1*max(a)])
    
    
    feature_matrix(:,[1 6])=feature_matrix(:,[6 1]);
    feature_matrix_index=length(feature_matrix(:,1));
    peaks_distances=zeros(1,length(feature_matrix(:,1))-1);
    for i=1:length(feature_matrix(:,1))-1
        peaks_distances(1,i)=abs(feature_matrix(i,6)-feature_matrix(i+1,6));
    end
    feature_matrix=feature_matrix(:,1:5);
    feature_vector=zeros(1,5*(feature_matrix_index) + length(peaks_distances));
    for i=1:(feature_matrix_index)
        feature_vector(1,((i*5)-4):(i*5))=feature_matrix(i,1:5);
    end
    feature_vector(1,((5*i) +1):end)=peaks_distances(1:end);
else
    feature_vector=0;
end


end
