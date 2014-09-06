function [reconstructed] =gauss_reconstruct(input)
a=input;
total_frames_processed=0;
total_gauss_added=0;
reconstructed=zeros(1,length(a));
for repeat=1:1
    for power=8:-0.1:5
        gauss_frame_size=round(2^power);
        numb_of_frames=length(a)/gauss_frame_size;
        total_frames_processed=total_frames_processed+round(numb_of_frames);
        for i=1:numb_of_frames
            original_frame=a(1,((i-1)*gauss_frame_size +1) : i*gauss_frame_size);
            sigma=gauss_frame_size/8;  % gauss std dev
            mu=gauss_frame_size/2;
            gauss=zeros(1,gauss_frame_size);
            for j=1:gauss_frame_size
                gauss(1,j)=(max(abs(a))/100)*exp(-(((j-mu)/(sigma*sqrt(2)))^2));
            end   %gauss created
            
            recent_error=sum(abs(reconstructed(1,((i-1)*gauss_frame_size +1) : i*gauss_frame_size)-original_frame));
            next_error=sum(abs(reconstructed(1,((i-1)*gauss_frame_size +1) : i*gauss_frame_size)+gauss-original_frame   ));
            while(recent_error>next_error)
                reconstructed(1,((i-1)*gauss_frame_size +1) : i*gauss_frame_size)=reconstructed(1,((i-1)*gauss_frame_size +1) : i*gauss_frame_size)+gauss;
                recent_error=next_error;
                next_error=sum(abs(reconstructed(1,((i-1)*gauss_frame_size +1) : i*gauss_frame_size)+gauss-original_frame   ))  ;
                total_gauss_added=total_gauss_added+1;
            end
            
            recent_error=sum(abs(reconstructed(1,((i-1)*gauss_frame_size +1) : i*gauss_frame_size)-original_frame));
            next_error=sum(abs(reconstructed(1,((i-1)*gauss_frame_size +1) : i*gauss_frame_size)-gauss-original_frame   ));
            while(recent_error>next_error)
                reconstructed(1,((i-1)*gauss_frame_size +1) : i*gauss_frame_size)=reconstructed(1,((i-1)*gauss_frame_size +1) : i*gauss_frame_size)-gauss;
                recent_error=next_error;
                next_error=sum(abs(reconstructed(1,((i-1)*gauss_frame_size +1) : i*gauss_frame_size)-gauss-original_frame   ))  ;
                total_gauss_added=total_gauss_added+1;
            end
        end
    end
end
total_frames_processed
total_gauss_added

end
