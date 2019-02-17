% This function takes in input a vector of numbers, and returns the subvector of peaks from it. A peak is defined as an 
% element which is preceeded and followed by lower elements than itself.

function peaks = extract_peaks(y)

n=length(y);
peaks=zeros(ceil(n/2),1);
j=0;
for i=2:n-1
    m = max(y(i-1), y(i+1));
    if y(i)>m
        j=j+1;
        peaks(j)=y(i);
    end
end
peaks = peaks(1:j);

end
