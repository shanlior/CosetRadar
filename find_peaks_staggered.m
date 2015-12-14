function [ locs ] = find_peaks_staggered( x, g_coset )

for i=1:size(x,1)
    x_peaks(i,:) = conv(abs(x(i,:)),g_coset.h,'same');
end
    for c = 1:size(x,1)
       [~,locs(c,:)] = findpeaks(abs(x(c,:),'MINPEAKDISTANCE',ceil(length(g_coset.h)/2));
    end


end

