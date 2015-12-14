function [targets] = classic_processing(g, x, sample_SubNyquist_factor, pulse_SubNyquist_factor)

if nargin == 2
    sample_SubNyquist_factor = 1;
    pulse_SubNyquist_factor = 1;
end
% decimate to Nyquist rate
L = round(g.Fs/(2*g.h_BW)* sample_SubNyquist_factor);
x_Nyquist = zeros(ceil(size(x,1)/L),round(g.P / pulse_SubNyquist_factor));
for p=1:size(x_Nyquist, 2)
    x_Nyquist(:,p) = decimate(x(:,p),L);
end
h_Nyquist = decimate(g.h,L);
if g.mf.use_windows
    h_Nyquist = h_Nyquist .* chebwin(length(h_Nyquist), 50);
end

% Create Time-Frequency map
map = xcorr2(x_Nyquist, h_Nyquist);
map = map(length(h_Nyquist):end-length(h_Nyquist),:);
if g.mf.use_windows
    map = map .* repmat(chebwin(size(map,2),50).', size(map,1), 1);
end
DopplerOversamplingRatio = g.Nyquist.df / g.CS.delta_f;
map = fft(map, DopplerOversamplingRatio*g.P, 2);
map_power = abs(map).^2;
if g.mf.debug_plot
    figure;
    t_axis = (0:(size(x_Nyquist,1)-1) - length(h_Nyquist)) * g.Nyquist.dt * sample_SubNyquist_factor;
    f_axis = (0:(DopplerOversamplingRatio*g.P-1)) * g.Nyquist.df / DopplerOversamplingRatio;
    imagesc(t_axis, f_axis, 10*log10(map_power.'));
    colormap(gray);
    ca = caxis;
    caxis([ca(2)-40 ca(2)]);
    colorbar;
    set(gca,'YDir','normal');
    hold on;
    title(sprintf('Noise floor = %.1f [dB]', 10*log10(median(map_power(:))/0.7)));
end

% Peak detection
[~,map_power_indexes] = sort(map_power(:),'descend');
targets.a = zeros(g.L, 1);
targets.t = zeros(g.L, 1);
targets.f = zeros(g.L, 1);
l = 0;
index_count = 1;
local_max = ordfilt2(map_power, 9, ones(3));
Q = [.5 -1 .5 ; -.5 0 .5 ; 0 1 0]; % for parabola interpolation
while l < g.L && index_count <= length(map_power_indexes)
    [t_index,f_index] = ind2sub(size(map_power), map_power_indexes(index_count));
    if local_max(t_index, f_index) == map_power(t_index, f_index)
        l = l + 1;
        targets.a(l) = map(t_index, f_index);
        time_bias_fix = 0.4;
        
        % time intepolation
        t_index_fix = 0;
        if t_index > 1 && t_index < size(map_power, 1)
            t_power = map_power(t_index+(-1:1), f_index);
            a = Q * t_power;
            assert(a(1) < 0);
            t_index_fix = -a(2) / 2 / a(1);
        end
        targets.t(l) = ((t_index - 1 + t_index_fix) * sample_SubNyquist_factor + time_bias_fix) * g.Nyquist.dt;
        assert(targets.t(l)>=0 && targets.t(l)<=g.tau);
        
        % frequency interpolation
        f_index_fix = 0;
        if f_index > 1 && f_index < size(map_power, 2)
            f_power = map_power(t_index, f_index+(-1:1));
            a = Q * f_power.';
            assert(a(1) < 0);
            f_index_fix = -a(2) / 2 / a(1);
        end
        targets.f(l) = (f_index - 1 + f_index_fix) * g.Nyquist.df / DopplerOversamplingRatio;
        assert(targets.f(l)>=0 && targets.f(l)<=1/g.tau);
        
        if g.mf.debug_plot
            h = plot(targets.t(l), targets.f(l), 'ro');
            set(h,'MarkerSize',10);
        end
    end
    index_count = index_count + 1;
end
