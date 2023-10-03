function [out] = NConvolved(times,timeseries,ffil,wtl,decf,wtcrt)
% convolution function for wavelet program
% specifically for applying filters

    it = 1:decf:(length(timeseries)-length(ffil)+1);

    out = zeros(length(it),1);
    
    if length(find(abs(timeseries)<50))/length(timeseries) > wtcrt
        P = 1;
        for h = it

            if ~isempty(times)
                
                out(P) = mean(times(h:h+length(ffil)-1));   % middle of analysis window
                % out(P) = times(h);                        % first time in analysis window
                
            else
                
                num = sum(ffil .* timeseries(h:h+length(ffil)-1) .* wtl(h:h+length(ffil)-1));      % multiplies and sums filter, data, and weights
                den = sum(wtl(h:h+length(ffil)-1) .* abs(ffil)) + 1e-13;
                flo = floor(sum(wtl(h:h+length(ffil)-1))/(wtcrt*length(ffil)));

                out(P) = num / den * flo;
            
            end
           
            P = P+1;
            
        end
        
    end

end



    


