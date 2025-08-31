function [n_curves,RA1, Dec]=seperate(ra,dec)
tol = 100;
            curve_no = 1;
            n_curves = 1;
            k = 0;
            ra_prev = ra(1);
            for i = 1:length(ra)
                if abs(ra(i) - ra_prev) > tol
                    curve_no = curve_no + 1;
                    n_curves = n_curves + 1;
                    k = 0;
                end
                k = k + 1;
                RA1{curve_no}(k) = ra(i);
                Dec{curve_no}(k) = dec(i);
                ra_prev = ra(i);
            end
end