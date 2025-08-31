function T1=sep(oo)

              event_no = 1;
                n_event = 1;
                n = 0;
                timed = datevec(oo(1));
                for i = 1:length(oo)
                    if abs(etime(datevec(oo(i)) , timed)) >=20
                        event_no = event_no + 1;
                        n_event= n_event + 1;
                        n = 0;
                    end
                    n = n + 1;

                    T1{event_no}(n) = oo(i);
                    timed = datevec(oo(i));
                end
            
end