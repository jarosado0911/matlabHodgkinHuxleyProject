classdef gates
    methods (Static)
        
        % this is the reaction term in the HH model
        % v  = input voltage vector
        % mm = input state m for potassium channels
        % hh = input state h for potassium channels
        % nn = input state n for sodium channels
        function result = X(v,ss,P)
            mm = ss(:,1); nn = ss(:,2); hh = ss(:,3);
            a = (-1/P.C) * (P.gk .* nn.^4 + ...
                            P.gna .* mm.^3 .* hh + ...
                            P.gl);
            
            b = (1/P.C) * (P.gk .* nn.^4 .* P.ek + ...
                           P.gna .* mm.^3 .* hh .* P.ena + ...
                           P.gl .* P.el);
            result = a.*v +b;
        end


        function yout = F(s,a,b)
        yout = a.*(1-s)-b.*s;
        end

        function out=an(vin)
        vin = vin.*1e3;
        out=(-0.032).*(vin-15)./(exp(-1.*(vin-15)./5)-1);
        out = out*1e3;
        end
        
        function out=bn(vin)
        vin = vin.*1e3;
        out=(0.5).*exp(-1.*(vin-10)./40);
        out = out*1e3;
        end
        
        function out=am(vin)
        vin = vin.*1e3;
        out=(-0.32).*(vin-13)./(exp(-1.*(vin-13)./4)-1);
        out = out*1e3;
        end
        
        function out=bm(vin)
        vin = vin.*1e3;
        out=(0.28).*(vin-40)./(exp((vin-40)./5)-1);
        out = out*1e3;
        end
        
        function out=ah(vin)
        vin = vin.*1e3;
        out=(0.128).*exp(-1.*(vin-17)./18);
        out = out*1e3;
        end
        
        function out=bh(vin)
        vin = vin.*1e3;
        out=4./(exp((40-vin)./5)+1);
        out = out*1e3;
        end
    end
end