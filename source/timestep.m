classdef timestep
    methods (Static)
        %--- Trapezoid method
        function result = tr(v,pp,dt,fun,P)
            % Check if the optional parameter P is provided
            if nargin < 5  % P is not provided
                result = timestep.trS(v, pp, dt, fun);
            else  % P is provided
                result = timestep.trX(v, pp, dt, fun, P);
            end 
        end
        function result = trS(v,pp,dt,fun)
            a = fun(v,0); b = -1.*fun(v,1);
            result = pp + 0.5.*dt.*(2.*a-(a+b).*pp);
            result = result ./ (1 + 0.5.*dt.*(a+b));
        end
        function result = trX(v,pp,dt,fun,P)
            b = fun(0,pp,P); a = fun(1,pp,P)-b;
            result = v + 0.5.*dt.*(fun(v,pp,P)+b);
            result = result ./ (1- 0.5.*dt.*(a-b));
        end

        %--- heuns method
        function result = hn(v,pp,dt,fun,P)
            % Check if the optional parameter P is provided
            if nargin < 5  % P is not provided
                result = timestep.hnS(v, pp, dt, fun);
            else  % P is provided
                result = timestep.hnX(v, pp, dt, fun, P);
            end 
        end
        function result = hnS(v,pp,dt,fun)
            result = pp + dt.*fun(v,pp);
            result = pp + 0.5.*dt.*(fun(v,pp)+fun(v,result));
        end
        function result = hnX(v,pp,dt,fun,P)
            result = v + dt.*fun(v,pp,P);
            result = v + 0.5.*dt.*(fun(v,pp,P)+fun(result,pp,P));
        end

        %--- backward euler methods
        function result = be(v,pp,dt,fun,P)
           % Check if the optional parameter P is provided
            if nargin < 5  % P is not provided
                result = timestep.beS(v, pp, dt, fun);
            else  % P is provided
                result = timestep.beX(v, pp, dt, fun, P);
            end 
        end
        function result = beS(v,pp,dt,fun)
            result = (pp+dt.*fun(v,0))./(1+dt.*(fun(v,0)-fun(v,1)));
        end
        function result = beX(v,pp,dt,fun,P)
            result = (v +dt.*fun(0,pp,P))./(1+dt.*(fun(0,pp,P)-fun(1,pp,P)));
        end

        %--- forward euler methods
        function result = fe(v,pp,dt,fun,P)
            % Check if the optional parameter P is provided
            if nargin < 5  % P is not provided
                result = timestep.feS(v, pp, dt, fun);
            else  % P is provided
                result = timestep.feX(v, pp, dt, fun, P);
            end
        end
        function result = feS(v,pp,dt,fun)
            result = pp + dt.*fun(v,pp);
        end
        function result = feX(v,pp,dt,fun,P)
            result = v+dt.*fun(v,pp,P);
        end
        
        %--- rk4 methods
        function result = rk4(v, pp, dt, fun, P)
            % Check if the optional parameter P is provided
            if nargin < 5  % P is not provided
                result = timestep.rk4S(v, pp, dt, fun);
            else  % P is provided
                result = timestep.rk4X(v, pp, dt, fun, P);
            end
        end
        function result = rk4S(v,pp,dt,fun)
            p1 = pp;
            p2 = pp + (0.5).*dt.*fun(v,p1);
            p3 = pp + (0.5).*dt.*fun(v,p2);
            p4 = pp + (1.0).*dt.*fun(v,p3);
            result = pp + (1/6).*dt.*(fun(v,p1)+ 2.*fun(v,p2)+ ...
                                   2.*fun(v,p3)+ fun(v,p4));
        end
        function result = rk4X(v,pp,dt,fun,P)
            p1 = v;
            p2 = v + (0.5).*dt.*fun(p1,pp,P);
            p3 = v + (0.5).*dt.*fun(p2,pp,P);
            p4 = v + (1.0).*dt.*fun(p3,pp,P);
            result = v + (1/6).*dt.*(fun(p1,pp,P) + 2.*fun(p2,pp,P) + ...
                                  2.*fun(p3,pp,P) + fun(p4,pp,P));
        end
        
        %--- sbdf2 method
        function bs = sbdf2(dt,u,upre,ss,sspre,func)
            bs = (4/3).*ss + (4/3).*dt.*func(u, ss);         
            bs = bs+(-1/3).*sspre-(2/3).*dt.*func(upre,sspre);
        end
        function b = sbdf2b(P, dt, n, m, h, u, npre, mpre, hpre, upre)
        % sbdf2 Compute the Adams–Bashforth/Moulton-style combination "b"
        % using Hodgkin–Huxley conductances for current and previous states.
        %
        % Inputs:
        %   P     : struct with fields C, gk, ek, gna, ena, gl, el
        %   n,m,h : current gate variables (vectors OK)
        %   u     : current membrane potential
        %   npre,mpre,hpre,upre : previous-step gate variables and potential
        %
        % Output:
        %   b     : combined RHS term
        %
        % Notes:
        %   - All operations are element-wise; vectors/matrices are supported.
        %   - dt is passed.
                
            % --- precompute constants ---
            invC    = 1 ./ P.C;
            gk_ek   = P.gk  .* P.ek;
            gna_ena = P.gna .* P.ena;
            gl_el   = P.gl  .* P.el;
        
            % --- powers used repeatedly ---
            n4   = n   .^ 4;   m3   = m   .^ 3;
            n4p  = npre.^ 4;   m3p  = mpre.^ 3;
        
            % --- helper: returns [a1, a2] for given (n^4, m^3, h) ---
            coeffs = @(n4v, m3v, hv) deal( ...
                invC .* (gk_ek   .* n4v + gna_ena .* m3v .* hv + gl_el), ...  % a1
               -invC .* (P.gk    .* n4v + P.gna    .* m3v .* hv + P.gl) ...   % a2
            );
        
            % current and previous coefficients
            [a1,  a2 ] = coeffs(n4,  m3,  h);
            [a1p, a2p] = coeffs(n4p, m3p, hpre);
        
            % --- final combination (unchanged math) ---
            b = (4/3).*u    + (4/3).*dt.*(a2  .* u    + a1 ) ...
              - (1/3).*upre - (2/3).*dt.*(a2p .* upre + a1p);
        end
    end
end