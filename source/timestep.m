classdef timestep
    methods (Static)
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
        
        %-----------------------------------------------------------------%

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