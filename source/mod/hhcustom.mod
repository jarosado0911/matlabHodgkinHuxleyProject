TITLE Custom Hodgkin-Huxley channel matching gates.m rate equations
: alpha/beta rate functions here are copied from source/gates.m (an,bn,am,bm,ah,bh)
: so that runsim_yale_neuron.m simulates the exact same ionic kinetics as the
: project's own SBDF2/Strang MATLAB solvers - only the numerical integration
: scheme differs (NEURON's built-in "hh" mechanism uses different, classic
: squid-axon rate constants and is not a fair cross-check of this project).

NEURON {
    SUFFIX hhcustom
    USEION na READ ena WRITE ina
    USEION k READ ek WRITE ik
    NONSPECIFIC_CURRENT il
    RANGE gnabar, gkbar, gl, el
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S)  = (siemens)
}

PARAMETER {
    gnabar = 0.05  (S/cm2)
    gkbar  = 0.005 (S/cm2)
    gl     = 0     (S/cm2)
    el     = -70   (mV)
}

STATE { m h n }

ASSIGNED {
    v   (mV)
    ena (mV)
    ek  (mV)
    ina (mA/cm2)
    ik  (mA/cm2)
    il  (mA/cm2)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ina = gnabar*m*m*m*h*(v-ena)
    ik  = gkbar*n*n*n*n*(v-ek)
    il  = gl*(v-el)
}

DERIVATIVE states {
    m' = alpham(v)*(1-m) - betam(v)*m
    h' = alphah(v)*(1-h) - betah(v)*h
    n' = alphan(v)*(1-n) - betan(v)*n
}

INITIAL {
    m = alpham(v)/(alpham(v)+betam(v))
    h = alphah(v)/(alphah(v)+betah(v))
    n = alphan(v)/(alphan(v)+betan(v))
}

FUNCTION alpham(v(mV)) (/ms) {
    alpham = (-0.32)*vtrap(v-13, -4)
}
FUNCTION betam(v(mV)) (/ms) {
    betam = (0.28)*vtrap(v-40, 5)
}
FUNCTION alphah(v(mV)) (/ms) {
    alphah = 0.128*exp(-(v-17)/18)
}
FUNCTION betah(v(mV)) (/ms) {
    betah = 4/(exp((40-v)/5)+1)
}
FUNCTION alphan(v(mV)) (/ms) {
    alphan = (-0.032)*vtrap(v-15, -5)
}
FUNCTION betan(v(mV)) (/ms) {
    betan = 0.5*exp(-(v-10)/40)
}

FUNCTION vtrap(x, y) {
    : removable-singularity-safe form of x/(exp(x/y)-1), used by alpham/alphan
    if (fabs(x/y) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    } else {
        vtrap = x/(exp(x/y) - 1)
    }
}
