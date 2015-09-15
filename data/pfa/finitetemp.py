from glob import glob
from math import pi, exp
import numpy as np
from scipy.integrate import quad
from scipy.special import zetac
from mpmath import polylog
from pyx import *

def extract_matsubara(line):
    '''extract Matsubara values from line

    The lines are assumed to start with "# n=". They should contain
    entries for n, value and possibly m. If m is not present, None
    is returned instead.

    '''
    assert line.startswith('# ')
    line = line[2:]
    entries = {}
    for entry in line.split(','):
        label, valstr = entry.strip().split('=')
        if label in ('m', 'n'):
            val = int(valstr)
        else:
            val = float(valstr)
        entries[label] = val
    return entries['n'], entries.get('m', None), entries['value']

def sum_over_m(m_terms):
    '''perform sum over m

    '''
    m_terms_sum = {}
    for k, v in m_terms.items():
        n, m = k
        if m == 0:
            v = 0.5*v
        m_terms_sum[n] = m_terms_sum.get(n, 0)+v
    return m_terms_sum

def check_matsubaras(m_terms, m_sums):
    '''check consistency of Matsubara sums

    '''
    m_terms_sum = sum_over_m(m_terms)
    for n in range(max(m_sums)+1):
        assert(abs(1-m_terms_sum[n]/m_sums[n]) < 1e-14)

def check_freeenergy(F, T, m_sums):
    '''check consistency of free energy

    The check is based on eq. (5.44) in M. Hartmann, master thesis

    '''
    f_from_sum = T/pi*(sum(m_sums.values())-0.5*m_sums[0])
    assert(abs(1-f_from_sum/F) < 1e-9)

def slurp(filename, mterms=False):
    '''read sphere-plane data from file

    A couple of consistency checks are done.

    Parameters:
    -----------
    filename : name of file to be read
    mterms : boolean indicating whether individual m terms should
             be returned (default: False)

    Returns:
    --------
    LbyR : float, surface-surface distance over sphere radius
    T : float, temperature
    F : float, free energy
    m_terms : dictionary containing the individual Matsubara terms
              with key (n, m)
    m_sums : dictionary where the sums over m have been carried out
             and the key is the Matsubara index n
    
    All values obey the scaling in eq. (5.41) in [1]

    References:
    -----------
    [1] M. Hartmann, master thesis (Univ. Augsburg, 2014),
        http://dx.doi.org/10.5281/zenodo.12476

    '''
    m_terms = {}
    m_sums = {}
    with open(filename, 'r') as fh:
        for line in fh:
            line = line.strip()
            empty = line == ''
            comment = line.startswith('#')
            values = line.startswith('# n=')
            if values:
                n, m, val = extract_matsubara(line)
                if m is None:
                    m_sums[n] = val
                else:
                    m_terms[(n, m)] = val
            elif not(empty or comment):
                vals = line.split(',')
                LbyR, T, F = map(float, vals[:3])
    if m_sums:
        check_matsubaras(m_terms, m_sums)
    else:
        m_sums = sum_over_m(m_terms)
    check_freeenergy(F, T, m_sums)
    if mterms:
        return LbyR, T, F, m_terms, m_sums
    else:
        return LbyR, T, F, m_sums

def fpfa(x):
    '''return free energy in PFA at T=0

    Parameters:
    -----------
    x : float, L/R

    Returns:
    --------
    f : float, free energy

    References:
    -----------
    [1] Eq. (7.7) in M. Hartmann, master thesis

    '''
    return -pi**3/720*(1+2*x)/(x**3+x**2)

def integrand(x, LbyR, T):
    """alternative integrand for PFA and finite temperature

    """
    n = 1
    sum = 0.5*(1+zetac(3))
    alpha1 = 2*T*x*LbyR/(1+LbyR)
    while True:
        alpha = alpha1*n
        arg = exp(-alpha)
        value = arg/(1-arg)*(1/n**3+alpha1/(n**2*(1-arg)))
        sum += value
        if value/sum < 1e-15:
            return sum/(x**2*LbyR)
        n += 1

def pfa(LbyR, T):
    """Calculate free energy according to PFA for perfect reflectors for L/R
    and T

    """
    I = quad(integrand, 1, 1+1/LbyR, args=(LbyR, T), epsrel=1e-12)
    return -T/(4*pi)*I[0]


def finitetemp(temperature, data, factor):
    '''determine a finite temperature free energy normalized by the
    zero temperature free energy in PFA

    Parameters:
    -----------
    temperature : base temperature relative to sphere radius
    data : list containing tuples (L/R, array of Matsubara terms)
    factor : float by which the base temperature is multiplied

    Returns:
    --------
    result : list containing tuples (L/R, free energy)

    '''
    return [(x, factor*temperature*(1+x)/pi*(np.sum(y[::factor])-0.5*y[0])
                / pfa(x, factor*temperature*(1+x)))
            for x, y in data]

if __name__ == '__main__':
    temperature = 0.8
    datadict = {}
    for filename in glob('slurm-*.out'):
        LbyR, T, F, m_sums = slurp(filename)
        if abs(temperature-T/(1+LbyR)) < 1e-10:
            datadict[LbyR] = m_sums
    print('all data have been read')
    data = sorted(datadict.items())
    data = [(x, np.array(sorted(d.items()))[:, 1]) for x, d in data]
    text.set(text.LatexRunner)
    g = graph.graphxy(width=10,
                      x=graph.axis.log(min=0.006, max=1, title='$L/R$'),
                      y=graph.axis.lin(min=0.5, max=1,
                          title='$\mathcal{F}/\mathcal{F}_\mathrm{PFA}(T)$'),
                      key=graph.key.key(pos="bl", dist=0.1,
                                        hdist=0.3*unit.v_cm,
                                        vdist=0.2*unit.v_cm)
                      )
    title_str = '$2\pi k_\mathrm{B}TR/\hbar c = %5.1f$'
    g.plot([graph.data.points(finitetemp(temperature, data, factor),
                    x=1, y=2, title=title_str % (temperature*factor))
                for factor in (1, 2, 5, 10, 20, 50, 100, 200, 500)],
           [graph.style.line([style.linestyle.solid,
                              color.gradient.ReverseRainbow])])
    g.writePDFfile()
