# Some utilities for the effective selection and stellar densities
import os, os.path
import numpy
from scipy import interpolate
from isodist import PadovaIsochrone
def load_spectral_types(filename=os.path.join(\
        os.path.dirname(os.path.realpath(__file__)),
        'EEM_dwarf_UBVIJHK_colors_Teff.txt')):
    names= ['SpT','Teff','logT','BCv','Mv','logL','B-V','Bt-Vt','U-B',
            'V-Rc','V-Ic','V-Ks','J-H','H-K','Ks-W1','Msun','logAge',
            'b-y','SpT_repeat','M_J','M_Ks','Mbol','i-z','z-Y','W1-W2']
    # Construct dtype
    dtype= ['S5']
    dtype.extend([numpy.float64 for ii in range(17)])
    dtype.append('S5')
    dtype.extend([numpy.float64 for ii in range(6)])
    final_dtype= []
    for name,dt in zip(names,dtype):
        final_dtype.append((name,dt))
    dtype= numpy.dtype(final_dtype)
    data= numpy.genfromtxt(filename,dtype=dtype,comments='#',
                           missing_values='...',filling_values=numpy.nan,
                           names=names,skip_footer=350)
    return data

sp= load_spectral_types()
sp_indx= numpy.array([(not 'O' in s)*(not 'L' in s)*(not 'T' in s)\
                          *(not 'Y' in s)\
                          *(not '.5V' in s)*(s != 'A8V') for s in sp['SpT']],
                     dtype='bool') #A8V has the same MJ as A9V
# Cut out the small part where the color decreases
sp_indx*= (numpy.roll((sp['JH']+sp['HK']),1)-(sp['JH']+sp['HK'])) <= 0.
ip_eems= interpolate.UnivariateSpline((sp['JH']+sp['HK'])[sp_indx],
                                      sp['M_J'][sp_indx],k=3,s=1.)
def main_sequence_cut_r(jk,low=False,tight=False):
    """Main-sequence cut, based on MJ, high as in low"""
    j_locus= ip_eems(jk)
    if low and tight:
        dj= 0.2-0.1*(j_locus-5.)
        dj[dj < 0.2]= 0.2
    elif low:
        dj= 0.2-0.25*(j_locus-5.)
        dj[dj < 0.2]= 0.2
    elif tight:
        djk= -(jk-0.6)/20.
        djk[djk>0.]= 0.
        j_locus= ip_eems(jk+djk)
        dj= 0.2-0.5*(j_locus-5.)
        dj[dj < 0.2]= 0.2
        dj[dj > 1.5]= 1.5
        dj*= -1.
    else:
        djk= -(jk-0.6)/5.
        djk[djk>0.]= 0.
        j_locus= ip_eems(jk+djk)
        dj= 1.-.8*(j_locus-5.)
        dj[dj < 0.2]= 0.2
        dj[dj > 2.5]= 2.5
        dj*= -1.
    return j_locus+dj

# Giant sequence from isochrone
iso= PadovaIsochrone(type='2mass-spitzer-wise',Z=0.017,parsec=True)
p= iso(logage=9.8,Z=0.017)
p_indx= (p['M_ini'] > 1.14)*(p['M_ini'] < 1.1845)\
    *(numpy.roll(p['J']-p['Ks'],-1) > p['J']-p['Ks'])
ip_giant= interpolate.UnivariateSpline((p['J']-p['Ks'])[p_indx],
                                        p['J'][p_indx],k=3,s=.01)
def giant_sequence_cut(jk,low=False,tight=False):
    out= numpy.empty_like(jk)
    if low and tight:
        jk_indx= jk < 0.45
        out[jk_indx]= main_sequence_cut_r(jk[jk_indx],low=False,tight=False)
        out[True-jk_indx]= ip_giant(jk[True-jk_indx]-0.05)+0.1
    elif low:
        jk_indx= jk < 0.45
        out[jk_indx]= main_sequence_cut_r(jk[jk_indx],low=False,tight=False)
        out[True-jk_indx]= ip_giant(jk[True-jk_indx]-0.1)+0.35
    elif tight:
        out= ip_giant(jk+0.1+0.05*numpy.exp(-(jk-0.58)**2./0.01))
    else:
        out= ip_giant(jk+0.25)
    return out
    
