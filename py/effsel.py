# Some utilities for the effective selection and stellar densities
import numpy
from scipy import interpolate
def load_spectral_types(filename='EEM_dwarf_UBVIJHK_colors_Teff.txt'):
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
        tjk= jk
        djk= -(jk-0.6)/5.
        djk[djk>0.]= 0.
        j_locus= ip_eems(jk+djk)
        dj= 1.-.8*(j_locus-5.)
        dj[dj < 0.2]= 0.2
        dj[dj > 2.5]= 2.5
        dj*= -1.
    return j_locus+dj
