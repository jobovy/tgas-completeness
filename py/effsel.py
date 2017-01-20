import numpy
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

                            