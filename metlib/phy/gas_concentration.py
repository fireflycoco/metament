from constants import *
from molecular_weight import *

def ppx_to_ugm3(spec='dry_air', conc=1.0, unit='ppb', p=p_std, T=T_room):
    return p * M[spec] * conc / (R * T) * units[unit] * 1E6

def ugm3_to_ppx(spec='dry_air', conc=1.0, unit='ppb', p=p_std, T=T_room):
    return R * T * conc / (p * M[spec]) / units[unit] * 1E-6

def mixratio_to_ugm3(spec='dry_air', conc=1.0, p=p_std, T=T_room):
    return ppx_to_ugm3(spec=spec, conc=conc, p=p, T=T, unit='one')

def ppm_to_ugm3(spec='dry_air', conc=1.0, p=p_std, T=T_room):
    return ppx_to_ugm3(spec=spec, conc=conc, p=p, T=T, unit='ppm')

def ppb_to_ugm3(spec='dry_air', conc=1.0, p=p_std, T=T_room):
    return ppx_to_ugm3(spec=spec, conc=conc, p=p, T=T, unit='ppb')
    
def ppt_to_ugm3(spec='dry_air', conc=1.0, p=p_std, T=T_room):
    return ppx_to_ugm3(spec=spec, conc=conc, p=p, T=T, unit='ppt')

def ugm3_to_mixratio(spec='dry_air', conc=1.0, p=p_std, T=T_room):
    return ugm3_to_ppx(spec=spec, conc=conc, p=p, T=T, unit='one')

def ugm3_to_ppm(spec='dry_air', conc=1.0, p=p_std, T=T_room):
    return ugm3_to_ppx(spec=spec, conc=conc, p=p, T=T, unit='ppm')

def ugm3_to_ppb(spec='dry_air', conc=1.0, p=p_std, T=T_room):
    return ugm3_to_ppx(spec=spec, conc=conc, p=p, T=T, unit='ppb')

def ugm3_to_ppt(spec='dry_air', conc=1.0, p=p_std, T=T_room):
    return ugm3_to_ppx(spec=spec, conc=conc, p=p, T=T, unit='ppt')

if __name__ == '__main__':
    print 'test ppm', ppm_to_ugm3(spec='O3', conc=0.12)
    print 'test ppb', ppb_to_ugm3(spec='O3', conc=0.12)
    print 'test ppt', ppt_to_ugm3(spec='O3', conc=0.12)
    print 'test ppm', ugm3_to_ppm(spec='O3', conc=235)
    print 'test ppb', ugm3_to_ppb(spec='O3', conc=235)
    print 'test ppt', ugm3_to_ppt(spec='O3', conc=235)
