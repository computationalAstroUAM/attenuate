# Esta función sirve para abrir un fichero de datos tipo .fits,
# toma como argumentos el nombre del fichero y la columna deseada
# y devuelve un array con los valores de dicha columna

from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table
from astropy.io import fits

def IO_fits(file,column):

    # file es el nombre del fichero y c es el nombre de la columna
    # deben escribirse ente '...' ej: IO('SAGE.fits','LOII')
    
    Data = get_pkg_data_filename(file)
    
    datos = Table.read(Data, hdu=1)
    
    # La función devuelve un array v con la columna deseada
    
    v = datos[column]

    return v

# Recuerda que si quieres guardarlo en tu Varible explorer,
# debes hacer v = IO_fits('SAGE.fits','LOII')