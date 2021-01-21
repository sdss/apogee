###############################################################################
# solarabundances.py: sets of solar abundances
###############################################################################
_ASPLUND05= {1:12.00,2:10.93,3:1.05,4:1.38,5:2.70,6:8.39,7:7.78,8:8.66,9:4.56,
             10:7.84,11:6.17,12:7.53,13:6.37,14:7.51,15:5.36,16:7.14,17:5.50,
             18:6.18,19:5.08,20:6.31,21:3.05,22:4.90,23:4.00,24:5.64,25:5.39,
             26:7.45,27:4.92,28:6.23,29:4.21,30:4.60,31:2.88,32:3.58,33:2.29,
             34:3.33,35:2.56,36:3.28,37:2.60,38:2.92,39:2.21,40:2.59,41:1.42,
             42:1.92,43:0.00,44:1.84,45:1.12,46:1.69,47:0.94,48:1.77,49:1.60,
             50:2.00,51:1.00,52:2.19,53:1.51,54:2.27,55:1.07,56:2.17,57:1.13,
             58:1.58,59:0.71,60:1.45,61:0.00,62:1.01,63:0.52,64:1.12,65:0.28,
             66:1.14,67:0.51,68:0.93,69:0.00,70:1.08,71:0.06,72:0.88,73:0.13,
             74:1.11,75:0.23,76:1.45,77:1.38,78:1.64,79:1.01,80:1.13,81:0.90,
             82:2.00,83:0.65,84:0.00,85:0.00,86:0.00,87:0.00,88:0.00,89:0.00,
             90:0.06}
def asplund05():
    """
    NAME:
       asplund05
    PURPOSE:
       return a dictionary with the Asplund et al. (2005) solar abundances
    INPUT:
       (none)
    OUTPUT:
       dictionary with the Asplund et al. (2005) abundances
    HISTORY:
       2015-04-16 - Written - Bovy (IAS)
    """
    return _ASPLUND05
