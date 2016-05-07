import sys
import numpy as np

def htmltab(plots, file=None, xtitle=None, ytitle=None, size=100, header=None) :
    """
    Given 2D list of plots, generates HTML that displays these in a table, with optional
    titling of rows and columns and top header

    Args:
         plots : 2D list of image files to display in table
         [file=] : name of file to write HTML to
         [xtitle=] : list of titles along the horizontal axis (first row)
         [ytitle=] : list of titles along the vertical axis (first column)  
         [header=] : header to write above the table
         [size=]   : width of each plot
    """
    p=np.array(plots)
    nx=p.shape[1]
    ny=p.shape[0]

    if file is not None : 
        f = open(file,'w') 
    else : 
        f=sys.stdout
    f.write('<HTML><BODY>\n')
    if header is not None :
        f.write(header+'<p>\n')
    f.write('<TABLE BORDER=2>\n')
    if xtitle is not None :
        f.write('<TR>\n')
        if ytitle is not None :
            f.write('<TD>\n')
        for ix in range(nx) :
            f.write('<TD>'+xtitle[ix]+'\n')
    for iy in range(ny) :
        f.write('<TR>\n')
        if ytitle is not None :
            f.write('<TD>'+ytitle[iy]+'\n')
        for ix in range(nx) :
            f.write('<TD>\n')
            f.write('<A HREF='+p[iy][ix]+'>'+
                    '<IMG SRC='+p[iy][ix]+' WIDTH='+str(size)+'%></A>\n')
            f.write('</TD>\n')
    f.write('</TABLE>\n')
    f.write('</BODY></HTML>\n')

def tab(tab,file=None,sorttable=None) :
    """
    Makes HTML table from input structured array
    """

    if file is not None : 
        f = open(file,'w') 
    else : 
        f=sys.stdout
    f.write('<HTML>\n')
    if sorttable is not None:
        f.write('<HEAD><script type=text/javascript src='+sorttable+'></script></head>')
    f.write('<BODY>\n')
    f.write('<TABLE BORDER=2>\n')
    f.write('<TR>\n')
    for name in tab.dtype.names :
        f.write('<TD>'+name+'\n')
    for i in range(len(tab)) :
        f.write('<TR>\n')
        for name in tab.dtype.names :
            f.write('<TD>'+tab[name][i]+'\n')
    f.write('</TABLE>')
    f.write('</BODY>')
    f.write('</HTML>')
