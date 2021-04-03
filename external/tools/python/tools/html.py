import sys
import numpy as np

def getfile(file=None,access='w') :
    if file is not None : 
        f = open(file,'w') 
    else : 
        f=sys.stdout
    return f

def head(file=None,sorttable=None) :
    f=getfile(file=file)
    f.write('<HTML>\n')
    if sorttable is not None:
        f.write('<HEAD><script type=text/javascript src='+sorttable+'></script></head>')
    f.write('<BODY>\n')
    return f

def tail(f) :
    f.write('</BODY>\n')
    f.write('</HTML>\n')
    if f is not sys.stdout : f.close()

def table(data, xtitle=None, ytitle=None, size=100, plots=True, formstr=':8.3f',border=2) :
    """ Returns HTML for a table of input figures
    """
    p=np.array(data)
    nx=p.shape[1]
    ny=p.shape[0]
    out='<TABLE BORDER={:d}>\n'.format(border)
    if xtitle is not None :
        out+='<TR>\n'
        if ytitle is not None :
            out+='<TD>\n'
        for ix in range(nx) :
            out+='<TD>'+xtitle[ix]+'\n'
    for iy in range(ny) :
        out+='<TR>\n'
        if ytitle is not None :
            out+='<TD>'+ytitle[iy]+'\n'
        for ix in range(nx) :
            out+='<TD>\n'
            if plots :
                out+=('<A HREF='+p[iy][ix]+'>'+
                    '<IMG SRC='+p[iy][ix]+' WIDTH='+str(size)+'%></A>\n')
            else :
                out+=('{'+formstr+'}').format(p[iy][ix])
            out+='</TD>\n'
    out+='</TABLE>\n'
    return out


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
    #p=np.array(plots)
    #nx=p.shape[1]
    #ny=p.shape[0]
    ny=len(plots)
    nx=0
    for row in plots : nx=np.max([len(row),nx])
    p=plots

    f=head(file=file)
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
        if len(p[iy]) == 1 and nx>1 :
            f.write('<TD COLSPAN={:d} bgcolor=lightblue>\n'.format(nx+1))
            f.write(p[iy][0])
            f.write('</TD>\n')
        else :
            if ytitle is not None :
                f.write('<TD>'+ytitle[iy]+'\n')
            for ix in range(nx) :
                f.write('<TD>\n')
                f.write('<A HREF='+p[iy][ix]+'>'+
                        '<IMG SRC='+p[iy][ix]+' WIDTH='+str(size)+'%></A>\n')
                f.write('</TD>\n')
    f.write('</TABLE>\n')
    tail(f)

def tab(tab,file=sys.stdout,sortable=False) :
    """
    Makes HTML table from input structured array
    """

    if sortable :
        file.write('Click on column headings to sort<br>')
        file.write('<TABLE BORDER=2 CLASS=sortable>\n')
    else :
        file.write('<TABLE BORDER=2>\n')
    file.write('<TR>\n')
    for name in tab.dtype.names :
        file.write('<TD>'+name+'\n')
    for i in range(len(tab)) :
        file.write('<TR>\n')
        for name in tab.dtype.names :
            file.write('<TD>'+str(tab[name][i])+'\n')
    file.write('</TABLE>')
