from tools import plots

def chip(a,row=150,ax=None) :
    """ Routine to plot 3 chips in 3 panels
    """
    if ax is None : fig,ax=plots.multi(1,3,hspace=0.3)
    plots.plotl(ax[0],a['a'][4].data[row,:],a['a'][1].data[row,:])
    plots.plotl(ax[1],a['b'][4].data[row,:],a['b'][1].data[row,:])
    plots.plotl(ax[2],a['c'][4].data[row,:],a['c'][1].data[row,:])

    try : return fig,ax
    except : return
