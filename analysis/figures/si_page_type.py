from matplotlib.text import Text
import matplotlib.pyplot as plt

def fit_si_type(fig, pad=None, *, ticks=9.8, labels=10.8, letters=12.0):
    if pad is None:
        pad = plt.rcParams['savefig.pad_inches']
    fig.canvas.draw()
    axes_labels = {text for ax in fig.axes for text in (ax.xaxis.label, ax.yaxis.label)}
    axes_labels.update(fig.texts)
    panels = {text for ax in fig.axes for text in ax.texts
              if text.get_text() in tuple('abcdefgh') and text.get_weight() == 'bold'}
    sizes = {t: (letters if t in panels else labels if t in axes_labels else ticks)
             for t in fig.findobj(match=Text) if t.get_visible() and t.get_text().strip()}
    for _ in range(4):
        bounds = fig.get_tightbbox(fig.canvas.get_renderer()).padded(pad)
        scale = min(432 / (bounds.width*72), 504 / (bounds.height*72))
        for text, size in sizes.items():
            text.set_fontsize(size / scale)
        fig.canvas.draw()
