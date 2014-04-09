"""Simple tools to style a talk presented from an IPython Notebook.

Author: Fernando Perez <fernando.perez@berkeley.edu>
"""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

# stdlib
import os

# Third party
import matplotlib.pyplot as plt

# From IPython
from IPython.display import (HTML, Image, display, YouTubeVideo, Math,
                             clear_output)

#from IPython.html.widgets.interact import (interact,
#                                           interactive as _interactive)

#-----------------------------------------------------------------------------
# Functions and classes
#-----------------------------------------------------------------------------

def interactive(**kw):
    def deco(f):
        display(_interactive(f, **kw))
    return deco


def prefix(url):
    prefix = '' if url.startswith('http') else 'http://'
    return prefix + url


def simple_link(url, name=None):
    name = url if name is None else name
    url = prefix(url)
    return '<a href="%s">%s</a>' % (url, name)


def html_link(url, name=None):
    return HTML(simple_link(url, name))


# Utility functions
def website(url, name='auto', width=800, height=450):
    html = []
    name = url if name == 'auto' else name
    if name:
        html.extend(['<div sytle="margin-bottom:10px">',
                     simple_link(url, name),
                     '</div>'] )

    html.append('<iframe src="%s"  width="%s" height="%s"></iframe>' % 
                (prefix(url), width, height))
    return HTML('\n'.join(html))


def nbviewer(url, name=None, width=800, height=450):
    return website('nbviewer.ipython.org/url/' + url, name, width, height)


def Audio(fname):
    """Provide a player widget for an audio file.
    
    Parameters
    ==========
    fname : string
      Filename to be played.
      
    Warning
    =======
    
    Browsers cache audio very aggressively. If you change an
    audio file on disk and are trying to listen to the  new version, you 
    may want to 
    """
    from IPython.display import HTML, display
    
    # Find out file extension and deduce MIME type for audio format
    ext = os.path.splitext(fname)[1].replace('.', '').lower()
    mimetype = 'audio/' + ('mpeg' if ext == 'mp3' else ext)
    
    tpl = """<p>{fname}:</p>
<audio controls>
    <source src="files/{fname}" type="{mimetype}">

Your browser does not support the Audio element; you can play 
<a href="files/{fname}">this file</a> manually.
</audio>
"""
    return HTML(tpl.format(**locals()))


def Video(fname):
    from IPython.display import HTML, display
    video = open(fname, "rb").read()
    video_encoded = video.encode("base64")
    video_tag = '''<video controls alt="test" src="data:video/x-m4v;base64,{0}">
    </video>'''.format(video_encoded)
    return HTML(data=video_tag)    

    
def plot_audio(fname):
    from scipy.io import wavfile
    rate, x = wavfile.read(fname)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    ax1.plot(x); ax1.set_title('Raw audio signal')
    ax2.specgram(x); ax2.set_title('Spectrogram');
    plt.show()

    
#-----------------------------------------------------------------------------
# Load and publish CSS
#-----------------------------------------------------------------------------
if __name__ == '__main__':
    display(HTML(open('style.css').read()))
