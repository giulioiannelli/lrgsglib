from IPython import get_ipython
from IPython.core.magic import register_cell_magic


@register_cell_magic
def skip_cell(line, cell):
    return


@register_cell_magic
def skip_cell_if(line, cell):
    if eval(line):
        return
    get_ipython().run_cell(cell)
