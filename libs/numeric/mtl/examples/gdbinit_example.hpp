python
import sys

# STL support and alike here

sys.path.insert(0, '/home/username/tools/gdb_printers/python')
from mtl.printers import register_mtl_printers
register_mtl_printers (None)

end
