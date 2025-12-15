# Generated automatically using the command :
# c++2py ../../c++/sc_expansion/sc_expansion.hpp -p --members_read_only -N sc_expansion -a sc_expansion -m sc_expansion_module -o sc_expansion_module --moduledoc="The sc_expansion python module" -C triqs --cxxflags="-std=c++17" --target_file_only
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "sc_expansion_module", doc = r"The sc_expansion python module", app_name = "sc_expansion")

# Imports

# Add here all includes
module.add_include("sc_expansion/sc_expansion.hpp")
module.add_include("sc_expansion/dimer_order_4.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/string.hpp>

using namespace sc_expansion;
""")

# The class toto
c = class_(
        py_type = "Order4",  # name of the python class
        c_type = "sc_expansion::order4",   # name of the C++ class
        doc = r"""A very useful and important class""",   # doc of the C++ class
        hdf5 = False
)

c.add_constructor("""(double U, double mu, double beta)""", doc = r"""""")

c.add_method("""double compute_sum_diagrams(std::vector<double> taus)""", doc = r""" """)

module.add_class(c)

c = class_(
        py_type = "Order2",  # name of the python class
        c_type = "sc_expansion::order2",   # name of the C++ class
        doc = r"""A very useful and important class""",   # doc of the C++ class
        hdf5 = False
)

c.add_constructor("""(double U, double mu, double beta)""", doc = r"""""")

c.add_method("""double compute_sum_diagrams(std::vector<double> taus)""", doc = r""" """)


module.generate_code()
