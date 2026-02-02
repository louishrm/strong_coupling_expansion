# Generated automatically using the command :
# c++2py ../../c++/sc_expansion/sc_expansion.hpp -p --members_read_only -N sc_expansion -a sc_expansion -m sc_expansion_module -o sc_expansion_module --moduledoc="The sc_expansion python module" -C triqs --cxxflags="-std=c++17" --target_file_only
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "sc_expansion_module", doc = r"The sc_expansion python module", app_name = "sc_expansion")

# Imports

# Add here all includes
module.add_include("sc_expansion/sc_expansion.hpp")
module.add_include("sc_expansion/dimer_order_4.hpp")
module.add_include("sc_expansion/diagram.hpp")

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

c.add_method("""double compute_sum_diagrams(std::vector<double> taus, bool infinite_U)""", doc = r""" """)

module.add_class(c)

c = class_(
        py_type = "Order2",  # name of the python class
        c_type = "sc_expansion::order2",   # name of the C++ class
        doc = r"""A very useful and important class""",   # doc of the C++ class
        hdf5 = False
)

c.add_constructor("""(double U, double mu, double beta)""", doc = r"""""")

c.add_method("""double compute_sum_diagrams(std::vector<double> taus, bool infinite_U)""", doc = r""" """)

module.add_class(c)

c = class_(
        py_type = "DiagramP",  # name of the python class
        c_type = "sc_expansion::Diagram",   # name of the C++ class
        doc = r"""A very useful and important class""",   # doc of the C++ class
        hdf5 = False
)

c.add_constructor("""(std::vector<std::vector<int>> adjacency_matrix, double U, double beta, double mu)""", doc = r"""""")

c.add_method("""double evaluate_at_taus(std::vector<double> taus, bool infinite_U)""", doc = r""" """)

module.add_class(c)


c = class_(
        py_type = "Order6",  # name of the python class
        c_type = "sc_expansion::order6",   # name of the C++ class
        doc = r"""A very useful and important class""",   # doc of the C++ class
        hdf5 = False
)

c.add_constructor("""(double U, double mu, double beta)""", doc = r"""""")

c.add_method("""double compute_sum_diagrams(std::vector<double> taus, bool infinite_U)""", doc = r""" """)

module.add_class(c)

c = class_(
        py_type = "Cumulant",
        c_type = "sc_expansion::CumulantHelper",
        doc = r"""Calculates cumulants for specific times with spin 0""",
        hdf5 = False
)

c.add_constructor("""(double U, double beta, double mu)""", doc = r"""""")

c.add_method("""double compute(std::vector<double> taus)""", doc = r""" """)

module.add_class(c)

module.generate_code()
