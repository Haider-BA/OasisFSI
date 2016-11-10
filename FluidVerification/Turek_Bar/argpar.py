import argparse
from argparse import RawTextHelpFormatter

def parse():
    parser = argparse.ArgumentParser(description="########################################"
    "\nImplementation of Turek CFD test cases\n######################################## \n \n"
    "-  The program automaticly stores experiment parameters and plots of lift and drag \n"
    "   in the experiment folder according to choice of experiment \n \n"
    "-  For details of numerical benchmark go to:\n   http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.550.1689&rep=rep1&type=pdf",\
     formatter_class=RawTextHelpFormatter, \
     epilog="############################################################################\n"
     "Example --> python cfd.py -case 1 -T 0.02 -dt 0.01 -v_deg 2 -p_deg 1 -solver Newton\n"
     "Example --> python cfd.py -case 2 -solver Newton2  -v_deg 2 -p_deg 1 -r  (Refines mesh one time, -rr for two etc.) \n"
     "############################################################################")
    group = parser.add_argument_group('Parameters')
    group.add_argument("-case",               help="Define which cfd case to run               --> Default=1", default="1")
    group.add_argument("-T",      type=float, help="Set end time                               --> Default=0.02", default=0.02)
    group.add_argument("-dt",     type=float, help="Set time step                              --> Default=0.01", default=0.01)
    group.add_argument("-p_deg",  type=int, help="Set degree of pressure                     --> Default=1", default=1)
    group.add_argument("-v_deg",  type=int, help="Set degree of velocity                     --> Default=2", default=2)
    group.add_argument("-theta",  type=float, help="Explicit, Implicit, Cranc-Nic (0, 1, 0.5)  --> Default=1", default=1)
    group.add_argument("-r", "--refiner", action="count", help="Mesh-refiner using built-in FEniCS method refine(Mesh)")
    group2 = parser.add_argument_group('Solvers')
    group2.add_argument("-solver", help="Newton   -- Fenics built-in module \n"
    "Newton2  -- Manuell implementation\n"
    "Default  --> Newton", default="Newton")
    parser.parse_args()
    return parser.parse_args()
