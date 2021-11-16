# Description
Calculate bond graph parameters from kinetic parameters for an individual system, or a composite system.

# Individual
Navigate to folder "individual" - this module is named "individual" and is important for the solver to work.

Run find_BG_parameters.py

# Composite
Navigate to parent folder.
Run find_BG_parameters_composite.py

The modules included in this system are specified by the user within the .py file. The default contains the module "individual".

# Outputs
Both python scripts output the bond graph parameters as a print out and as a cellml.txt file to be copied into OpenCOR. 
The user specifies all species' initial molar amounts and ODE reactions and fluxes of each reaction manually.
