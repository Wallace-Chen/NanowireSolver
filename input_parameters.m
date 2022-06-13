% this file sets the various input parameters to be used.
global l0;
global C;
global include_ex;
global draw_ex;
global is_pinned;
global pinned_fermi_level;
global dos;
global critical_height;
global pol_spont;
global pol_piezo;
global n_D;
global degree_of_polygon;
global vector_of_side_lengths;
global vector_of_V0;
global vector_of_masses;
global vector_of_eps;
global number_of_triangles;
global number_of_interfaces;
global outer_shell_side_length;
global out_folder;

% characteristic length scale, in the unit of nm
% note: if l0 is changed, you need to also change doping density in the
% file set_doping_density.m and surface donor density below.
l0 = 10;

C = hbar^2 / (2*m0*(l0*nm)^2) / eV; % scaling factor for energies, in the unit of eV

% indicator: whether we want to include exchange potentials when solving
% schrodinger equation
include_ex = false;
% indicator whether we will draw exchange potential and consider exchange
% potential when solving for fermi energy and defining depletion region. Th
% exchange potential will be taken as diag(V_ex), where V_ex is in the real
% space.
draw_ex = false;

% indicate whether we want to pin the fermi level, if not, fermi level is
% determined by the charge neutrality constraint
is_pinned = true;

% fermi level to be pinned, in the unit of eV, note we will assume the
% shell the zero potential. In this case, the shell is equi-potential, to
% apply charge neutrality, we apply boundary conditions to determine the
% empty surface states (surface holes) on the shell.
% Otherwise, we set the minimum potential as the zero potential.
pinned_fermi_level = -1.65;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define your surface donor density and the critical surface barrier height
% (highest occupied surface state)
% see: Luke Gordon et al 2010 J. Phys. D: Appl. Phys. 43 505501
% only used for triangular geometry atm
% use symbolic expression, x: energy in the unit of eV, the function should
% return the density corresponding to that energy, in the unit of 10^12
% cm^-2eV-1, if l0=10, otherwise, change accordingly.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms x;
critical_height = -1;
n0 = 100;
dos = n0 / (exp((x-critical_height)*10000 ) + 1);


% a number indicating the degree of the polygon, e.g. 3 means triangle, 6
% means hexagon, 0 means circle
% more specifically, 3 is the Ga-face triangular, -3 is the N-face triangle
degree_of_polygon = 3;

% spontaneous polarization
% You must provide if the geometry is triangular; The vector will only be
% used for the triagular case.
% a vector containing the spontaneous polarizations of each of the
% interfaces in the units of C/m2. Take care of the sign, relative to the
% [0001] direction. This vector is ordered as follows: the
% first element contains the valus of the outermost shell, and the
% subsequent elements contains of the remaining interfaces as
% one "works inwards" towards the innermost core region
pol_spont = [-0.0446, -0.029];

% piezoelectric polarization
% You must provide if the geometry is triangular; The vector will only be
% used for the triagular case.
% a vector containing the piezoelectric polarizations of each of the
% interfaces in the units of C/m2. Take care of the sign, relative to the
% [0001] direction. This vector is ordered as follows: the
% first element contains the value of the outermost shell, and the
% subsequent elements contains of the remaining interfaces as
% one "works inwards" towards the innermost core region
% this should be the name of a function to compute the longitudinal piezo
% polarization of that layer, given a theta between the natural c-axis of
% the crystal and the prime surface normal in the crystal growth diretion,
% in radians
pol_piezo = ["piezo_analytic", "piezo_zeros"];

% doping density (in units of 10^18 cm^-3, or 10^24/m^3 if l0=10,
% otherwise change unit accordingly).
n_D = 0.1;

% a vector containing the side lengths of each of the interfaces in units
% of l0 nm [see constat.m for l0]. This vector is ordered as follows: the
% first element contains the side length of the outermost shell, and the
% subsequent elements contains side lengths of the remaining interfaces as
% one "works inwards" towards the innermost core region
% vector_of_side_lengths = [4.5 3];
vector_of_side_lengths = [7 4];
% vector_of_side_lengths = [11 6];

% a vector containing the bandedge energies in units of eV [see constat.m
% for eV]. This vector is ordered as follows: the first element contains
% the bandedge energy of the outermost region, and the subsequent elements
% contain bandedge energies of the remaining regions as one "works inwards"
% towards the innermost core region
vector_of_V0 = [0.5 0.0];

% a vector containing the effective mass of each of the nanowire regions in
% units of m_0 [see constat.m for m0]. This vector is ordered as follows:
% the first element contains the effective mass of the outermost region,
% and the subsequent elements contain effective masses of the remaining
% regions as one "works inwards" towards the innermost core region
vector_of_masses = [0.2-0.12*0.3 0.2];

% a vector containing the static dielectric of each of the nanowire
% regions. This vector is ordered as follows: the first element contains
% the static dielectric of the outermost region, and the subsequent
% elements contain static dielectrics of the remaining regions as one
% "works inwards" towards the innermost core region
vector_of_eps = [9.28-0.61*0.3 9.28];

% number of triangles to be used in the finite element mesh
number_of_triangles = 1000;

% the folder under which output files are saved
out_folder = sprintf('%d_length_%f_%f_dope_%f_n0_%f_Hc_%f_ex_%d_grids_%d', degree_of_polygon, vector_of_side_lengths(1), vector_of_side_lengths(2), n_D, n0, critical_height, include_ex==true, number_of_triangles);
number_of_interfaces = length(vector_of_side_lengths);
outer_shell_side_length = vector_of_side_lengths(1);