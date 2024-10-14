import openmc

def lmap(func, *iterables):
    return list(map(func, *iterables))

#db = "endfb8"

f_pitch_met, f_pitch_uo2 = 2.917, 1.837

d_iso_He     = {"He4":1}

d_iso_H2O    = {"H1" :6.66500E-01,
                "H2" :1.00000E-04,
                "O16":3.33400E-01}

d_iso_Al6060 = {"Al27":9.83697E-01,
                "Si28":5.33392E-03,
                "Si29":2.71049E-04,
                "Si30":1.80715E-04,
                "Fe54":8.53577E-05,
                "Fe56":1.33637E-03,
                "Fe57":2.85411E-05,
                "Fe58":4.67491E-06,
                "Cu63":2.92655E-04,
                "Cu65":1.33481E-04,
                "Mn55":4.92979E-04,
                "Mg24":5.28454E-03,
                "Mg25":6.72047E-04,
                "Mg26":7.40078E-04,
                "Cr50":1.08450E-05,
                "Cr52":2.19000E-04,
                "Cr53":2.55789E-05,
                "Cr54":5.02111E-06,
                "Ti46":4.66785E-05,
                "Ti47":4.21088E-05,
                "Ti48":4.16858E-04,
                "Ti49":3.05980E-05,
                "Ti50":2.92823E-05,
                "Zn64":2.99670E-04,
                "Zn66":1.74085E-04,
                "Zn67":2.54901E-05,
                "Zn68":1.18127E-04,
                "Zn70":3.73019E-06}

d_iso_Umet   = {"U235":9.58998E-03,
                "U238":9.90410E-01}

d_iso_UO2    = {"U235":6.09556E-03,
                "U238":3.27238E-01,
                "O16" :6.66666E-01}


Helium = openmc.Material(name='Helium-4',temperature = 300)
Helium.set_density('g/cm3', 1.6422e-4)
Helium.add_nuclide(nuclide="He4",percent=1,percent_type='ao')

Umet = openmc.Material(name='Metallic Uranium',temperature = 300)
Umet.set_density('g/cm3', 16.67655) #true density = 18.67655
Umet.add_nuclide(nuclide="U235",percent=9.58998E-03,percent_type='ao')
Umet.add_nuclide(nuclide="U238",percent=9.90410E-01,percent_type='ao')

Water = openmc.Material(name = 'H2O',temperature= 300)
Water.set_density('g/cm3', 0.9983)
for isotope, abundance in d_iso_H2O.items():
    Water.add_nuclide(nuclide= isotope,percent= abundance)

Aluminium = openmc.Material(name = 'Alu6060',temperature= 300)
Aluminium.set_density('g/cm3', 2.702)
for isotope, abundance in d_iso_Al6060.items():
    Aluminium.add_nuclide(nuclide= isotope,percent= abundance)

UO2 = openmc.Material(name = 'Oxide Uranium',temperature= 300)
UO2.set_density('g/cm3', 10)#true density = 10.55553
for isotope, abundance in d_iso_UO2.items():
    UO2.add_nuclide(nuclide= isotope,percent= abundance)

MA = [Helium,Umet,Water,Aluminium,UO2]
materials = openmc.Materials(MA)
materials.export_to_xml()
#add density later on... 

#Geometry 

def cyl(x=0,y=0,z=0,h=None,r=None, bl="transmission") : return openmc.model.RightCircularCylinder((x,y,z), height=h, radius=r, boundary_type=bl)

rec = openmc.model.RectangularParallelepiped

def make_lc_pin(z, h, l_r, l_mat_name):
    l_s = lmap(lambda r: cyl(z=z, h=h,r=r), l_r)
    l_c = [openmc.Cell(region=-l_s[0], fill=l_mat_name[0])]
    for i in range(1, len(l_r)):
        l_c += [openmc.Cell(region=+l_s[i-1] & -l_s[i], fill=l_mat_name[i])]
    l_c += [openmc.Cell(region=+l_s[i], fill=(l_mat_name[i+1]))]
    return l_c

s_cyl_ext = cyl(z=-5, h=110, r=65, bl="vacuum")

r_uo2     = -rec(-f_pitch_uo2*11, f_pitch_uo2*11, -f_pitch_uo2*3,  f_pitch_uo2*3 , 0, 100) | \
            -rec(-f_pitch_uo2*9 , f_pitch_uo2*9 , -f_pitch_uo2*6,  f_pitch_uo2*6 , 0, 100) | \
            -rec(-f_pitch_uo2*6 , f_pitch_uo2*6 , -f_pitch_uo2*9,  f_pitch_uo2*9 , 0, 100) | \
            -rec(-f_pitch_uo2*3 , f_pitch_uo2*3 , -f_pitch_uo2*11, f_pitch_uo2*11, 0, 100)
r_met     = -rec(-f_pitch_met*7 , f_pitch_met*7 , -f_pitch_met*2,  f_pitch_met*2 , 0, 100) | \
            -rec(-f_pitch_met*6 , f_pitch_met*6 , -f_pitch_met*4,  f_pitch_met*4 , 0, 100) | \
            -rec(-f_pitch_met*4 , f_pitch_met*4 , -f_pitch_met*6,  f_pitch_met*6 , 0, 100) | \
            -rec(-f_pitch_met*2 , f_pitch_met*2 , -f_pitch_met*7,  f_pitch_met*7,  0, 100)
        
c_reg_uo2 = openmc.Cell(region=              r_uo2, fill= Water)
c_reg_met = openmc.Cell(region=-s_cyl_ext & ~r_met, fill= Water)

univ_ocean = openmc.Universe( cells=[openmc.Cell(fill=Water)])
univ_UO2   = openmc.Universe( cells=make_lc_pin(0,100, [0.526, 0.545,  0.63  ], [UO2,  Helium, Aluminium, Water]))
univ_met   = openmc.Universe( cells=make_lc_pin(0,100, [0.85,  0.8675, 0.9675], [Umet, Helium, Aluminium, Water]))

lattice_uo2 = openmc.RectLattice()
lattice_uo2.lower_left = [-11*f_pitch_uo2, -11*f_pitch_uo2]
lattice_uo2.pitch = [f_pitch_uo2, f_pitch_uo2]
lattice_uo2.outer = univ_ocean
lattice_uo2.universes = [[univ_UO2]*22]*22
c_reg_uo2.fill = lattice_uo2

lattice_met = openmc.RectLattice()
lattice_met.lower_left = [-10*f_pitch_met, -10*f_pitch_met]
lattice_met.pitch = [f_pitch_met, f_pitch_met]
lattice_met.outer = univ_ocean
f , um, uo = univ_ocean, univ_met, univ_UO2
lattice_met.universes = [[f , f , f , f , f , f , f , um, um, um, um, um, um, f , f , f , f , f , f , f ],
                            [f , f , f , f , um, um, um, um, um, um, um, um, um, um, um, f , f , f , f , f ],
                            [f , f , f , um, um, um, um, um, um, um, um, um, um, um, um, um, um, f , f , f ],
                            [f , f , um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, f , f ],
                            [f , um, um, um, f , um, um, um, um, um, um, um, um, um, um, um, um, um, f , f ],
                            [f , um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, f ],
                            [um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, f ],
                            [um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um],
                            [um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um],
                            [um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um],
                            [um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um],
                            [um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um],
                            [um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um],
                            [f , um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um],
                            [f , um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, f ],
                            [f , f , um, um, um, um, um, um, um, um, um, um, um, um, um, f , um, um, um, f ],
                            [f , f , um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, um, f , f ],
                            [f , f , f , um, um, um, um, um, um, um, um, um, um, um, um, um, um, f , f , f ],
                            [f , f , f , f , f , um, um, um, um, um, um, um, um, um, um, um, f , f , f , f ],
                            [f , f , f , f , f , f , f , um, um, um, um, um, um, f , f , f , f , f , f , f ]]
c_reg_met.fill = lattice_met
c_inter = openmc.Cell(region=r_met & ~r_uo2, fill=Water)
geometry = openmc.Geometry(openmc.Universe(cells=[c_reg_uo2, c_reg_met, c_inter]))
geometry.export_to_xml()

plot_xy = openmc.Plot()
plot_xy.filename = 'plot_xy'
plot_xy.origin = (0, 0, 50)  # Center of the plot (x, y, z)
plot_xy.width = (70, 70)  # Width of the plot in the x and y directions
plot_xy.pixels = (1000, 1000)  # Resolution of the plot
plot_xy.color_by = 'material'  # Color by material
plot_xy.basis = 'xy'  # XY plane

plot_xyc = openmc.Plot()
plot_xyc.filename = 'plot_xyc'
plot_xyc.origin = (0, 0, 50)  # Center of the plot (x, y, z)
plot_xyc.width = (70, 70)  # Width of the plot in the x and y directions
plot_xyc.pixels = (1000, 1000)  # Resolution of the plot
plot_xyc.color_by = 'cell'  # Color by material
plot_xyc.basis = 'xy'  # XY plane


plots = openmc.Plots([plot_xy, plot_xyc])
plots.export_to_xml()


settings = openmc.Settings()
settings.run_mode = 'fixed source'
settings.retina = { 'mt_numbers' : [18], 'nuclide_ids' : [92235, 92238],  \
                   'max_particles' : 100000, 'mcpl' : False}                   
settings.batches = 10
settings.inactive = 1
settings.particles = 100
settings.export_to_xml()

#show settings.xml





