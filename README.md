# SimulateScatteringSignal
Simulate radial distributions of X-ray Scattering Signals (monochromatic beam) of a polyatomic given its atoms names and positions 

Inputs:

  Ekev   - the X-ray photon energy in KeV (scalar)
  Name   - the names of the atoms in the polyatomic (cell array or chars)
  xyz    - the (x,y,z) positions of each atom in angstrom (nx3 matrix)

 Outputs:

  out            - A Matlab structure that contains the following:
  out.q          - the scattering vector.
  out.S0         - the total angle avg scattering signal of the molecule (S_0(q)) 
  out.pair       - the angle avg scattering signal of each atomic pair
  out.f0         - the atomic form factor of each atom (f0(q)).
  out.fafb       - the form factor of each atom pair (fa(q)*fb(q)).
  out.Rab        - the atom pair distances norm(R_a-R_b).
  out.pairlabels - the atom pair names.

for example:
```
 Ekev = 25;
 Name = {'C';'I';'I';'H';'H'};
 xyz  = [    0.142149937     -0.113392611     -0.010383833   %C
             1.804868889     -1.135903532      0.001393572   %I1
            -1.821501675     -1.141064312     -0.001561729   %I2
             0.191382262      0.584389635      0.898176095   %H1
            -0.052475117      0.636542526     -0.844064941]; %H2         

 out = SSS_polyatomic(Ekev,Name,xyz)
```
