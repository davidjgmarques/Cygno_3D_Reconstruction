# CYGNO - 3D analysis

This includes: 
 - 3D
 - PMT-CAM association
 - Longitudinal diffusion studies
 - Alpha analysis
 - and other small comparison with Camera data.

I will start with the Americium dataset since we should have straight tracks thus easy to associate and or determine the 3D profile.


![3D example](examples/3D_vector_example.gif)

### Requirements

A lot, to be updated soon.


### Running the script

```
make

./alpha_3d.out data/reco_run22100_3D.root data/reco_run22100_3D_pmt.root outfile 11
```