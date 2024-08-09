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

### Typical terminal output:

```c++
./alpha_3d.out debug <path_to_file>/reco_run41507_3D.root <path_to_file>/reco_run41507_3D.root out3D_507 295


CAM Reco data file openend: reco_run41507_3D.root
PMT Reco data file openend: reco_run41507_3D.root



************   Analysis  CAMERA    ************




	==> Cam run: 41507; event: 295; cluster ID: 0


--> The particle in this cluster was identified as an alpha: true

Track information: 

--> Position barycenter: x: 1098.27; y: 1147.38
--> Quadrant: 4
--> Angle: -81.9575 degrees.
--> Length (cm): 7.69849


	==> Cam run: 41507; event: 295; cluster ID: 1


--> The particle in this cluster was identified as an alpha: true

Track information: 

--> Position barycenter: x: 134.329; y: 630.402
--> Quadrant: 4
--> Angle: -166.481 degrees.
--> Length (cm): 4.78113


************   Analysis  PMT    ************




	==> PMT run: 41507; event: 295; trigger: 0; sampling: 1024


**PMT Track information: 

--> The TOT20/TOT30 ratios were: 1.15909 * 1.2 * 1.17073 * 1.2 * 
--> The TOT20 lengths were: 51 * 48 * 48 * 54 * 
--> The average number of peaks per waveform was: 1.5
--> The particle in this trigger was identified as an alpha: false


	==> PMT run: 41507; event: 295; trigger: 1; sampling: 1024


**PMT Track information: 

--> The TOT20/TOT30 ratios were: 1.05537 * 1.05161 * 1.04207 * 1.03822 * 
--> The TOT20 lengths were: 324 * 326 * 322 * 326 * 
--> The average number of peaks per waveform was: 0
--> The particle in this trigger was identified as an alpha: true
--> Ambiguous. Score: 26.6514
--> The average travelled Z (cm) is: 1.77534
--> The track is in the quadrant: 4

Launching BAT script ...
../BAT_PMTs/./runfit.out -i bat_files/input_for_bat.txt -o bat_files/output_from_bat.txt -s 0 -e 10000 -m association > bat_files/bat_system_out.txt
-> Script executed successfully.


	==> PMT run: 41507; event: 295; trigger: 2; sampling: 1024


**PMT Track information: 

--> The TOT20/TOT30 ratios were: 1.23913 * 1.38806 * 1.23256 * 1.17333 * 
--> The TOT20 lengths were: 171 * 186 * 159 * 176 * 
--> The average number of peaks per waveform was: 0.5
--> The particle in this trigger was identified as an alpha: true
--> Moving towards the GEMs with score: -53.5258
--> The average travelled Z (cm) is: 0.946483
--> The track is in the quadrant: 4

Launching BAT script ...
../BAT_PMTs/./runfit.out -i bat_files/input_for_bat.txt -o bat_files/output_from_bat.txt -s 0 -e 10000 -m association > bat_files/bat_system_out.txt
-> Script executed successfully.


************   COMBINED Analysis   ************


Associations for run 41507, picture 295
# CAM alphas: 2;  # PMT alphas: 2
==> We'll perform > 2 < associations.

---------------------------------------------------

*** Association # 0:

 -> CAM alpha # = 0, cluster = 0
 -> PMT alpha # = 0, trigger 1
 -> Distance: 71.9377

 ==> Event in run 41507, event 295, trigger 1, quadrant 4, with Alpha-PID = 1

  ** 3D Alpha track information: ** 

--> Position, X: 16.5679; Y: 20.981
--> Travelled XY: 7.69849
--> Angle XY (#phi): -81.9575
--> Travelled Z: 1.77534
--> Direction in Z: 0 at 26.6514 score
--> Angle Z (#theta): 12.9859
--> 3D alpha length (cm): 7.90054

---------------------------------------------------

---------------------------------------------------

*** Association # 1:

 -> CAM alpha # = 1, cluster = 1
 -> PMT alpha # = 1, trigger 2
 -> Distance: 176.946

 ==> Event in run 41507, event 295, trigger 2, quadrant 4, with Alpha-PID = 1

  ** 3D Alpha track information: ** 

--> Position, X: 2.52698; Y: 9.88266
--> Travelled XY: 4.78113
--> Angle XY (#phi): -166.481
--> Travelled Z: 0.946483
--> Direction in Z: -1 at 53.5258 score
--> Angle Z (#theta): 11.1976
--> 3D alpha length (cm): 4.87392

---------------------------------------------------
Deleting non-alpha events from root file...
**Finished**
```

### Output file

The root output file will contain a TTree with all the relevant information.

If necessary, you can turn on the option `bool save_everything = true;` to save all the interesting plots.

`Debug mode` automatically turns this option on.

The output file will look like this:

<div style="display: flex;">
  <div style="flex: 50%; padding: 5px;">
    <figure>
      <img src="examples/root_file.png" alt="Image 1" style="width: 100%;">
      <figcaption>(a) Root file overview.</figcaption>
    </figure>
  </div>
  <div style="flex: 50%; padding: 5px;">
    <figure>
      <img src="examples/folder_event.png" alt="Image 2" style="width: 100%;">
      <figcaption>(b) Example of saved alpha information.</figcaption>
    </figure>
  </div>
</div>