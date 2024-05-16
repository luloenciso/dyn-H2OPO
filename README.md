dyn-H2PO
========
dynamics of water (H2O) in Paper-Oil insulation

This repository contains a time-dependent model to simulate the water content dynamics in cellulose-oil insulation. The model considers the diffusion of water in cellulose as if it were a single time constant (for small thicknesses). Only mineral oil is considered.

Overview
---------
Key features:

* Simulate the water content dynamics in cellulose-oil insulation over time
* Consider the temperature to which the insulation is subjected at each time step
* Simplified simulation that considers water diffusion in cellulose as a single time constant
* Only mineral oil is considered
* Use different types of cellulose (kraft paper or pressboard)
* Use different conditions of cellulose (new or aged)
* Include different qualities of dielectric oil modeled with oil acidity [mg KOH/g]
* Include different types of oil modeled with aromatic content [%]
* Use oil-cellulose water content relative saturation equilibrium curves

## Example of simulation
The next gif is generated with matplotlib with a simulation resulted of using this model. 
The modeled insulation it has the next characteristics:
- initial cellulose water content: 2.1%
- Pressboard cellulose 
- Condition cellulose: aged
- Oil Acidity: 0.05 mg KOH/g
- Oil Aromatic content: 15%

![Water dynamics in paper-oil](example_gif_simulated/water_dinamic_in_paper_oil.gif)

## License:

This repository is licensed under the GNU GPLv3.

## Contact:

If you have any questions, please contact luloenciso@gmail.com.

I hope this is helpful! Let me know if you have any other questions.
