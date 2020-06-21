# TBETA (Tight-Binding Electronic Transport Application) for graphene nanoribbon junctions

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/KrisCer/Transport_app/master?urlpath=voila%2Frender%2FWeb_app_title_page.ipynb)

<b>TBETA</b> gives you a simple graphical interface that let's you design a two-terminal Graphene nanoribbon (GNR) junction with a 60°, 120° or 180°, angle between the leads. Afterwards calculate it's transport and electronic properties based on nearest-neighbour tight-binding model. Conductance is expressed within Landauer formalism and calculated using Green's function methods.

Current version: 1.0.1

## Usage
To design the junction structure:
1)Select the desired junction angle, width of GNR leads and three additional parameters called distance, shift and chirality.
2)Atoms can be removed from the initial scattering region by entering the corresponding atom number(s) and then pressing the 'Remove atoms' button
3)Initial structure template with no atoms removed can be reset by pressing 'Reset system'
The electronic structure and transport properties are calculated:
1)Conductance and density of states (DOS) are calculated in a selected energy range by adjusting the 'min' and 'max' energy values and pressing the Calculate button
2)Local density of states, local current map and wavefunctions in the lead channels are displayed for a selected single energy value defined in terms of the tight-binding hopping integral t

## License
MIT. The terms of the license can be found in the LICENSE file.

## Contact
kristians.cernevics@epfl.ch




