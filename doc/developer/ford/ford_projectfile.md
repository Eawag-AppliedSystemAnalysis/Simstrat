src_dir: ../../../src/
output_dir: ./ford_doc
project: Simstrat
project_github: https://github.com/Eawag-AppliedSystemAnalysis/Simstrat
project_website: https://www.eawag.ch/en/department/surf/projects/simstrat/
summary: Simstrat: a one-dimensional physical lake model
author: Eawag - Department Surface Waters - Research and Management
author_description:
email: martin.schmid@eawag.ch
fpp_extensions: fpp
media_dir: ./media
docmark_alt: *
predocmark: >
predocmark_alt: <
display: public
         protected
         private
source: false
graph: true
search: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: ../../../lib/json_fortran/src/
            csv_module: ../../../lib/csv_fortran/src/
license: by-nc
extra_filetypes: sh #

Simstrat is a one-dimensional physical lake model for the simulation of stratification and mixing in deep stratified lakes. The model was originally developed by Goudsmit et al. (2002) and has been successfully applied to lakes with different physical properties. A k-ε model is used to model turbulent mixing including energy transfer of internal seiches. River or groundwater inflow can be added at specific depths or as density-dependent intrusions. The newest version of Simstrat (see below) can also simulate ice/snow covers.

@Note
You can include any notes (or bugs, warnings, or todos) like so.