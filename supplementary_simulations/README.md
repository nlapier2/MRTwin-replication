## Replication code for simulation in the supplementary

See the scripts/ directory.

`mrtrio_controlPC.R` contains the pipeline for genotype and trait simulation, together with statistics computation for MRTwin, Brumpton and IVW.
You can easily run the script using
```
sh runmrtrio.sh
```
To change to parameters, you can modify the `params.txt` file. The explanation of each variable can be found in `mrtrio_controlPC.R`.


## Example dataset contains the following informations

### Data Description:
external.X: genotype information for external population

trio.child: genotype information for trio-child 

trio.father: genotype information for trio-father

trio.mother: genotype information for trio-mother


external.exposure: exposure trait for external population

exposure: exposure trait for the trios 

external.outcome: outcome trait for external population

outcome: outcome trait for the trios

### Folder description:
ps 1 null: no population structure, null model (no causal effect)

ps 2 null: 2 populations, null model (no causal effect)

ps 2 alter: 2 populations, Exposure has causal effect to Outcome


