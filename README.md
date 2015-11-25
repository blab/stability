## Analysis of stabilizing and destabilizing mutations in influenza evolution

## Pipeline

### Run nextflu pipeline and FoldX

From 'augur/' directory run nextflu pipeline to filter, align, build a tree, do ancestral state reconstruction, annotate trunk branches and list mutations. Also feeds mutations to FoldX to get the resulting ddG from the mutations listed. 

```
python run_foldxPipeline.py 4WE4 1 30
```

#### Parameters

First parameter "4WE4" specifies the protein structure to be used by foldx to calculate ddGs. Second "1" and third "30" parameters specify viruses sampled per month and years respectively. 

#### Changes to Code when changing protein structure

Before running the full pipeline, need to make sure that the nextflu pipeline is configured correctly for the HA protein structure that will be used. 

##### Outgroup
In [`H3N2_process.py`](src/H3N2_process.py) need to make the outgroup equal to the protein structure that will be used. Have added outgroups for 4WE4, 4WE5, 4WE6 and 4WE9. Just uncomment the one you're going to use and comment out the others. 

##### Amino Acid Sites

In [`tree_mutations.py`](src/tree_mutations.py) need to specify what range of amino acid sites the used protein structure covers since the structure doesn't automatically cover all sites. In calculating ddGs we have to ignore mutations that are not within the range. Have added ranges for 4WE4, 4WE5, 4WE6, 4WE9, uncomment out whichever one you're using. Change lowerRange and upperRange.

##### Pdb file repaired and formatted
Pipeline assumes that the pdb file to be used is properly repaired and formatted for foldX. Can use  [`repair_format_structure.py`](repair_format_structure.py) to repair and format the structure. From the 'augur/' directory run the following command, give the function the code for a pdb structure (ex. 4WE4)

```
python repair_format_structure.py 4WE4
```

### Generate list of mutations

From `augur/` directory run nextflu pipeline to filter, align, build a tree and do ancestral state reconstruction on this tree:

```
python src/H3N2_process.py -v 1 -y 10 --stop refine
```

Run `tree_refine` to annotate tree with trunk vs side branch etc and list mutations:

```
python src/H3N2_process.py --start mutations --stop mutations
```

### Visual data check

From `augur/`, run full pipeline:

```
python src/H3N2_process.py -v 1 -y 10`
```

Then, from `auspice/`, run:

```
jekyll serve
```

The nextflu tree will be visible at `http://localhost:4000/`.

## Notes

Change the parameters `-v 1 -y 10` to something sensible in the final analysis. This will give something easier to run during testing however.

## Citation

Please cite nextflu as:

[Neher RA, Bedford T. 2015. nextflu: real-time tracking of seasonal influenza virus evolution in humans. Bioinformatics DOI: 10.1093/bioinformatics/btv381.](http://dx.doi.org/10.1093/bioinformatics/btv381)

## License and copyright

This software is based on the nextflu pipeline.

nextflu is copyright 2014-2015 Trevor Bedford and Richard Neher.

Source code to nextflu is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). nextflu is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
