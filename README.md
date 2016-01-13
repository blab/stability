## Analysis of stabilizing and destabilizing mutations in influenza evolution

## Pipeline

### Run H3N2_process to get new sequences

From 'augur/' directory run augur pipeline to filter, align, build a tree, do ancestral state reconstruction, annotate trunk branches and list sequences that have not yet had ddG calculated for them. It checks the stability table in dynamodb to see if there has been a stability value calculated. New sequences are put in new_eq_file.txt in /stability-data. 

```
python src/H3N2_proccess.py --stop mutations
```

### Calculate ddG for new sequences

Copy the /stability-data/ directory onto the cluster. Change the batch file and foldx file to be executable. Define AWS keys to access the stability table in dynamodb. 

```
chmod 755 stability_sbatch_call.sh
cd foldx_essentials/
chmod 755 foldx3b6
cd ..
export AWS_ACCESS_KEY_ID=''
export AWS_SECRET_ACCESS_KEY=''
export AWS_DEFAULT_REGION='us-west-2'
```

Then call the sbatch file

```
sbatch ./stability_sbatch_call.sh
```
Wait till its completed.

### Complete H3N2_process

After calculating a stability value for all sequences in the current run, starting from the stability step, finish the rest of the augur pipeline. You should now be able to visualize the tree from the auspice directory. Currently change in stability values from the outgroup are put in epitope values place. 

```
python src/H3N2_proccess.py --start stability
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
